"""
Description: Core implementation of oblique plasma wave dispersion relation solver.
Features a composite velocity distribution function combining:
1. Product-bi-Kappa distribution (anisotropic kappa)
2. Kappa-Maxwellian hybrid distribution  
3. Bi-Maxwellian distribution (temperature anisotropy)
This hybrid approach enables modeling of non-Maxwellian space plasmas.
Ref. Astfalk, P., and F. Jenko (2017), LEOPARD: A grid-based dispersion 
relation solver for arbitrary gyrotropic distributions, 
J. Geophys. Res. Space Physics, 122, 89–101, 
% doi:10.1002/2016JA023522.
Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
Date: 2023-08-1
Last Modified: 2026-02-05
"""


using CSV
using DataFrames
# using DelimitedFiles   # This provides writedlm and readdlm
using ProgressMeter
# using HDF5  # for .h5 files
# using JLD2  # for .jld2 files
# using MAT   # for .mat files

const PROJECT_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "tools", "create_data_folder.jl"))
include(joinpath(PROJECT_ROOT, "tools", "save_results.jl"))
include(joinpath(PROJECT_ROOT, "tools", "read_bopbk.jl"))
include(joinpath(PROJECT_ROOT, "tools", "display_runtime.jl"))
include(joinpath(PROJECT_ROOT, "src", "BOPBK.jl"))
using .BOPBK: select_solver, PlasmaParameters

if !@isdefined(phys)
    const phys = BOPBK.PhysicalConstants.SI
end


##################
# Input Parameters
##################
# J-pole
"""
# J=2.1:  Huba2009, 2-pole
# J=2.2:  Martin1979, 2-pole, J=2, I=3
# J=3:    Martin1980, 3-pole, J=3, I=3
# J=4.1:  Martin1980, 4-pole, J=4, I=5 
# J=4.2:  new calculation, J=4, I=5
# J=8.1:  Ronnmark1982, 8-pole for Z function
# J=8.2:  J=8, I=8
# J=8.3:  J=8, I=10
# J=8.4:  optimized J=8 pole from Xie 2024
# J=10.1  J=10 (I=12,K=8) pole from Xie 2024
# J=10.2  J=10 (I=14,K=6) pole from Xie 2024
# J=12.1: J=12; I=16; 2014;
# J=12.2: J=12; I=16; 2018
# J=12.3: J=12; I=12;
# J=16.1: J=16; I=18; 
# J=16.2: J=16 (I=16,K=20) pole from Xie 2024 
# J=16.3: J=16 (J=24,I=8) pole from Xie 2024
# J=20:   J=20 (I=23,K=17) pole from Xie 2024
# J=24.1: (J=24,I=24) pole from BO code
# J=24.2: J=24 (I=24,K=24) pole from Xie 2024
"""
J_opt = 8.1

# select_eigsolver=0: 'eigen()'; select_eigsolver=1: sparse 'eigs()'
# For matrices with size n ≤ 8,000:
# directly calling eigen(Matrix(M)) is the fastest and most accurate method
# memory usage is completely fine on any modern computer (2025)
select_eigsolver = 0
#
EPS0 = 1.0e-2  ## for integrals S1,S2,...,S6 when \lambda_s \rightarrow 0
B0 = 1.0e-1  ## background magnetic field in z direction 

# angle theta between k and B0
nth = 60;
theta_deg = range(1.0e-3, 90.0-1.0e-2, nth) # in degree
theta_rad = theta_deg * pi / 180.0  # in radian

#
nk = 40
kk0 = range(1.0e-3, 0.8, nk)

# load bopbk.in
df, par_data = read_bopbk("./bopbk.in")

# Get plasma parameters
plasmaParams = PlasmaParameters(B0, par_data, phys, J_opt)
(; S, Ns, J_opt) = plasmaParams.meta
(; is_pbk, is_bm) = plasmaParams.flags
(; Tsz, Tsx, vtsz, vtsx, us0) = plasmaParams.thermo
(; wps, wcs, rhocs, lambdaDs, ms, ns0) = plasmaParams.physics
(; kappasz, kappasx, sgms) = plasmaParams.dist_params

mi = ms[2]
ni = ns0[2]
vA_p = B0/sqrt(phys.mu0*mi*ni) # Alfven speed of proton
di = vA_p / wcs[2]   # ion inertial length di=vA/Omega_i

# kDs = 1.0./lambdaDs
# kn = sqrt(kDs[1]^2 + kDs[2]^2)
# wn = sqrt(wps[1]^2 + wps[2]^2)

# Thermal speed of the sth component is vs = (Tsparallel/ms)^0.5
# vs = @. (phys.kB * Tsz / ms)^0.5

# gyroradius of the ion is ai = vi/Omegai*(Tiperp/Tiparallel)^0.5
# ai = vs[3]/wcs[3]*sqrt(Tsx[3]/Tsz[3])

###############################
# Precompute k and trig factors
###############################
k_scale = 1 / di
kvals = kk0 .* k_scale
sθ = sin.(theta_rad)
cθ = cos.(theta_rad)
sθ1 = (sθ isa AbstractArray) ? sθ[1] : sθ
cθ1 = (cθ isa AbstractArray) ? cθ[1] : cθ

#
solver = select_solver(plasmaParams, select_eigsolver, EPS0, phys)
###
p = Progress(nk, desc="Simulation Progress: ", dt=0.5, barglyphs=BarGlyphs("[=> ]"))
start_time = time()
# -------------------------------
# Pre-determine the number of modes (nw) by running the solver once
# -------------------------------
kval_test = kvals[1]         # Wavevector magnitude for the first k-point
kx_test = kval_test * sθ1    # x-component of the wavevector
kz_test = kval_test * cθ1    # z-component of the wavevector

# Call the solver once to find out how many eigenvalues/modes it returns
nw = length(solver(kx_test, kz_test, theta_rad[1])[1])
# Pre-allocate the matrix that will store all frequencies for all k-points
#   - size: (nk × nth × nw)  →  one row per k-point, one column per theta, one column per mode
ww = zeros(ComplexF64, nk, nth, nw)    # Pre-allocation (avoids scope issues)

# -------------------------------
# Main loop over all k-points
# -------------------------------
for ik = 1:nk
    for ith = 1:nth
        # Current wavevector components
        kx = kvals[ik] * sθ[ith]
        kz = kvals[ik] * cθ[ith]
        # Solve the eigenvalue problem for this (kx, kz)
        w, _ = solver(kx, kz, theta_rad[ith])   # w contains the complex frequencies
        # Store the frequencies for this k-point in the corresponding row
        ww[ik,ith, :] .= w
    end
    # Progress indicator
    next!(p)
end
runtime = time() - start_time
display_runtime(runtime)
println("ww type: ", typeof(ww))
println("ww size: ", size(ww))
# wr = real(ww)
# wi = imag(ww)

# -------------------------------
# Save the results
# -------------------------------
prefix::String = "bopbk"
dat_dir = create_data_folder(pwd(), "output")
save_results(prefix, dat_dir; formats=[:jld2],
    theta_deg, kk0, kvals, B0, Tsz, Tsx, us0, 
    ww, wps, wcs, rhocs, lambdaDs, ms, ns0, Ns, runtime)


# save_results(prefix, dat_dir; formats=[:h5, :mat, :jld2],
#     theta_deg, kk0, kk, B0, Tsz, Tsx, us0,
#     ww, wps, wcs, rhocs, lambdaDs, ms, ns0, Ns, runtime)

# -------------------------------
# Read the data back
# data_dict = h5open(h5_file, "r") do file
#     # Read all keys into a dictionary or assign to local variables
#     # This reads the datasets into memory
#     ww_loaded = read(file, "ww")
#     theta_loaded = read(file, "theta")
#     runtime_loaded = read(file, "runtime")

#     # Return as a tuple or dict
#     return (ww_loaded, theta_loaded, runtime_loaded)
# end

# (ww, theta, runtime) = data_dict