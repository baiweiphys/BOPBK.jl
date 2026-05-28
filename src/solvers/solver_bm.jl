"""
Description: Optimized solver for BM plasma. (In-place).
Uses optimized M_maxwell! functions to assemble the global matrix for the BM plasma case.
Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
Created: 2025-12-16
Changelog:
2026-02-08: Optimized solver_bm. 
"""

const maxwell_path = abspath(joinpath(@__DIR__, "..", "maxwell"))
include(joinpath(maxwell_path, "maxwell_csnj.jl"))
include(joinpath(maxwell_path, "maxwell_b12snj.jl"))
include(joinpath(maxwell_path, "maxwell_b34snj.jl"))
include(joinpath(maxwell_path, "maxwell_b56snj.jl"))
#
include(joinpath(maxwell_path, "maxwell_bx1.jl"))
include(joinpath(maxwell_path, "maxwell_bx1snj.jl"))
include(joinpath(maxwell_path, "maxwell_bx2.jl"))
include(joinpath(maxwell_path, "maxwell_bx2snj.jl"))
include(joinpath(maxwell_path, "maxwell_bx3.jl"))
include(joinpath(maxwell_path, "maxwell_bx3snj.jl"))

include(joinpath(maxwell_path, "maxwell_by1.jl"))
include(joinpath(maxwell_path, "maxwell_by1snj.jl"))
include(joinpath(maxwell_path, "maxwell_by2.jl"))
include(joinpath(maxwell_path, "maxwell_by2snj.jl"))
include(joinpath(maxwell_path, "maxwell_by3.jl"))
include(joinpath(maxwell_path, "maxwell_by3snj.jl"))
#
include(joinpath(maxwell_path, "maxwell_bz1.jl"))
include(joinpath(maxwell_path, "maxwell_bz1snj.jl"))
include(joinpath(maxwell_path, "maxwell_bz2.jl"))
include(joinpath(maxwell_path, "maxwell_bz2snj.jl"))
include(joinpath(maxwell_path, "maxwell_bz3.jl"))
include(joinpath(maxwell_path, "maxwell_bz3snj.jl"))

# matrices
const matrices_path = abspath(joinpath(@__DIR__, "..", "matrices"))
include(joinpath(matrices_path, "M_maxwell.jl"))


function solver_bm!(
    cache,
    kx::Float64, 
    kz::Float64,
    theta::Float64, 
    PlasmaParameters,
    select_eigsolver::Int, 
    EPS0::Float64,
    phys::PlasmaConfig)

       
    # Parameters
    (; S, Ns, J_opt) = PlasmaParameters.meta
    # (; is_pbk, is_bm) = PlasmaParameters.flags
    (; Tsz, Tsx, vtsz, vtsx, us0) = PlasmaParameters.thermo
    (; wps, wcs, rhocs, lambdaDs, ms, ns0) = PlasmaParameters.physics
    (; kappasz, kappasx, sgms) = PlasmaParameters.dist_params
    J = floor(Int, J_opt)
    lambdas = 0.5 * kx^2 * rhocs .^ 2  # for argument of the modified bessel function
    #
    Ns_bm = Ns

    mask = abs.(lambdas) .< 1.0e-40
    lambdas[mask] .= 1.0e-40
    # lambdas[abs(lambdas)<1.0e-40] = 1.0e-40  # 2024.0916, to avoid singular when k_perp=0

    # BM functions
    @inline csnj(s,n,jj) = maxwell_csnj(s,n,jj,kz,J_opt,vtsz,wcs,us0)
    # bsnj = @(s,n,jj) maxwell_bsnj(s,n,jj,kz,J_opt,vts_parallel,Ts_parallel,Ts_perp,wcs)
    @inline b12snj(s,n,jj) = maxwell_b12snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0)
    @inline b34snj(s,n,jj) = maxwell_b34snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0)
    @inline b56snj(s,n,jj) = maxwell_b56snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0)

    # for x-component (bi-Maxwell) 
    bx1_bm = maxwell_bx1(S,Ns,J,b12snj,csnj,wps,phys)
    @inline bx1snj(s,n,jj) = maxwell_bx1snj(s,n,jj,b12snj,csnj,wps,phys)
    bx2_bm = maxwell_bx2(S,Ns,J,b34snj,csnj,wps,phys)
    @inline bx2snj(s,n,jj) = maxwell_bx2snj(s,n,jj,b34snj,csnj,wps,phys)
    bx3_bm = maxwell_bx3(S,Ns,J,theta,b12snj,csnj,wps,phys)
    @inline bx3snj(s,n,jj) = maxwell_bx3snj(s,n,jj,theta,b12snj,csnj,wps,wcs,phys)

    # for y-component (bi-Maxwell) 
    by1_bm = maxwell_by1(S,Ns,J,b34snj,csnj,wps,phys)
    @inline by1snj(s,n,jj) = maxwell_by1snj(s,n,jj,b34snj,csnj,wps,phys)
    by2_bm = maxwell_by2(S,Ns,J,b56snj,csnj,wps,phys)
    @inline by2snj(s,n,jj) = maxwell_by2snj(s,n,jj,b56snj,csnj,wps,phys)
    by3_bm = maxwell_by3(S,Ns,J,theta,b34snj,csnj,wps,phys)
    @inline by3snj(s,n,jj) = maxwell_by3snj(s,n,jj,theta,b34snj,csnj,wps,wcs,phys)

    # for z-component (bi-Maxwell) 
    bz1_bm = maxwell_bz1(S,Ns,J,theta,b12snj,csnj,wps,phys)
    @inline bz1snj(s,n,jj) = maxwell_bz1snj(s,n,jj,theta,b12snj,csnj,wps,wcs,phys)
    bz2_bm = maxwell_bz2(S,Ns,J,theta,b34snj,csnj,wps,phys)
    @inline bz2snj(s,n,jj) = maxwell_bz2snj(s,n,jj,theta,b34snj,csnj,wps,wcs,phys)
    bz3_bm = maxwell_bz3(S,Ns,J,theta,b12snj,csnj,wps,phys)
    @inline bz3snj(s,n,jj) = maxwell_bz3snj(s,n,jj,theta,b12snj,csnj,wps,wcs,phys)


    # MEMORY ALLOCATION FOR THE GLOBAL MATRIX
    if !cache.is_initialized
        # 1. Calc BM size
        len_M_bm = J * (2 * sum(Ns_bm) + length(Ns_bm)) + 1
        # 2. Total structure
        total_rows = 3 * len_M_bm + 9
        total_cols = total_rows

        cache.M = zeros(ComplexF64, total_rows, total_cols)
        cache.is_initialized = true
    end

    M = cache.M
    fill!(M, 0.0 + 0.0im)  # Reset the matrix to zero before filling it with new values.


    # ASSEMBLE THE GLOBAL MATRIX
    """
    # ExyzNo = 8 for Jx1_bm
    # ExyzNo = 7 for Jy1_bm
    # ExyzNo = 6 for Jz1_bm
    # ExyzNo = 5 for Ex
    # ExyzNo = 4 for Ey
    # ExyzNo = 3 for Ez
    # BxyzNo = 2 for Bx
    # BxyzNo = 1 for By
    # BxyzNo = 0 for Bz
    """
    current_row = 1  # To keep track of the current row index for filling the matrix


    ## Step 1.1 Mx BM (MatrixNo=1)
    by20 = 1im*phys.epsilon0*sum(wps.^2)
    h_x_bm = M_maxwell!(M, current_row, S, Ns, J, csnj, bx1snj, bx2snj, bx3snj, 1, 5, 4, 3)  # Matrix No.1
    current_row += h_x_bm

    # Step 1.2 My BM (MatrixNo=2)
    h_y_bm = M_maxwell!(M, current_row, S, Ns, J, csnj, by1snj, by2snj, by3snj, 2, 5, 4, 3)  # Matrix No.2
    current_row += h_y_bm

    # Step 1.3 Mz BM (MatrixNo=3)
    h_z_bm = M_maxwell!(M, current_row, S, Ns, J, csnj, bz1snj, bz2snj, bz3snj, 3, 5, 4, 3)  # Matrix No.3
    current_row += h_z_bm
    

    ## Step 2: The perturbed currents of first term of BM
    # Incorporation the first-term of BM (Jx1_bm, Jy1_bm, Jz1_bm).
    # BM perturbations (Jx1)
    M[end-8,end-5] = M[end-8,end-5] + bx1_bm  # the first term of dEx
    M[end-8,end-4] = M[end-8,end-4] + bx2_bm  # the second term of dEy
    M[end-8,end-3] = M[end-8,end-3] + bx3_bm  # the third term of dEz

    # BM perturbations (Jy1)
    M[end-7,end-5] = M[end-7,end-5] + by1_bm  # the first term of dEx
    M[end-7,end-4] = M[end-7,end-4] + by2_bm  # the second term of dEy
    M[end-7,end-3] = M[end-7,end-3] + by3_bm  # the third term of dEz

    # BM perturbations (Jz1)
    M[end-6,end-5] = M[end-6,end-5] + bz1_bm  # the first term of dEx
    M[end-6,end-4] = M[end-6,end-4] + bz2_bm  # the second term of dEy
    M[end-6,end-3] = M[end-6,end-3] + bz3_bm  # the third term of dEz

    
    ### Step 3: Maxwell's equation
    # Maxwell's equations for Ex
    idx_Jx_bm = getIndexOfBlkMatrix_bm(Ns,J,1)
    M[end-5,end-1] = M[end-5,end-1] + phys.c2*kz  # for the first item of Ex
    M[end-5,end-8] = M[end-5,end-8] - 1im/phys.epsilon0  # for dJx1
    M[end-5,idx_Jx_bm[1:end-1]] = M[end-5,idx_Jx_bm[1:end-1]] .- 1im/phys.epsilon0  # for dJx2

    # Maxwell's equations for Ey
    idx_Jy_bm = getIndexOfBlkMatrix_bm(Ns,J,2);
    M[end-4,end-2] = M[end-4,end-2] - phys.c2*kz  # for the first item of Ey
    M[end-4,end] = M[end-4,end] + phys.c2*kx  # for the second item of Ey
    M[end-4,end-7] = M[end-4,end-7] - 1im/phys.epsilon0  # for dJy1
    M[end-4,idx_Jy_bm[1:end-1]] = M[end-4,idx_Jy_bm[1:end-1]] .- 1im/phys.epsilon0  # for dJy2

    # Maxwell's equations Ez
    idx_Jz_bm = getIndexOfBlkMatrix_bm(Ns,J,3)
    M[end-3,end-1] = M[end-3,end-1] - phys.c2*kx  # for the first item of Ez
    M[end-3,end-6] = M[end-3,end-6] - 1im/phys.epsilon0 # for dJz1
    M[end-3,idx_Jz_bm[1:end-1]] = M[end-3,idx_Jz_bm[1:end-1]] .- 1im/phys.epsilon0 # for dJz2

    # Maxwell's equations for Bx, By and Bz
    M[end-2,end-4] = M[end-2,end-4] - kz  # for Bx
    M[end-1,end-5] = M[end-1,end-5] + kz  # for the first item of By
    M[end-1,end-3] = M[end-1,end-3] - kx  # for the second item of By
    M[end,end-4] = M[end,end-4] + kx      # for Bz

    ##Solver##
    eigenvalues, eigenvectors = eigensolver(M, select_eigsolver)
    # eigenvalues, eigenvectors = eigensolver(M, select_eigsolver=1, nev=10, sort_by=:real_desc)
    # eigenvalues, eigenvectors = eigensolver(M, select_eigsolver=1, nev=10)

    return (eigenvalues, eigenvectors)

end