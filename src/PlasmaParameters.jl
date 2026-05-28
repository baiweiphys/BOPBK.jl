"""
Description: To obtain the plasma parameters.
Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
Date: 2022-01-05
Last Modified: 2026-02-03
"""

using .PhysicalConstants: PlasmaConfig

function PlasmaParameters(B0::Float64, par_data::Matrix{Float64}, phys::PlasmaConfig, J_opt=nothing)
    S = size(par_data, 1)

    # 1. Alias columns with descriptive names (using @views to ensure zero-copy)
    @views begin
        qs_raw = par_data[:, 1];
        ms_raw = par_data[:, 2];
        ns0 = par_data[:, 3];
        Tz_ev  = par_data[:, 4];  
        Tx_ev  = par_data[:, 5];  
        u_drft = par_data[:, 6];
        sgms   = par_data[:, 7];  
        kappasz = par_data[:, 8];
        kappasx = par_data[:, 9];
        kappasz_threshold = par_data[:, 10]
        Ns = par_data[:, 11];
    end

    # 2. Base physical unit conversions (Vectorized)
    qs = qs_raw .* phys.qe
    ms = ms_raw .* phys.me
    Tsz = Tz_ev .* (phys.qe / phys.kB)
    Tsx = Tx_ev .* (phys.qe / phys.kB)
    us0 = u_drft .* sqrt(phys.c2)
    # 3. Calculate thermal velocities
    vtsz = similar(Tsz)
    vtsx = similar(Tsx)
    for s in 1:S
        if kappasz[s] < kappasz_threshold[s]
            # Thermal velocity for loss-cone PBK distribution
            vtsz[s] = sqrt((2.0 - 1.0 / kappasz[s]) * phys.kB * Tsz[s] / ms[s])
            vtsx[s] = sqrt((2.0 - 2.0 / kappasx[s]) * phys.kB * Tsx[s] / ms[s] / (1.0 + sgms[s]))
        else
            # Thermal velocity for loss-cone bi-Maxwellian distribution
            vtsz[s] = sqrt(2.0 * phys.kB * Tsz[s] / ms[s])
            vtsx[s] = sqrt(2.0 * phys.kB * Tsx[s] / ms[s] / (1.0 + sgms[s]))
        end

        # (M. S. dos Santos, 2016 PoP)
        # for loss-cone PBK and bi-Maxwellian distribution
        # vtsz[s] = sqrt(2.0*phys.kB*Tsz[s]/ms[s])
        # vtsx[s] = sqrt(2.0*phys.kB*Tsx[s]/ms[s]/(1.0+sgms[s]))
    end

    # 4. Derived physical quantities
    wps = sqrt.(ns0 .* qs .^ 2 ./ (ms .* phys.epsilon0)) # plasma frequency
    wcs = B0 .* qs ./ ms # cyclotron frequency
    rhocs = vtsx ./ abs.(wcs) # cyclotron radius
    lambdaDs = sqrt.(phys.epsilon0 .* phys.kB .* Tsz ./ (ns0 .* qs .^ 2)) # Debye length
    # kDs = 1.0./lambdaDs
    betasz = 2 * phys.mu0 * phys.kB .* ns0 .* Tsz ./ B0^2  # Parallel Beta
    betasx = 2 * phys.mu0 * phys.kB .* ns0 .* Tsx ./ B0^2  # Perpendicular Beta 
    vA = B0 / sqrt(phys.mu0 * sum(ms .* ns0))  # Alfvén Velocity, where ρ = sum(m_s * n_s).
    cS = sqrt(2 * minimum(phys.kB .* Tsz) / maximum(ms)) # Characteristic Sound Speed 


    # 5. Return a Nested NamedTuple
    return (
        meta = (
            S = Int(S),
            Ns = Int.(Ns), 
            J_opt = J_opt
        ),
        
        flags = (
            # index of PBK plasmas for s-th species
            is_pbk = kappasz .< kappasz_threshold, 
            # index of bi-Maxwellian plasmas for s-th species
            is_bm  = kappasz .>= kappasz_threshold
        ),
        
        thermo = (
            Tsz = Tsz, 
            Tsx = Tsx, 
            vtsz = vtsz, 
            vtsx = vtsx, 
            us0 = us0
        ),
        
        physics = (
            wps = wps, 
            wcs = wcs, 
            rhocs = rhocs, 
            lambdaDs = lambdaDs,
            ms = ms,
            ns0 = ns0,
            betasz =  betasz,
            betasx = betasx,
            vA = vA,
            cS = cS    
        ),
        
        dist_params=(
            kappasz = Int.(kappasz),
            kappasx = Int.(kappasx),
            sgms = sgms,
            kappasz_threshold = kappasz_threshold
        ),
    )
end

# End of PlasmaParameters.jl