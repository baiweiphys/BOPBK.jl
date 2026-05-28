"""
Description: Calculate the coefficients of bz2 for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2025-09-03
Last Modified: 2026-02-10
"""

function maxwell_bz2(
    S_bm::Int,
    Ns_bm::AbstractVector{Int},
    J::Int,
    theta::Float64,
    b34snj, csnj,
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)


    coef = 0.0
    @inbounds for s = 1:S_bm
        N_bm = Ns_bm[s]
        for (_, n) in enumerate(-N_bm:N_bm)
            for jj = 1:J
                coef = coef + wps[s]^2 * n * b34snj(s, n, jj) / csnj(s, n, jj)
            end
        end
    end
    bz2 = phys.epsilon0 * tan(theta) * coef

    return bz2
end
