"""
Description: Calculate the coefficients of bz1 for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2026-02-10
"""

function maxwell_bz1(
    S_bm::Int,
    Ns_bm::AbstractVector{Int},
    J::Int,
    theta::Float64,
    b12snj, csnj,
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    coef1 = sum(wps .^ 2)
    coef2 = 0.0

    @inbounds for s = 1:S_bm
        N_bm = Ns_bm[s]
        for (_, n) in enumerate(-N_bm:N_bm)
            for jj = 1:J
                coef2 = coef2 + wps[s]^2 * n^2 * b12snj(s, n, jj) / csnj(s, n, jj)
            end
        end
    end
    bz1 = -1im * phys.epsilon0 * tan(theta) * (coef1 + coef2)

    return bz1

end