"""
Description: Calculate the coefficients of by2 for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2026-02-10
"""

function maxwell_by2(
    S_bm::Int,
    Ns_bm::AbstractVector{Int},
    J::Int,
    b56snj, csnj,
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)


    coef1 = sum(wps .^ 2)
    coef2 = 0.0

    @inbounds for s = 1:S_bm
        N_bm = Ns_bm[s]
        for (_, n) in enumerate(-N_bm:N_bm)
            for jj = 1:J
                coef2 = coef2 + wps[s]^2 * b56snj(s, n, jj) / csnj(s, n, jj)
            end
        end
    end
    by2 = 1im * phys.epsilon0 * (coef1 + coef2)

    return by2
end