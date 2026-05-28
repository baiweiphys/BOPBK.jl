"""
Description: Calculate the coefficients of bz1snj for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2026-02-03
"""

function maxwell_bz1snj(
    s::Int, n::Int, jj::Int,
    theta::Float64,
    b12snj, csnj,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)

    coef = wps[s]^2 * n * (1 - n * wcs[s] / csnj(s, n, jj)) * b12snj(s, n, jj) / wcs[s]
    bz1snj = -1im * phys.epsilon0 * tan(theta) * coef

    return bz1snj
end