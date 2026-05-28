"""
Description: Calculate the coefficients of bz2snj for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2026-12-03
Last Modified: 2026-02-03
"""


function maxwell_bz2snj(
    s::Int, n::Int, jj::Int,
    theta::Float64,
    b34snj, csnj,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)

    coef = wps[s]^2 * (1 - n * wcs[s] / csnj(s, n, jj)) * b34snj(s, n, jj) / wcs[s]
    bz2snj = phys.epsilon0 * tan(theta) * coef
    return bz2snj
end