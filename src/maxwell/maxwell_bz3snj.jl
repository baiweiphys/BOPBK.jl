"""
Description: Calculate the coefficients of bz3snj for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2026-02-03
"""


function maxwell_bz3snj(
    s::Int, n::Int, jj::Int,
    theta::Float64,
    b12snj, csnj,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)


    coef = wps[s]^2 * b12snj(s, n, jj) * (csnj(s, n, jj) / wcs[s]^2 - 2 * n / wcs[s] + n^2 / csnj(s, n, jj))
    bz3snj = -1im * phys.epsilon0 * tan(theta)^2 * coef
    return bz3snj

end