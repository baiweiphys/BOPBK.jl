"""
Description: Calculate the coefficients of by3snj for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2026-02-03
"""


function maxwell_by3snj(
    s::Int, n::Int, jj::Int,
    theta::Float64,
    b34snj, csnj,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)

    by3snj = -1.0 * phys.epsilon0 * tan(theta) * wps[s]^2 *
             (1 - n * wcs[s] / csnj(s, n, jj)) * b34snj(s, n, jj) / wcs[s]

    return by3snj

end
