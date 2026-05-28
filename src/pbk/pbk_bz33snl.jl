"""
Description: Calculate the coefficients of bz33snl for the z-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bz33snl(
    s::Int, n::Int, l::Int,
    theta::Float64,
    csn, b1snl, b2snl,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)

    bz33snl = -1im * phys.epsilon0 * tan(theta)^2 * wps[s]^2 *
              (2 * b1snl(s, n, l) * (csn(s, n) - n * wcs[s]) + b2snl(s, n, l)) / wcs[s]^2

    return bz33snl
end