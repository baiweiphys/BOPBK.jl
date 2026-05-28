"""
Description: Calculate the coefficients of bz34snl for the z-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bz34snl(
    s::Int, n::Int, l::Int,
    theta::Float64,
    b1snl,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)

    bz34snl = -1im * phys.epsilon0 * tan(theta)^2 * wps[s]^2 * b1snl(s, n, l) / wcs[s]^2

    return bz34snl
end
