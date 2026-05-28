"""
Description: Calculate the coefficients of bx32snl for the x-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bx32snl(
    s::Int, n::Int, l::Int,
    theta::Float64,
    csn, b2snl,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)


    bx32snl = -1im * phys.epsilon0 * tan(theta) * wps[s] .^ 2 * n * 
            b2snl(s, n, l) * (csn(s, n) - n * wcs[s]) / wcs[s]

    return bx32snl
end
