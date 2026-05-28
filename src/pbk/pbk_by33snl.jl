"""
Description: Calculate the coefficients of by33snl for the y-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_by33snl(
    s::Int, n::Int, l::Int,
    theta::Float64,
    b3snl,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)


    by33snl = -1.0 * phys.epsilon0 * tan(theta) * wps[s] .^ 2 * b3snl(s, n, l) / wcs[s]

    return by33snl
end
