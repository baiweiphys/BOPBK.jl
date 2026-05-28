"""
Description: Calculate the coefficients of by31snl for the y-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_by31snl(
    s::Int, n::Int, l::Int,
    theta::Float64,
    csn,
    b3snl, b4snl,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)

    by31snl = -1.0 * phys.epsilon0 * tan(theta) * wps[s] .^ 2 *
              (b3snl(s, n, l) * (csn(s, n) - n * wcs[s]) + b4snl(s, n, l)) / wcs[s]

    return by31snl
end