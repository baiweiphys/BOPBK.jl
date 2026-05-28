"""
Description: Calculate the coefficients of by32snl for the y-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_by32snl(
    s::Int, n::Int, l::Int,
    theta::Float64,
    csn, b4snl,
    wps::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    phys::PlasmaConfig)

    by32snl = -1.0 * phys.epsilon0 * tan(theta) * wps[s] .^ 2 *
              b4snl(s, n, l) * (csn(s, n) - n * wcs[s]) / wcs[s]

    return by32snl

end