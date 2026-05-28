"""
Description: Calculate the coefficients of by21snl for the y-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_by21snl(
    s::Int, n::Int, l::Int,
    b5snl,
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    by21snl = -1im * phys.epsilon0 * wps[s] .^ 2 * b5snl(s, n, l)

    return by21snl
end