"""
Description: Calculate the coefficients of bx21snl for the x-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bx21snl(
    s::Int, n::Int, l::Int,
    b3snl, wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    bx21snl = phys.epsilon0 * wps[s] .^ 2 * n * b3snl(s, n, l)

    return bx21snl
end