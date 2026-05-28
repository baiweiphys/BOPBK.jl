"""
Description: Calculate the coefficients of bx11snl for the x-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bx11snl(
    s::Int, n::Int, l::Int,
    b1snl, wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    bx11snl = -1im * phys.epsilon0 * wps[s] .^ 2 * n^2 * b1snl(s, n, l)

    return bx11snl
end