"""
Description: Calculate the coefficients of bx12snl for the x-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bx12snl(
    s::Int, n::Int, l::Int,
    b2snl, wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    bx12snl = -1im * phys.epsilon0 * wps[s] .^ 2 * n^2 * b2snl(s, n, l)

    return bx12snl
end
