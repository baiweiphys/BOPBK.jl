"""
Description: Calculate the coefficients of bx22snl for the x-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bx22snl(
    s::Int, n::Int, l::Int,
    b4snl, wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    bx22snl = phys.epsilon0 * wps[s] .^ 2 * n * b4snl(s, n, l)

    return bx22snl
end
