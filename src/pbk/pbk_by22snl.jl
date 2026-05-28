"""
Description: Calculate the coefficients of by22snl for the y-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_by22snl(
    s::Int, n::Int, l::Int,
    b6snl, wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    by22snl = -1im * phys.epsilon0 * wps[s] .^ 2 * b6snl(s, n, l)

    return by22snl

end
