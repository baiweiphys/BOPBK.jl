"""
Description: Calculate the coefficients of by1snj for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2026-02-03
"""

function maxwell_by1snj(
    s::Int, n::Int, j::Int,
    b34snj, csnj,
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    by1snj = -1.0 * phys.epsilon0 * wps[s]^2 * n * b34snj(s, n, j) / csnj(s, n, j)

    return by1snj
end