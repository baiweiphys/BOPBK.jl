"""
Description: Calculate the coefficients of by2snj for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2026-02-03
"""

function maxwell_by2snj(
    s::Int, n::Int, jj::Int,
    b56snj, csnj,
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    by2snj = -1im * phys.epsilon0 * wps[s]^2 * b56snj(s, n, jj) / csnj(s, n, jj)

    return by2snj
end