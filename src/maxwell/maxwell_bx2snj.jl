"""
Description: Calculate the coefficients of bx2snj for the oblique 
plasma wave model with a bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2026-02-03
"""


function maxwell_bx2snj(
    s::Int, n::Int, jj::Int,
    b34snj, csnj,
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    coef = wps[s]^2 * n * b34snj(s, n, jj) / csnj(s, n, jj)
    bx2snj = phys.epsilon0 * coef

    return bx2snj

end