"""
Description: Calculate the coefficients of bsnl for the oblique 
plasma wave model with a bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2025-12-06
"""

include("../utils/func_Jpole.jl")

function maxwell_bsnj(
    s::Int, n::Int, jj::Int, 
    kz::Float64, 
    J_opt::Float64,
    vtsz::AbstractVector{Float64},
    Tsz::AbstractVector{Float64}, 
    Tsx::AbstractVector{Float64}, 
    wcs::AbstractVector{Float64})

    # J-pole
    (bj,cj) = func_Jpole(J_opt)
    
    bsnj = bj[jj]*cj[jj]*kz*vtsz[s]*Tsx[s]/Tsz[s] + bj[jj]*n*wcs[s]

    return bsnj

end