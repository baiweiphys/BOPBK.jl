"""
Description: Calculate the coefficients of bx1snj for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2026-02-03
"""

# include("./constants.jl")
# using .PhysicalConstants

function  maxwell_bx1snj(
    s::Int, n::Int, jj::Int, 
    b12snj, csnj, 
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    bx1snj = -1im*phys.epsilon0*wps[s]^2*n^2*b12snj(s,n,jj)/csnj(s,n,jj)
    return bx1snj
end
