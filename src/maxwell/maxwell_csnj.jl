"""
Description: Calculate the coefficients for the oblique plasma 
wave model with a bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
LastEditTime: 2026-02-5
"""

const utils_path = abspath(joinpath(@__DIR__, "..", "utils"))
include(joinpath(utils_path, "func_Jpole.jl"))


function maxwell_csnj(
    s::Int, n::Int, jj::Int,
    kz::Float64,
    J_opt::Float64,
    vtsz::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    us0::AbstractVector{Float64})

    # J-pole
    (_, cj) = func_Jpole(J_opt)
    csnj = n * wcs[s] + kz * us0[s] + cj[jj] * kz * vtsz[s]

    return csnj
end
