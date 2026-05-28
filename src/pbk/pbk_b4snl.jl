"""
Description: Calculate the coefficients of b4snl for the oblique 
plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-11-20
Last Modified: 2026-02-06
"""


# local
include("pbk_bsl.jl")

const utils_path = abspath(joinpath(@__DIR__, "..", "utils"))
include(joinpath(utils_path, "funcS_pbk.jl"))



function pbk_b4snl(
    s::Int, n::Int, l::Int,
    kz::Float64,
    kappasz::AbstractVector{Int},
    kappasx::AbstractVector{Int},
    vtsz::AbstractVector{Float64},
    vtsx::AbstractVector{Float64},
    sgms::AbstractVector{Float64},
    lambdas::AbstractVector{Float64},
    EPS0::Float64)

    # FuncS_pbk = (Jnum,dJnum,num,den) -> funcS_pbk(n,kappax_val,lambda_val,sgm_val,Jnum,dJnum,num,den,EPS0)
    # S1_pbk = FuncS_pbk(2,0,1,2)
    # S2_pbk = FuncS_pbk(2,0,1,1)
    # S3_pbk = FuncS_pbk(1,1,2,2)
    # S4_pbk = FuncS_pbk(1,1,2,1)
    # S5_pbk = FuncS_pbk(0,2,3,2)
    # S6_pbk = FuncS_pbk(0,2,3,1)
    # S7_pbk = FuncS_pbk(2,0,-1,1)
    # S8_pbk = FuncS_pbk(1,1,0,1)
    # S9_pbk = FuncS_pbk(0,2,1,1)

    kappax_val = kappasx[s]
    lambda_val = lambdas[s]
    sgm_val = sgms[s]

    S4 = funcS_pbk(n, kappax_val, lambda_val, sgm_val, 1.0, 1.0, 2.0, 1.0, EPS0)
    bsl = pbk_bsl(s, l, kappasz, vtsz, vtsx, kz)
    # b4snl = S4 * bsl

    return S4 * bsl
end
