"""
Description: Calculate the coefficients of b3snl for the oblique 
plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-11-20
Last Modified: 2026-02-05
"""

# local
include("pbk_bsnl.jl")

const utils_path = abspath(joinpath(@__DIR__, "..", "utils"))
include(joinpath(utils_path, "funcS_pbk.jl"))


function pbk_b3snl(
    s::Int, n::Int, l::Int,
    kz::Float64,
    kappasz::AbstractVector{Int},
    kappasx::AbstractVector{Int},
    vtsz::AbstractVector{Float64},
    sgms::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
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

    S3 = funcS_pbk(n, kappax_val, lambda_val, sgm_val, 1.0, 1.0, 2.0, 2.0, EPS0)
    S8 = funcS_pbk(n, kappax_val, lambda_val, sgm_val, 1.0, 1.0, 0.0, 1.0, EPS0)
    S38 = (kappax_val + sgm_val + 1.0) / kappax_val * S3 - 2.0 * sgm_val * lambda_val * S8

    bsnl = pbk_bsnl(s, n, l, kappasz, vtsz, wcs, kz)
    # b3snl = S38 * bsnl

    return S38 * bsnl

end
