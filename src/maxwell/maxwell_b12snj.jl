"""
Description: Calculate the coefficients of b12snj for the oblique 
plasma wave model with a loss-cone bi-Maxwellian distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-03
Last Modified: 2025-12-11
"""

using SpecialFunctions: besselix

# local
include("maxwell_bsnj.jl")

const utils_path = abspath(joinpath(@__DIR__, "..", "utils"))
include(joinpath(utils_path, "safe_besselix.jl"))
include(joinpath(utils_path, "func_Jpole.jl"))
include(joinpath(utils_path, "funcS_maxwell.jl"))


function maxwell_b12snj(
    s::Int, n::Int, jj::Int, 
    kz::Float64, 
    J_opt::Float64, 
    vtsz::AbstractVector{Float64},
    Tsz::AbstractVector{Float64},
    Tsx::AbstractVector{Float64},
    sgms::AbstractVector{Float64}, 
    wcs::AbstractVector{Float64}, 
    lambdas::AbstractVector{Float64},
    EPS0::Float64)

    # FuncS_max = (Jnum,dJnum,num) -> funcS_max(n,lambda_val,sgm_val,Jnum,dJnum,num,EPS0)
    # S1_max = S2_max = FuncS_max(2,0,1)
    # S3_max = S4_max = FuncS_max(1,1,2)
    # S5_max = S6_max = FuncS_max(0,2,3)
    # S7_max = FuncS_max(2,0,-1)
    # S8_max = FuncS_max(1,1,0)
    # S9_max = FuncS_max(0,2,1)

    lambda_val = lambdas[s]
    sgm_val = sgms[s]
 
    # J-pole
    (bj,cj) = func_Jpole(J_opt)
    bsnj = maxwell_bsnj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,wcs)

    # Scaled Mod. Bessel Functions of the First Kind for n.
    # e.g., e^{-|z|} \cdot I_\nu(z)
    # scale = 1.0
    In = besselix_safe(n, lambda_val)

    if sgm_val==0
        b12snj = bsnj*In/lambda_val
    else
        S1_max = funcS_maxwell(n, lambda_val, sgm_val, 2.0, 0.0, 1.0, EPS0)  # S1_mx=S2_mx
        S7_max = funcS_maxwell(n, lambda_val, sgm_val, 2.0, 0.0, -1.0, EPS0)
        S17_max = S1_max - 2.0*sgm_val*lambda_val*S7_max

        b12snj = bj[jj]*n*wcs[s]*S17_max +
            bj[jj]*cj[jj]*kz*vtsz[s]*S1_max*Tsx[s]/Tsz[s]/(1+sgm_val)
    end

    return b12snj

end
