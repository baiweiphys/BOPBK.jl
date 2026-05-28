"""
Description: Calculate the coefficients of csn(s) for the oblique plasma 
waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2025-12-10
"""

function pbk_csn(
    s::Int, n::Int,
    kz::Float64,
    kappasz_pbk::AbstractVector{Int},
    vtsz_pbk::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    us0::AbstractVector{Float64})

    csn = n * wcs[s] + kz * us0[s] - 1im * sqrt(kappasz_pbk[s]) * kz * vtsz_pbk[s]

    return csn

end
