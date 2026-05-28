"""
Description: Calculate the coefficients of bsl for the oblique 
plasma waves with a product-bi-kappa distribution (for PBK of type 2).
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2025-12-07
"""

function pbk_bsl(
    s::Int, l::Int,
    kappasz::AbstractVector{Int},
    vtsz::AbstractVector{Float64},
    vtsx::AbstractVector{Float64},
    kz::Float64)

    # Use gammaln to avoid factorial overflow
    log_csl = loggamma(kappasz[s] + 1.0) + loggamma(2.0 * kappasz[s] - l + 2.0) - 
        loggamma(kappasz[s] - l + 2.0) - loggamma(2.0 * kappasz[s] + 1.0) + 
        (l - 1.0) * log(2im * sqrt(kappasz[s]) * kz * vtsz[s])

    bsl = -0.5 * l * kz^2.0 * vtsx[s] .^ 2.0 * exp(log_csl)

    return bsl
end
