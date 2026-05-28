"""
Description: Calculate the coefficients of bsnl for the oblique 
plasma waves with a loss-cone product-bi-kappa distribution (for PBK of type 2).
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2025-12-07
"""

function pbk_bsnl(
    s::Int, n::Int, l::Int,
    kappasz::AbstractVector{Int},
    vtsz::AbstractVector{Float64},
    wcs::AbstractVector{Float64},
    kz::Float64)

    # Use gammaln to avoid factorial overflow
    log_csl = loggamma(kappasz[s]+1.0) + loggamma(2.0*kappasz[s]-l+2.0) - 
        loggamma(kappasz[s]-l+2.0) - loggamma(2.0*kappasz[s]+1.0) + 
        (l-1.0)*log(2im*sqrt(kappasz[s])*kz*vtsz[s])

    bsnl = -1.0*n*wcs[s] * exp(log_csl)

    return bsnl
end