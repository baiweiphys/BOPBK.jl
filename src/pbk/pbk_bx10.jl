"""
Description: Calculate the coefficients of bx10 for the x-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bx10(wps::AbstractVector{Float64}, phys::PlasmaConfig)

    # coefs = (sgms+1).*kappasx./(kappasx+sgms)
    # bx10 = 1im*epsilon0*sum(wps.^2.*coefs)

    # FIX: Removed redundant coefs (baiwei 2025-06-16)
    bx10 = 1im * phys.epsilon0 * sum(wps .^ 2)

    return bx10

end