"""
Description: Calculate the coefficients of bz30 for the z-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bz30(
    theta::Float64,
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    bz30 = 1im * phys.epsilon0 * sum(wps .^ 2) * tan(theta) .^ 2
    return bz30
end
