"""
Description: Calculate the coefficients of bx30 for the x-component 
of oblique plasma waves with a loss-cone PBK distribution.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-09-24
Last Modified: 2026-02-03
"""

function pbk_bx30(
    theta::Float64, 
    wps::AbstractVector{Float64},
    phys::PlasmaConfig)

    bx30 = -1im * phys.epsilon0 * tan(theta) * sum(wps .^ 2)

    return bx30
end
