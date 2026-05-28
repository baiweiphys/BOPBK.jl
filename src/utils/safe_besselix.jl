"""
Description: Safe wrapper for scaled modified Bessel I_n (besselix).
Falls back to an asymptotic expansion for large arguments to avoid
AMOS overflow/accuracy errors.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2026-02-04
"""

using SpecialFunctions: besselix

const BESSELIX_ASYMPTOTIC_Z = 700.0

@inline function besselix_safe(n::Integer, z::Float64)
    if z > BESSELIX_ASYMPTOTIC_Z
        # Asymptotic for scaled I_n(z): e^{-z} I_n(z)
        # I_n(z) ~ e^{z}/sqrt(2π z) * (1 - (4n^2-1)/(8z) + (4n^2-1)(4n^2-9)/(2!(8z)^2) + ...)
        invsqrt = inv(sqrt(2.0 * pi * z))
        mu = 4.0 * (n^2)
        term1 = 1.0 - (mu - 1.0) / (8.0 * z)
        term2 = (mu - 1.0) * (mu - 9.0) / (2.0 * (8.0 * z)^2)
        return invsqrt * (term1 + term2)
    else
        return besselix(n, z)
    end
end
