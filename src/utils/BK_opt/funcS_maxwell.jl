"""
Description: Compute the integral of funcS_maxwell for 
maxwellian of loss-cone PBK model.
Author: Bai Wei (baiweiphys@gmail.com)
Created: 2025-12-20
Changelog: 
2026-02-28: Optimized fast_powe: Add fast and safe handling for Float-typed integers.

example:
funcS = @(Jnum,dJnum,num) funcS_maxwell(s,n,lambdas,sgms,Jnum,dJnum,num,EPS0)
S1_maxwell = S2_maxwell = funcS(2,0,1)
S3_maxwell = S4_maxwell = funcS(1,1,2)
S5_maxwell = S6_maxwell = funcS(0,2,3)
S7_maxwell = funcS(2,0,-1)
S8_maxwell = funcS(1,1,0)
S9_maxwell = funcS(0,2,1)

Jnum: the Jnum-th power of Bessel function
Jpnum: the Jpnum-th power of the derivative of Bessel function
num: numerator number
"""

using Bessels: besselj        # 3-10x faster than SpecialFunctions for Bessel evaluations
using QuadGK         # High-performance adaptive integration
using SpecialFunctions: loggamma  # Stable gamma calculation


# Optimized fast_powe: Add fast and safe handling for Float-typed integers (e.g., 3.0)
@inline function fast_pow(x, p::Float64)
    p == 0.0 && return 1.0
    p == 1.0 && return x
    p == 2.0 && return x * x
    p == 3.0 && return x * x * x
    p == 4.0 && return (x2 == x * x; x2 * x2)

    # In physical models, exponents are often given as floats like 2.0 or 3.0.
    # Casting to Int prevents Domain Error when evaluating x^Float64 for x<0.
    isinteger(p) && return x^Int(p)
    # Fall back for non-integer powers:
    # Prevents Domain error by catching negative x, return NaN if undefined.
    return x > 0.0 ? x^p : (x == 0.0 ? 0.0 : NaN)
end



function funcS_maxwell(
    n::Float64,        # Order of Bessel function
    lambda::Float64, 
    sgm::Float64, 
    Jnum::Float64, 
    dJnum::Float64, 
    num::Float64, 
    EPS0::Float64;
    rtol=1.0e-10, atol=0.0) # Use 1e-10 for speed; 1e-16 is usually overkill

    # 1. Precompute the global coefficient S0
    # Use loggamma to prevent overflow and improve numerical stability
    lS0 = log(4.0) - (sgm + 2.0) * log(2.0 * lambda) - loggamma(sgm + 1.0)
    S0 = exp(lS0)

    # Precompute exponents for the integrand
    exp_x = 2.0 * sgm + num

    if lambda > EPS0
        # Case 1: lambda > EPS0
        inv_2lambda = 0.5 / lambda

        f_integrand = (x) -> begin
            # Avoid singularity at x=0; most physics kernels go to 0 here
            if x < 1e-20
                return 0.0
            end

            # Optimization: Calculate Jn and use recurrence for derivative
            # Original: 0.5*(J_{n-1} - J_{n+1}) needs 2 calls. 
            # With Jn, total was 3 calls. Now we only need 2 calls total: Jn and J_{n+1}.
            jn = besselj(n, x)
            # if n == 0
            #     djn = -besselj(1, x)
            # else
            #     djn = (n / x) * jn - besselj(n + 1, x)
            # end
            # [Speed & Accuracy Secret]: Conditional derivative calculation
            # For extremely small x (< 1e-8), the standard formula `(n/x)*jn` can cause 
            # division-by-zero or massive precision loss. We use the safe, division-free 
            # formula here. For all other x, we use the faster standard formula to save 
            # one expensive `besselj` evaluation per integration step!
            djn = if x < 1e-8
                n == 0 ? -besselj(1, x) : 0.5 * (besselj(n - 1, x) - besselj(n + 1, x))
            else
                (n / x) * jn - besselj(n + 1, x)
            end

            # Core kernel calculation
            term1 = fast_pow(x, exp_x)
            term2 = fast_pow(jn, Jnum)
            term3 = fast_pow(djn, dJnum)
            term4 = exp(-x^2 * inv_2lambda)
            return term1 * term2 * term3 * term4
        end
        res, _ = quadgk(f_integrand, 0.0, Inf, rtol=rtol, atol=atol)
    else
        # Case 2: lambda <= EPS0
        sqrt_2lam = sqrt(2.0 * lambda)

        f_integrand = (x) -> begin
            if x < 1e-20
                return 0.0
            end

            arg = x * sqrt_2lam
            jn = besselj(n, arg)
            # if n == 0
            #     djn = -besselj(1, arg)
            # else
            #     # Recurrence for derivative: Jn'(arg)
            #     djn = (n / arg) * jn - besselj(n + 1, arg)
            # end
            # [Crucial Fix]: Division-free derivative formula.
            # Avoids `(n/x)*jn - besselj(n+1,x)` which causes underflow/division-by-zero near x=0.
            djn = n == 0 ? -besselj(1, arg) : 0.5 * (besselj(n - 1, arg) - besselj(n + 1, arg))


            term1 = fast_pow(arg, exp_x) * sqrt_2lam
            term2 = fast_pow(jn, Jnum)
            term3 = fast_pow(djn, dJnum)
            term4 = exp(-x^2)
            return term1 * term2 * term3 * term4
        end

        # Integration limit management
        # Note: quadgk handles Inf well, but scaling the limit can improve stability
        res, _ = quadgk(f_integrand, 0.0, Inf, rtol=rtol, atol=atol)
    end

    return S0 * res
end
