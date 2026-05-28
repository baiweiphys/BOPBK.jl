"""
Description: Compute the integral of funcS_maxwell for 
maxwellian of loss-cone PBK model.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-11-20
Last Modified: 2025-12-04

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


# Helper function to avoid the expensive generic ^ operator for common exponents
@inline function fast_pow(x, p)
    if p == 0.0 return 1.0 end
    if p == 1.0 return x end
    if p == 2.0 return x * x end
    return x^p
end


function funcS_maxwell(
    n::Float64,        # Order of Bessel function
    lambda::Float64, 
    sgm::Float64, 
    Jnum::Float64, 
    dJnum::Float64, 
    num::Float64, 
    EPS0::Float64;
    rtol=1e-10, atol=1e-12) # Use 1e-10 for speed; 1e-16 is usually overkill

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
            if n == 0
                djn = -besselj(1, x)
            else
                djn = (n / x) * jn - besselj(n + 1, x)
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
            if n == 0
                djn = -besselj(1, arg)
            else
                # Recurrence for derivative: Jn'(arg)
                djn = (n / arg) * jn - besselj(n + 1, arg)
            end

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
