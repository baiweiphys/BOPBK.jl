"""
Description: Compute the integral of funcS_pbk for loss-cone PBK model.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-11-20
Last Modified: 2025-12-05

example:
funcS = @(Jnum,dJnum,num,den) funcS_pbk(s,n,kappas_perp,lambdas,sgms,Jnum,dJnum,num,den,EPS0)
S1_pbk = FuncS_pbk(2.0, 0.0, 1.0, 2.0)
S2_pbk = FuncS_pbk(2.0, 0.0, 1.0, 1.0)
S3_pbk = FuncS_pbk(1.0, 1.0, 2.0, 2.0)
S4_pbk = FuncS_pbk(1.0, 1.0, 2.0, 1.0)
S5_pbk = FuncS_pbk(0.0, 2.0, 3.0, 2.0)
S6_pbk = FuncS_pbk(0.0, 2.0, 3.0, 1.0)
S7_pbk = FuncS_pbk(2.0, 0.0, -1.0, 1.0)
S8_pbk = FuncS_pbk(1.0, 1.0, 0.0, 1.0)
S9_pbk = FuncS_pbk(0.0, 2.0, 1.0, 1.0)

% Jnum: the Jnum-th power of Bessel function
% Jpnum: the Jpnum-th power of the derivative of Bessel function
% num: numerator number
% den: denominator number
"""

using Bessels: besselj        # Optimized Bessel functions (faster than SpecialFunctions)
using QuadGK         # Adaptive Gauss-Kronrod integration
using SpecialFunctions: loggamma # Only need loggamma from here


# Helper function to handle small integer powers much faster than generic x^p
@inline function fast_pow(x, p)
    if p == 0.0 return 1.0 end
    if p == 1.0 return x end
    if p == 2.0 return x * x end
    return x^p
end


function funcS_pbk(
    n::Int, 
    kappax::Int, 
    lam::Float64, 
    sgm::Float64, 
    p1::Float64, p2::Float64, 
    p3::Float64, p4::Float64, 
    EPS0::Float64;
    rtol=1e-10, atol=1e-12) # Slightly relaxed tolerances for significant speedup

  
    # 1. Precompute Gamma ratio and coefficients (invariant during integration)
    # Using loggamma prevents overflow for large arguments
    lgamma_ratio = loggamma(kappax + sgm + 1.0) - loggamma(kappax) - loggamma(sgm + 1.0)

    # Precompute exponents for the integrand
    exp_num = 2.0 * sgm + p3
    exp_den = kappax + sgm + p4

    # Pre-determine which branch to use to avoid 'if' overhead inside the integrand
    if lam > EPS0
        # Precompute constants for the lam > EPS0 case
        S0 = 4.0 * (2.0 * lam)^(-sgm - 2.0) * kappax^(-sgm - 1.0) * exp(lgamma_ratio)
        inv_lam_kappa = 0.5 / (lam * kappax)

        # Define the integrand closure
        f_integrand = (x) -> begin
            # Handle x=0 to avoid division by zero or NaN in Bessel/Powers
            if x < 1e-20
                return 0.0
            end

            # Calculate Jn and use recurrence for derivative
            jn = besselj(n, x)
            djn = (n / x) * jn - besselj(n + 1, x)

            # Efficiency: fast_pow is significantly faster for small p values
            num = fast_pow(x, exp_num) * fast_pow(jn, p1) * fast_pow(djn, p2)
            den = (1.0 + x^2 * inv_lam_kappa)^exp_den
            return num / den
        end
    else
        # Precompute constants for the lam <= EPS0 case
        sqrt_2lam = sqrt(2.0 * lam)
        S0 = 4.0 * (2.0 * lam)^(0.5 * p3 - 1.5) * kappax^(-sgm - 1.0) * exp(lgamma_ratio)
        inv_kappa = 1.0 / kappax

        f_integrand = (x) -> begin
            if x < 1e-20
                return 0.0
            end

            arg = x * sqrt_2lam
            jn = besselj(n, arg)
            # Derivative chain rule: d/dx [Jn(a*x)] = a * J'n(ax)
            # Original code used dJn(x*sqrt(2lam)), implying derivative of the function itself
            djn = (n / arg) * jn - besselj(n + 1, arg)

            num = fast_pow(x, exp_num) * fast_pow(jn, p1) * fast_pow(djn, p2)
            den = (1.0 + x^2 * inv_kappa)^exp_den
            return num / den
        end
    end


    # 2. Perform Numerical Integration
    # QuadGK handles (0, Inf) by internal variable transformation
    # Reduced tolerance (e.g., 1e-10) is usually sufficient for physical models
    res, _ = quadgk(f_integrand, 0.0, Inf, rtol=rtol, atol=atol)

    return S0 * res
end

# Example usage:
# S1_pbk = funcS_pbk_optimized(2, 2.0, 1.0, 1.0, 2.0, 0.0, 1.0, 2.0, 1e-16)
