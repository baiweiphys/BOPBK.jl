using SpecialFunctions
using QuadGK
using Plots
# using LinearAlgebra

include("../src/funcS_pbk.jl")


# function logspace(start, stop, n=50; base=10)
#     start_exp = log(base, start)
#     stop_exp = log(base, stop)
#     exponents = range(start_exp, stop_exp, length=n)
#     return base .^ exponents
# end


num = 100
lambda_min = 1.0e-5
lambda_max = 1.0e-1
lambdas = range(lambda_min, lambda_max, length=num)
# lambda = logspace(lambda_min, lambda_max, num)

s = 1
n = 5
sgm = 1.0e-9
kappasx = [10]  # kappax >=3
kappax = kappasx[s]
EPS0 = 1.0e-2

S1_pbk = lambda -> funcS_pbk(n,kappax,lambda,sgm,2.0,0.0,1.0,1.0,EPS0);
S2_pbk = lambda -> funcS_pbk(n,kappax,lambda,sgm,2.0,0.0,1.0,0.0,EPS0);
S3_pbk = lambda -> funcS_pbk(n,kappax,lambda,sgm,1.0,1.0,2.0,1.0,EPS0);
S4_pbk = lambda -> funcS_pbk(n,kappax,lambda,sgm,1.0,1.0,2.0,0.0,EPS0);
S5_pbk = lambda -> funcS_pbk(n,kappax,lambda,sgm,0.0,2.0,3.0,1.0,EPS0);
S6_pbk = lambda -> funcS_pbk(n,kappax,lambda,sgm,0.0,2.0,3.0,0.0,EPS0);
S7_pbk = lambda -> funcS_pbk(n,kappax,lambda,sgm,2.0,0.0,-1.0,0.0,EPS0);
S8_pbk = lambda -> funcS_pbk(n,kappax,lambda,sgm,1.0,1.0,0.0,0.0,EPS0);
S9_pbk = lambda -> funcS_pbk(n,kappax,lambda,sgm,0.0,2.0,1.0,0.0,EPS0);

S1 = Vector{Float64}(undef, num)
S2 = Vector{Float64}(undef, num)
S3 = Vector{Float64}(undef, num)
S4 = Vector{Float64}(undef, num)
S5 = Vector{Float64}(undef, num)
S6 = Vector{Float64}(undef, num)
S7 = Vector{Float64}(undef, num)
S8 = Vector{Float64}(undef, num)
S9 = Vector{Float64}(undef, num)


println("Starting calculation...")
for i in 1:num
    S1[i] = S1_pbk(lambdas[i])
    S2[i] = S2_pbk(lambdas[i])
    S3[i] = S3_pbk(lambdas[i])
    S4[i] = S4_pbk(lambdas[i])
    S5[i] = S5_pbk(lambdas[i])
    S6[i] = S6_pbk(lambdas[i])
    S7[i] = S7_pbk(lambdas[i])
    S8[i] = S8_pbk(lambdas[i])
    S9[i] = S9_pbk(lambdas[i])
end


# gr()  # MATLAB风格
# plots
h = plot(size=(600, 600), dpi=100)
plot!(fontfamily="Arial", fontsize=12)

plot!(lambdas, real.(S1), linewidth=3, label="S1")
plot!(lambdas, real.(S2), linewidth=3, label="S2")
plot!(lambdas, real.(S3), linewidth=3, label="S3")
plot!(lambdas, real.(S4), linewidth=3, label="S4")
plot!(lambdas, real.(S5), linewidth=3, label="S5")
plot!(lambdas, real.(S6), linewidth=3, label="S6")
plot!(lambdas, real.(S7), linewidth=3, label="S7")
plot!(lambdas, real.(S8), linewidth=3, label="S8")
plot!(lambdas, real.(S9), linewidth=3, label="S9")

plot!(legend=:best, legendfontsize=10)
xlabel!("λ")
ylabel!("Re(S_j)")

# display(h)
savefig(h, "plt_funcS_pbk.pdf")

