using SpecialFunctions
using QuadGK
using Plots
# using LinearAlgebra

include("../src/funcS_maxwell.jl")


num = 100
lambda_min = 1.0e-2
lambda_max = 1.0e-1
lambdas = range(lambda_min, lambda_max, length=num)

s = 1
n = 6
sgm = 0.0
EPS0 = 1e-1

S1_max = lambda -> funcS_maxwell(n, lambda, sgm, 2.0, 0.0, 1.0, EPS0)
S3_max = lambda -> funcS_maxwell(n, lambda, sgm, 1.0, 1.0, 2.0, EPS0)
S5_max = lambda -> funcS_maxwell(n, lambda, sgm, 0.0, 2.0, 3.0, EPS0)
S7_max = lambda -> funcS_maxwell(n, lambda, sgm, 2.0, 0.0, -1.0, EPS0)
S8_max = lambda -> funcS_maxwell(n, lambda, sgm, 1.0, 1.0, 0.0, EPS0)
S9_max = lambda -> funcS_maxwell(n, lambda, sgm, 0.0, 2.0, 1.0, EPS0)


S1 = Vector{Float64}(undef, num)
S3 = Vector{Float64}(undef, num)
S5 = Vector{Float64}(undef, num)
S7 = Vector{Float64}(undef, num)
S8 = Vector{Float64}(undef, num)
S9 = Vector{Float64}(undef, num)

println("Starting calculation...")
for i in 1:num
    S1[i] = S1_max(lambdas[i])
    S3[i] = S3_max(lambdas[i])
    S5[i] = S5_max(lambdas[i])
    S7[i] = S7_max(lambdas[i])
    S8[i] = S8_max(lambdas[i])
    S9[i] = S9_max(lambdas[i])
end


gr()  # MATLAB风格

# plots
h = plot(size=(600, 600), dpi=100)  
plot!(fontfamily="Arial", fontsize=12)

plot!(lambdas, real.(S1), linewidth=3, label="S1")
plot!(lambdas, real.(S3), linewidth=3, label="S3")
plot!(lambdas, real.(S5), linewidth=3, label="S5")
plot!(lambdas, real.(S7), linewidth=3, label="S7")
plot!(lambdas, real.(S8), linewidth=3, label="S8")
plot!(lambdas, real.(S9), linewidth=3, label="S9")

plot!(legend=:best, legendfontsize=10)
xlabel!("λ")
ylabel!("Re(S_j)")

display(h)

savefig(h, "plt_funcS_max.pdf")

