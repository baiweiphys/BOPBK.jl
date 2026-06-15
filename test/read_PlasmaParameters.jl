
using CSV 
using DataFrames


include("../src/constants.jl")
using .PhysicalConstants

# par = importdata("./bopbk.in", " ", 1)


B0 = 5e-5 # Tesla

include("../src/getPlasmaParameters.jl")


# using CSV, DataFrames
df = CSV.read("./bopbk.in", DataFrame;
    delim=' ',
    header=true,
    ignorerepeated=true,
    types=Float64)

par_data = Matrix(df)

(S, Ns, index_pbk, index_bm, 
kappasz, kappasx, vtsz, vtsx,
Tsz, Tsx, sgms, wps, wcs,us0, 
rhocs, lambdaDs) = getPlasmaParameters(B0, par_data)
# println(plasma_params[1])
println(index_pbk)
println(typeof(index_pbk))
# println(Ns)
# println(typeof(Ns))
# println("成功读取数据:")
# println("行数: ", nrow(df))
# println("列数: ", ncol(df))
# println("列名: ", names(df))

# println(df[1, 11])
# println(par_data)
# show(first(df, 2), allcols=true)



