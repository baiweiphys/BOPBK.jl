# Include necessary function files
include("../src/func_Jpole.jl")

J_opt = 8.1
(bzj, czj) = func_Jpole(J_opt)
println("bjz = $bzj")
println("cjz = $czj")
println("sizeOfbzj = ", size(bzj))
println("lengthOfbzj = ", length(bzj))


include("../src/constants.jl")
using .PhysicalConstants

energy = 0.5 * me * c2
println("Energy = ", energy)

include("../src/getLen_SNJ2oneDim.jl")

S_bm = 3
Ns_bm = [0, 1, 2]
J = 2
result_bm = getLen_SNJ2oneDim(S_bm, Ns_bm, J)
println("result: $result_bm")


include("../src/getLen_SNLJ2oneDim.jl")
S_pbk = 2
Ns_pbk = [0, 1]
kappasz_pbk = [2, 4]

result_pbk = getLen_SNLJ2oneDim(S_pbk, Ns_pbk, kappasz_pbk)
println("result: $result_pbk")


include("../src/getIndexOfBlkMatrix_maxwell.jl")

MatrixNo = 2
S_bm = 2
Ns_bm = [0, 1]
J = 3
idx_bm = getIndexOfBlkMatrix_maxwell(S_bm, Ns_bm, J, MatrixNo)
println("idx_bm: $idx_bm")


include("../src/getIndexOfBlkMatrix_mixed.jl")

MatrixNo = 5
S_pbk = 2
Ns_pbk = [0, 1]
kappasz_pbk = [2, 4]
S_bm = 2
Ns_bm = [0, 1]
J = 3
idx_mixed =  getIndexOfBlkMatrix_mixed(S_pbk, S_bm, Ns_pbk, Ns_bm,kappasz_pbk,J,MatrixNo)
println("idx_mixed: $idx_mixed")


include("../src/getIndexOfBlkMatrix_pbk.jl")

MatrixNo = 3
S_pbk = 2
Ns_pbk = [0, 1]
kappasz_pbk = [2, 4]
idx_pbk =  getIndexOfBlkMatrix_pbk(S_pbk, Ns_pbk, kappasz_pbk, MatrixNo)
println("idx_pbk: $idx_pbk")







