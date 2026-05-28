"""
Description: To obtain the size of the complete matrix.
Filename: getSizeOfAllMatrix.m
Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
Date: 2024-09-20
Last Modified: 2026-01-30
"""

include("getIndexOfBlkMatrix_pbk.jl")
include("getIndexOfBlkMatrix_bm.jl")
include("getIndexOfBlkMatrix_mixed.jl")

function getSizeOfAllMatrix(
    Ns::AbstractVector{Int},
    J::Int,
    pbk_mask::AbstractVector{Bool},
    bm_mask::AbstractVector{Bool},
    kappasz::AbstractVector{Int})

    S_pbk = count(pbk_mask)
    S_bm = count(bm_mask)
    Ns_pbk = Ns[pbk_mask]
    Ns_bm = Ns[bm_mask]
    kappasz_pbk = kappasz[pbk_mask]
    len = 0

    if S_bm == 0
        MatrixNo_pbk = 3
        tmp_pbk = getIndexOfBlkMatrix_pbk(Ns_pbk, kappasz_pbk, MatrixNo_pbk)
        len = tmp_pbk[end] + 6
    elseif S_pbk == 0
        MatrixNo_bm = 3
        tmp_bm = getIndexOfBlkMatrix_bm(Ns_bm, J, MatrixNo_bm)
        len = tmp_bm[end] + 9
    else
        MatrixNo_mixed = 6
        tmp_mixed = getIndexOfBlkMatrix_mixed(Ns_pbk, Ns_bm, kappasz_pbk, J, MatrixNo_mixed)
        len = tmp_mixed[end] + 9
    end

    return len

end






