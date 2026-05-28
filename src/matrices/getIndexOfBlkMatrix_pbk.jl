"""
Description: To obtain the index of a subblock matrix within a composite 
matrix that exhibits a PBK plasma distribution.
Filename: getIndexOfBlkMatrix_pbk.m
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2025-12-03
Last Modified: 2026-01-30
"""

function getIndexOfBlkMatrix_pbk(
    Ns_pbk::AbstractVector{Int},
    kappasz_pbk::AbstractVector{Int},
    MatrixNo::Int)

    # Calculate matrix size for PBK species
    len_M_pbk = sum((2n + 1) * (k + 4) * (k + 1) for (n, k) in zip(Ns_pbk, kappasz_pbk)) ÷ 2 + 1


    if MatrixNo == 1
        # for dJx of BPK
        idx = 1:len_M_pbk
    elseif MatrixNo == 2
        # for dJy of BPK
        idx = len_M_pbk+1:2*len_M_pbk
    elseif MatrixNo == 3
        #  for dJz of BPK     
        idx = 2*len_M_pbk+1:3*len_M_pbk
    else
        error("The PBK MatrixNo must be an integer in the range 1 to 3.")
    end

    return idx

end
