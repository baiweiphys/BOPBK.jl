"""
Description: To obtain the index of a subblock matrix within a composite 
matrix that exhibits a mixed plasma distribution of PBK and bi-Maxwellian.
Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
Date: 2025-12-03
Last Modified: 2026-02-01
"""

function getIndexOfBlkMatrix_mixed(
    Ns_pbk::AbstractVector{Int},
    Ns_bm::AbstractVector{Int},
    kappasz_pbk::AbstractVector{Int},
    J::Int, 
    MatrixNo::Int)
    
    
    # Calculate matrix size for PBK species
    len_M_pbk = sum((2n + 1) * (k + 4) * (k + 1) for (n, k) in zip(Ns_pbk, kappasz_pbk)) ÷ 2 + 1

    # Calculate matrix size for bi-Maxwellian species
    len_M_bm = J * (2 * sum(Ns_bm) + length(Ns_bm)) + 1

    # The indices 1 to 3 correspond to PBK matrix
    # The indices 4 to 6 correspond to bi-Maxwellian matrix
    if MatrixNo == 1
        # for dJx of BPK
        idx = 1 : len_M_pbk
    elseif  MatrixNo == 2
        # for dJy of BPK
        idx = len_M_pbk+1 : 2*len_M_pbk
    elseif MatrixNo == 3
        # for dJz of BPK
        idx = 2*len_M_pbk+1 : 3*len_M_pbk
    elseif MatrixNo == 4
        # for dJx of maxwellian
        idx = 3*len_M_pbk+1 : 3*len_M_pbk+len_M_bm
    elseif MatrixNo == 5
        # for dJy of maxwellian
        idx = 3*len_M_pbk+len_M_bm+1 : 3*len_M_pbk+2*len_M_bm
    elseif MatrixNo == 6
        # for dJz of maxwellian
        idx = 3*len_M_pbk+2*len_M_bm+1 : 3*len_M_pbk+3*len_M_bm
    else
        error("The mixed MatrixNo must be an integer in the range 1 to 6.")
    end

    return idx
end