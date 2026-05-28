"""
Description: To obtain the index of a subblock matrix within a composite 
matrix that exhibits a bi-Maxwellian plasma distribution.
Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
Date: 2025-12-03
Last Modified: 2026-01-30
"""

function getIndexOfBlkMatrix_bm(
    Ns_bm::AbstractVector{Int},
    J::Int, 
    MatrixNo::Int)


    # Calculate matrix size for bi-Maxwellian species
    len_M_bm = J * (2 * sum(Ns_bm) + length(Ns_bm)) + 1

    if MatrixNo == 1
        # for Mx_maxwell
        idx = 1:len_M_bm
    elseif MatrixNo == 2
        # for My_maxwell
        idx = len_M_bm+1:2*len_M_bm
    elseif MatrixNo == 3
        # for Mz_maxwell
        idx = 2*len_M_bm+1:3*len_M_bm
    else
        error("The bi-Maxwellian MatrixNo must be an integer in the range 1 to 3.")
    end

    return idx
end
