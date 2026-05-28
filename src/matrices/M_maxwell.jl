"""
Description: Calculate matrix M_maxwell for the oblique 
plasma wave model with a loss-cone bi-Maxwellian plasmas.
Optimized solver for Mixed BM plasma (In-place).
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
Created: 2025-12-03
Changelog:
2026-02-08: Optimized M_maxwell (In-place). 
"""

include("SNJ2oneDim.jl")
include("getIndexOfBlkMatrix_bm.jl")

function M_maxwell!(
    M_out::AbstractMatrix{ComplexF64},   # Global pre-allocated matrix
    row_start::Int,   # Start row for this BM block
    S_bm::Int,
    Ns_bm::AbstractVector{Int},
    J::Int,
    csnj, bxyz1snj, bxyz2snj, bxyz3snj,
    MatrixNo::Int,
    ExNo::Int,
    EyNo::Int,
    EzNo::Int)

    @inline idx_M_bm(s, idx_n, j) = SNJ2oneDim(Ns_bm, J, s, idx_n, j)

    # Step 1
    # Calculate block height for BM part
    len_M_bm = J * (2 * sum(Ns_bm) + length(Ns_bm)) + 1

    # Step 2: Get column offset in the GLOBAL matrix
    idx_BlkMatrix = getIndexOfBlkMatrix_bm(Ns_bm, J, MatrixNo)
    firstIdx = idx_BlkMatrix[1] - 1

    # Step 3: Global column indices for fields (Ex, Ey, Ez)
    total_cols = size(M_out, 2)
    col_ex = total_cols - ExNo
    col_ey = total_cols - EyNo
    col_ez = total_cols - EzNo

    # The summation row (last row of this BM block)
    row_summation_global = row_start + len_M_bm - 1

    
    # Step 4: Assemble Matrix
    @inbounds for s = 1:S_bm
        N_bm = Ns_bm[s]
        for (idx_n, n) in enumerate(-N_bm:N_bm)
            for jj = 1:J
                # local indices
                snj = idx_M_bm(s, idx_n, jj)
                # global indices
                current_row = row_start + snj - 1
                target_col = firstIdx + snj
                # Cache physics values
                csnj_val = csnj(s, n, jj)
                bxyz1snj_val = bxyz1snj(s, n, jj)
                bxyz2snj_val = bxyz2snj(s, n, jj)
                bxyz3snj_val = bxyz3snj(s, n, jj)

                # update main row (xyz_{snj})
                M_out[current_row, target_col] += csnj_val
                M_out[current_row, col_ex] += bxyz1snj_val
                M_out[current_row, col_ey] += bxyz2snj_val
                M_out[current_row, col_ez] += bxyz3snj_val

                # update. summation row (dJxyz2)
                M_out[row_summation_global, target_col] += csnj_val
                M_out[row_summation_global, col_ex] += bxyz1snj_val
                M_out[row_summation_global, col_ey] += bxyz2snj_val
                M_out[row_summation_global, col_ez] += bxyz3snj_val
            end
        end
    end

    return len_M_bm
end
