"""
Description: Assembles the M_maxwell_mixed block within the global Mixed matrix. 
This function is designed to compute the M matrix used in the analysis of oblique 
plasma waves in a mixed distribution of loss-cone PBK and bi-Maxwellian plasmas. 
Optimized solver for Mixed BM plasma (In-place).
Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
Created: 2025-12-03
Changelog:
2026-02-08: Optimized M_maxwell_mixed block. 
"""

# include("SNJ2oneDim.jl")
include("getIndexOfBlkMatrix_mixed.jl")


function M_maxwell_mixed!(
    M_out::AbstractMatrix{ComplexF64},   # Global pre-allocated matrix
    row_start::Int,   # Start row for this BM block
    S_bm::Int,
    Ns_pbk::AbstractVector{Int},
    Ns_bm::AbstractVector{Int},
    kappasz_pbk::AbstractVector{Int},
    J::Int,
    csnj, bxyz1snj, bxyz2snj, bxyz3snj,
    MatrixNo::Int,
    ExNo::Int,
    EyNo::Int,
    EzNo::Int)

    # Helper for BM internal indexing
    @inline idx_M_bm(s, idx_n, j) = SNJ2oneDim(Ns_bm, J, s, idx_n, j)

    # Step 1: Calculate BM block height (number of rows)
    len_M_bm = J * (2 * sum(Ns_bm) + length(Ns_bm)) + 1

    # Step 2: Get column offset in the GLOBAL matrix
    idx_BlkMatrix = getIndexOfBlkMatrix_mixed(Ns_pbk, Ns_bm, kappasz_pbk, J, MatrixNo)
    firstIdx = idx_BlkMatrix[1] - 1

    # Step 3: Pre-calculate constant column indices
    total_cols = size(M_out, 2)
    col_ex = total_cols - ExNo
    col_ey = total_cols - EyNo
    col_ez = total_cols - EzNo

    # The summation row (last row of this BM block)
    row_summation_global = row_start + len_M_bm - 1
 

    # Step 4: Create Matrix
    @inbounds for s = 1:S_bm
        N_bm = Ns_bm[s]
        for (idx_n, n) in enumerate(-N_bm:N_bm)
            for jj = 1:J
                # Local index within the block
                snj = idx_M_bm(s, idx_n, jj)
                # Global indices
                current_row = row_start + snj - 1
                target_col = firstIdx + snj
                # Cache physics values
                csnj_val = csnj(s, n, jj)
                bxyz1snj_val = bxyz1snj(s, n, jj)
                bxyz2snj_val = bxyz2snj(s, n, jj)
                bxyz3snj_val = bxyz3snj(s, n, jj)

                # Update Main Row (xyz_{snj})
                M_out[current_row, target_col] += csnj_val
                M_out[current_row, col_ex] += bxyz1snj_val
                M_out[current_row, col_ey] += bxyz2snj_val
                M_out[current_row, col_ez] += bxyz3snj_val

                # Update summation Row (Jxyz2)
                M_out[row_summation_global, target_col] += csnj_val
                M_out[row_summation_global, col_ex] += bxyz1snj_val
                M_out[row_summation_global, col_ey] += bxyz2snj_val
                M_out[row_summation_global, col_ez] += bxyz3snj_val
            end
        end
    end

    return len_M_bm
end
