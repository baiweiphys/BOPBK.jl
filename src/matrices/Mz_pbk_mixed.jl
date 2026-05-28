"""
Description: Assembles the Mz_pbk_mixed block within the global Mixed matrix. 
This function is designed to compute the M matrix used in the analysis of oblique 
plasma waves in a mixed distribution of loss-cone PBK and bi-Maxwellian plasmas. 
Optimized solver for Mixed PBK plasma (In-place).
Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
Created: 2025-12-06
Changelog:
2026-02-08: Optimized Mz_pbk_mixed block. 
"""

# include("SNLJ2oneDim.jl")
include("getIndexOfBlkMatrix_mixed.jl")


function Mz_pbk_mixed!(
    M_out::AbstractMatrix{ComplexF64},   # Global pre-allocated matrix
    row_start::Int,
    S_pbk::Int,
    Ns_pbk::AbstractVector{Int},
    Ns_bm::AbstractVector{Int},
    kappasz_pbk::AbstractVector{Int},
    J::Int,
    csn_pbk,
    bz11snl, bz21snl, bz31snl,
    bz13snl, bz23snl, bz33snl,
    bz12snl, bz22snl, bz32snl,
    bz34snl,
    by20::Complex{Float64},
    MatrixNo::Int,
    ExNo::Int, EyNo::Int, EzNo::Int)


    # ExyzNo = 5 for Ex
    # ExyzNo = 4 for Ey
    # ExyzNo = 3 for Ez

    @inline idx_M_pbk(s, idx_n, l, j) = SNLJ2oneDim(Ns_pbk, kappasz_pbk, s, idx_n, l, j)

    # Step 1: Calculate PBK block height (number of rows)
    len_M_pbk = sum((2n + 1) * (k + 4) * (k + 1) for (n, k) in zip(Ns_pbk, kappasz_pbk)) ÷ 2 + 1


    # Step 2: Global column offsets
    idx_BlkMatrix = getIndexOfBlkMatrix_mixed(Ns_pbk, Ns_bm, kappasz_pbk, J, MatrixNo)
    firstIdx = idx_BlkMatrix[1] - 1


    # Step 3: Pre-calculate column indices
    total_cols = size(M_out, 2)
    col_ex = total_cols - ExNo
    col_ey = total_cols - EyNo
    col_ez = total_cols - EzNo

    # The summation row (last row of this PBK block)
    row_summation_global = row_start + len_M_pbk - 1


    # Step 4: Assemble Matrix
    @inbounds for s = 1:S_pbk
        kappaz_pbk = kappasz_pbk[s]
        N_pbk = Ns_pbk[s]
        for (idx_n, n) in enumerate(-N_pbk:N_pbk)
            csn_val = csn_pbk(s, n)
            for l = 1:kappaz_pbk+1
                # Cache physics values
                bz11snl_val = bz11snl(s, n, l)
                bz21snl_val = bz21snl(s, n, l)
                bz31snl_val = bz31snl(s, n, l)
                bz12snl_val = bz12snl(s, n, l)
                bz22snl_val = bz22snl(s, n, l)
                bz32snl_val = bz32snl(s, n, l)
                for jj = 1:l+1
                    # Local index within the block
                    snlj = idx_M_pbk(s, idx_n, l, jj)
                    # Global indices
                    current_row = row_start + snlj - 1
                    target_col = firstIdx + snlj
                    target_col_next = target_col + 1
                    # update Diagonnal-like term
                    M_out[current_row, target_col] += csn_val
                    if (jj < l + 1)
                        # snljp1 = idx_Mpbk(s,idx_n,l,jj+1);
                        M_out[current_row, target_col_next] += 1.0
                    end
                    if (jj == l)
                        if (l <= kappaz_pbk - 1)
                            M_out[current_row, col_ez] += bz34snl(s, n, l + 2)
                            M_out[current_row, col_ex] += bz13snl(s, n, l + 1)
                            M_out[current_row, col_ey] += bz23snl(s, n, l + 1)
                            M_out[current_row, col_ez] += bz33snl(s, n, l + 1)
                            M_out[current_row, col_ex] += bz11snl_val
                            M_out[current_row, col_ey] += bz21snl_val
                            M_out[current_row, col_ez] += bz31snl_val
                        elseif (l == kappaz_pbk)
                            M_out[current_row, col_ex] += bz13snl(s, n, l + 1)
                            M_out[current_row, col_ey] += bz23snl(s, n, l + 1)
                            M_out[current_row, col_ez] += bz33snl(s, n, l + 1)
                            M_out[current_row, col_ex] += bz11snl_val
                            M_out[current_row, col_ey] += bz21snl_val
                            M_out[current_row, col_ez] += bz31snl_val
                        elseif (l == kappaz_pbk + 1)
                            M_out[current_row, col_ex] += bz11snl_val
                            M_out[current_row, col_ey] += bz21snl_val
                            M_out[current_row, col_ez] += bz31snl_val
                        end
                    elseif (jj == l + 1)
                        M_out[current_row, col_ex] += bz12snl_val
                        M_out[current_row, col_ey] += bz22snl_val
                        M_out[current_row, col_ez] += bz32snl_val
                    end
                end
            end
        end
    end

    # step 5: jxyz summation row
    @inbounds for s = 1:S_pbk
        kappaz_pbk = kappasz_pbk[s]
        N_pbk = Ns_pbk[s]
        for idx_n = 1:(2*N_pbk+1)
            for l = 1:kappaz_pbk+1
                snl1 = idx_M_pbk(s, idx_n, l, 1)
                M_out[row_summation_global, firstIdx+snl1] += 1.0
            end
        end
    end

    # for z-dirction:
    # by20 = 1i*\epsilon_0*sum_{s}\omega_{ps}^2
    M_out[row_summation_global, col_ez] += by20

    return len_M_pbk
end