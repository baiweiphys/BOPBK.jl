"""
Description: Calculate the 1-dimensional PBK array index corresponding 
to the given parameters (N_pbk,kappasz_pbk,s,idx_n,l,j)
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2025-02-06
Last Modified: 2026-02-01
"""

function SNLJ2oneDim(
    Ns_pbk::AbstractVector{Int},
    kappasz_pbk::AbstractVector{Int},
    s::Int, idx_n::Int, l::Int, j::Int)
    
    if (idx_n < 1 || idx_n > 2 * Ns_pbk[s] + 1)
        error("idx_n out of range")
    elseif (l < 1 || l > kappasz_pbk[s]+1)
        error("l out of range")
    elseif (j < 1 || j > l + 1)
        error("j out of range")
    end

    
    # 1. Count the elements at the species s-1.
    soff_pbk = 0
    for idx_s in 1:(s-1)
        soff_pbk += (2 * Ns_pbk[idx_s] + 1) * (kappasz_pbk[idx_s] + 1) * (kappasz_pbk[idx_s] + 4) ÷ 2
    end

    # 2. Count the elements at the index n-1.
    noff_pbk = (idx_n - 1) * (kappasz_pbk[s] + 4) * (kappasz_pbk[s] + 1) ÷ 2

    # 3. Calculate the offset for the index l-1.
    loff_pbk = (l + 2) * (l - 1) ÷ 2

    idx = soff_pbk + noff_pbk + loff_pbk + j

    # println("1-dimensional array index: $idx")
    return idx

end
