"""
Description: Calculate the 1-dimensional bi-Maxwellian array index corresponding 
to the given parameters (N_max, J, s, idx_n, j)
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2025-02-06
Last Modified: 2026-02-01
"""

function SNJ2oneDim(
    Ns_bm::AbstractVector{Int}, 
    J::Int,
    s::Int, idx_n::Int, j::Int)

    # Check input parameters
    if (idx_n < 1 || idx_n > 2*Ns_bm[s] + 1)
        error("n must be in the range 1 to 2 * Ns_bm[s] + 1")
    elseif (j < 1 || j > J)
        error("j must be in the range 1 to J")
    end


    # 1. Count the elements at the species s-1.
    soff_bm = 0
    for idx_s in 1:(s-1)
        soff_bm += (2 * Ns_bm[idx_s] + 1) * J
    end

    # 2. Calculate the elements at the idx_n-1
    noff_bm = (idx_n - 1) * J

    idx = soff_bm + noff_bm + j

    # println("1-dimensional array index: $idx")
    return idx

end