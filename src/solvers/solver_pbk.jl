"""
Description: Optimized solver for PBK plasma. (In-place).
Uses optimized M_maxwell! functions to assemble the global matrix for the PBK plasma case.
Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
Created: 2025-12-16
Changelog:
2026-02-08: Optimized solver_pbk. 
"""


const pbk_path = abspath(joinpath(@__DIR__, "..", "pbk"))
include(joinpath(pbk_path, "pbk_csn.jl"))
#
include(joinpath(pbk_path, "pbk_b1snl.jl"))
include(joinpath(pbk_path, "pbk_b2snl.jl"))
include(joinpath(pbk_path, "pbk_b3snl.jl"))
include(joinpath(pbk_path, "pbk_b4snl.jl"))
include(joinpath(pbk_path, "pbk_b5snl.jl"))
include(joinpath(pbk_path, "pbk_b6snl.jl"))
#
include(joinpath(pbk_path, "pbk_bx10.jl"))
include(joinpath(pbk_path, "pbk_bx11snl.jl"))
include(joinpath(pbk_path, "pbk_bx12snl.jl"))
include(joinpath(pbk_path, "pbk_bx21snl.jl"))
include(joinpath(pbk_path, "pbk_bx22snl.jl"))
include(joinpath(pbk_path, "pbk_bx31snl.jl"))
include(joinpath(pbk_path, "pbk_bx32snl.jl"))
include(joinpath(pbk_path, "pbk_bx33snl.jl"))
#
include(joinpath(pbk_path, "pbk_by21snl.jl"))
include(joinpath(pbk_path, "pbk_by22snl.jl"))
include(joinpath(pbk_path, "pbk_by31snl.jl"))
include(joinpath(pbk_path, "pbk_by32snl.jl"))
include(joinpath(pbk_path, "pbk_by33snl.jl"))
#
include(joinpath(pbk_path, "pbk_bz31snl.jl"))
include(joinpath(pbk_path, "pbk_bz32snl.jl"))
include(joinpath(pbk_path, "pbk_bz33snl.jl"))
include(joinpath(pbk_path, "pbk_bz34snl.jl"))

# matrices
const matrices_path = abspath(joinpath(@__DIR__, "..", "matrices"))
include(joinpath(matrices_path, "Mxy_pbk.jl"))
include(joinpath(matrices_path, "Mz_pbk.jl"))


function solver_pbk!(
    cache,
    kx::Float64,
    kz::Float64,
    theta::Float64,
    PlasmaParameters,
    select_eigsolver::Int,
    EPS0::Float64,
    phys::PlasmaConfig)


    # Parameters
    (; S, Ns, J_opt) = PlasmaParameters.meta
    # (; is_pbk, is_bm) = PlasmaParameters.flags
    (; Tsz, Tsx, vtsz, vtsx, us0) = PlasmaParameters.thermo
    (; wps, wcs, rhocs, lambdaDs, ms, ns0) = PlasmaParameters.physics
    (; kappasz, kappasx, sgms) = PlasmaParameters.dist_params
    J = floor(Int, J_opt)
    lambdas = 0.5 * kx^2 * rhocs .^ 2  # for argument of the modified bessel function
    #
    S_pbk = S
    Ns_pbk = Ns
    kappasz_pbk = kappasz


    # PBK functions
    @inline csn(s, n) = pbk_csn(s, n, kz, kappasz, vtsz, wcs, us0)
    @inline b1snl(s, n, l) = pbk_b1snl(s, n, l, kz, kappasz, kappasx, vtsz, sgms, wcs, lambdas, EPS0)
    @inline b2snl(s, n, l) = pbk_b2snl(s, n, l, kz, kappasz, kappasx, vtsz, vtsx, sgms, lambdas, EPS0)
    @inline b3snl(s, n, l) = pbk_b3snl(s, n, l, kz, kappasz, kappasx, vtsz, sgms, wcs, lambdas, EPS0)
    @inline b4snl(s, n, l) = pbk_b4snl(s, n, l, kz, kappasz, kappasx, vtsz, vtsx, sgms, lambdas, EPS0)
    @inline b5snl(s, n, l) = pbk_b5snl(s, n, l, kz, kappasz, kappasx, vtsz, sgms, wcs, lambdas, EPS0)
    @inline b6snl(s, n, l) = pbk_b6snl(s, n, l, kz, kappasz, kappasx, vtsz, vtsx, sgms, lambdas, EPS0)

    # for x-component
    bx10_pbk = pbk_bx10(wps, phys)
    # bx30_pbk = pbk_bx30(theta,wps)
    @inline bx11snl(s, n, l) = pbk_bx11snl(s, n, l, b1snl, wps, phys)
    @inline bx12snl(s, n, l) = pbk_bx12snl(s, n, l, b2snl, wps, phys)
    @inline bx21snl(s, n, l) = pbk_bx21snl(s, n, l, b3snl, wps, phys)
    @inline bx22snl(s, n, l) = pbk_bx22snl(s, n, l, b4snl, wps, phys)
    @inline bx31snl(s, n, l) = pbk_bx31snl(s, n, l, theta, csn, b1snl, b2snl, wps, wcs, phys)
    @inline bx32snl(s, n, l) = pbk_bx32snl(s, n, l, theta, csn, b2snl, wps, wcs, phys)
    @inline bx33snl(s, n, l) = pbk_bx33snl(s, n, l, theta, b1snl, wps, wcs, phys)

    # for y-component
    @inline by11snl(s, n, l) = -1.0 * bx21snl(s, n, l)
    @inline by12snl(s, n, l) = -1.0 * bx22snl(s, n, l)
    by20_pbk = 1im * phys.epsilon0 * sum(wps .^ 2)
    @inline by21snl(s, n, l) = pbk_by21snl(s, n, l, b5snl, wps, phys)
    @inline by22snl(s, n, l) = pbk_by22snl(s, n, l, b6snl, wps, phys)
    @inline by31snl(s, n, l) = pbk_by31snl(s, n, l, theta, csn, b3snl, b4snl, wps, wcs, phys)
    @inline by32snl(s, n, l) = pbk_by32snl(s, n, l, theta, csn, b4snl, wps, wcs, phys)
    @inline by33snl(s, n, l) = pbk_by33snl(s, n, l, theta, b3snl, wps, wcs, phys)


    # for z-component
    # bz10_pbk = bx30_pbk;
    bz11snl = bx31snl
    bz12snl = bx32snl
    bz13snl = bx33snl
    @inline bz21snl(s, n, l) = -1.0 * by31snl(s, n, l)
    @inline bz22snl(s, n, l) = -1.0 * by32snl(s, n, l)
    @inline bz23snl(s, n, l) = -1.0 * by33snl(s, n, l)
    # bz30_pbk = pbk_bz30(theta,wps);
    @inline bz31snl(s, n, l) = pbk_bz31snl(s, n, l, theta, csn, b1snl, b2snl, wps, wcs, phys)
    @inline bz32snl(s, n, l) = pbk_bz32snl(s, n, l, theta, csn, b2snl, wps, wcs, phys)
    @inline bz33snl(s, n, l) = pbk_bz33snl(s, n, l, theta, csn, b1snl, b2snl, wps, wcs, phys)
    @inline bz34snl(s, n, l) = pbk_bz34snl(s, n, l, theta, b1snl, wps, wcs, phys)


    # MEMORY ALLOCATION FOR THE GLOBAL MATRIX
    if !cache.is_initialized
        # 1. Calc PBK size
        len_M_pbk = sum((2n + 1) * (k + 4) * (k + 1) for (n, k) in zip(Ns_pbk, kappasz_pbk)) ÷ 2 + 1
        # 2. Total structure
        total_rows = 3 * len_M_pbk + 6
        total_cols = total_rows

        cache.M = zeros(ComplexF64, total_rows, total_cols)
        cache.is_initialized = true
    end

    M = cache.M
    fill!(M, 0.0 + 0.0im)  # Reset the matrix to zero before filling it with new values.

    # ASSEMBLE THE GLOBAL MATRIX
    ##################
    # ExyzNo = 5 for Ex
    # ExyzNo = 4 for Ey
    # ExyzNo = 3 for Ez
    # BxyzNo = 2 for Bx
    # BxyzNo = 1 for By
    # BxyzNo = 0 for Bz
    ##################  
    current_row = 1  # To keep track of the current row index for filling the matrix

    # Step 1: Fill the PBK blocks
    # Mx PBK (MatrixNo=1)
    h_x_pbk = Mxy_pbk!(M, current_row, S_pbk, Ns_pbk, kappasz_pbk,
        csn, bx11snl, bx21snl, bx31snl,
        bx12snl, bx22snl, bx32snl, bx33snl,
        bx10_pbk, 1, 5, 4, 3, 5)
    current_row += h_x_pbk

    # My PBK (MatrixNo=2)
    h_y_pbk = Mxy_pbk!(M, current_row, S_pbk, Ns_pbk, kappasz_pbk,
        csn, by11snl, by21snl, by31snl,
        by12snl, by22snl, by32snl, by33snl,
        by20_pbk, 2, 5, 4, 3, 4)
    current_row += h_y_pbk

    # Mz PBK (MatrixNo=3)
    h_z_pbk = Mz_pbk!(M, current_row, S_pbk, Ns_pbk, kappasz_pbk,
        csn, bz11snl, bz21snl, bz31snl,
        bz13snl, bz23snl, bz33snl, bz12snl, bz22snl, bz32snl, bz34snl,
        by20_pbk, 3, 5, 4, 3)
    current_row += h_z_pbk


    # Step 2: Maxwell's equations for Ex
    idx_Jx_pbk = getIndexOfBlkMatrix_pbk(Ns_pbk, kappasz_pbk, 1)
    M[end-5, end-1] = M[end-5, end-1] + phys.c2 * kz # for the first item of Ex
    M[end-5, idx_Jx_pbk[end]] = M[end-5, idx_Jx_pbk[end]] - 1im / phys.epsilon0  # for the second item of Ex

    # Maxwell's equations for Ey
    idx_Jy_pbk = getIndexOfBlkMatrix_pbk(Ns_pbk, kappasz_pbk, 2)
    M[end-4, end-2] = M[end-4, end-2] - phys.c2 * kz  # for the first item of Ey
    M[end-4, end] = M[end-4, end] + phys.c2 * kx  # for the second item of Ey
    M[end-4, idx_Jy_pbk[end]] = M[end-4, idx_Jy_pbk[end]] - 1im / phys.epsilon0 # for the third item of Ey

    # Maxwell's equations Ez
    idx_Jz_pbk = getIndexOfBlkMatrix_pbk(Ns_pbk, kappasz_pbk, 3)
    M[end-3, end-1] = M[end-3, end-1] - phys.c2 * kx  # for the first item of Ez
    M[end-3, idx_Jz_pbk[end]] = M[end-3, idx_Jz_pbk[end]] - 1im / phys.epsilon0  # for the second item of Ez

    # Maxwell's equations for Bx, By and Bz
    M[end-2, end-4] = M[end-2, end-4] - kz   # for Bx
    M[end-1, end-5] = M[end-1, end-5] + kz   # for the first item of By
    M[end-1, end-3] = M[end-1, end-3] - kx   # for the second item of By
    M[end, end-4] = M[end, end-4] + kx       # for Bz

    ##Solver##
    eigenvalues, eigenvectors = eigensolver(M, select_eigsolver)
    # eigenvalues, eigenvectors = eigensolver(M, select_eigsolver=1, nev=10, sort_by=:real_desc)
    # eigenvalues, eigenvectors = eigensolver(M, select_eigsolver=1, nev=10)

    return (eigenvalues, eigenvectors)
end
