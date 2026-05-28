"""
Description: Optimized solver for mixed plasma. (In-place).
Uses optimized M_xy_pbk_mixed!, Mz_pbk_mixed!, and M_maxwell_mixed! functions 
to assemble the global matrix for the mixed plasma case.
Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
Created: 2025-12-16
Changelog:
2026-02-08: Optimized solver_mixed. 
"""

# pbk
const pbk_path = abspath(joinpath(@__DIR__, "..", "pbk"))
include(joinpath(pbk_path, "pbk_csn.jl"))
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

# maxwell
const maxwell_path = abspath(joinpath(@__DIR__, "..", "maxwell"))
include(joinpath(maxwell_path, "maxwell_csnj.jl"))
include(joinpath(maxwell_path, "maxwell_b12snj.jl"))
include(joinpath(maxwell_path, "maxwell_b34snj.jl"))
include(joinpath(maxwell_path, "maxwell_b56snj.jl"))
#
include(joinpath(maxwell_path, "maxwell_bx1.jl"))
include(joinpath(maxwell_path, "maxwell_bx1snj.jl"))
include(joinpath(maxwell_path, "maxwell_bx2.jl"))
include(joinpath(maxwell_path, "maxwell_bx2snj.jl"))
include(joinpath(maxwell_path, "maxwell_bx3.jl"))
include(joinpath(maxwell_path, "maxwell_bx3snj.jl"))

include(joinpath(maxwell_path, "maxwell_by1.jl"))
include(joinpath(maxwell_path, "maxwell_by1snj.jl"))
include(joinpath(maxwell_path, "maxwell_by2.jl"))
include(joinpath(maxwell_path, "maxwell_by2snj.jl"))
include(joinpath(maxwell_path, "maxwell_by3.jl"))
include(joinpath(maxwell_path, "maxwell_by3snj.jl"))
#
include(joinpath(maxwell_path, "maxwell_bz1.jl"))
include(joinpath(maxwell_path, "maxwell_bz1snj.jl"))
include(joinpath(maxwell_path, "maxwell_bz2.jl"))
include(joinpath(maxwell_path, "maxwell_bz2snj.jl"))
include(joinpath(maxwell_path, "maxwell_bz3.jl"))
include(joinpath(maxwell_path, "maxwell_bz3snj.jl"))

# matrices
const matrices_path = abspath(joinpath(@__DIR__, "..", "matrices"))
include(joinpath(matrices_path, "Mxy_pbk_mixed.jl"))
include(joinpath(matrices_path, "Mz_pbk_mixed.jl"))
include(joinpath(matrices_path, "M_maxwell_mixed.jl"))

# utils
# const utils_path = abspath(joinpath(@__DIR__, "..", "utils"))
# include(joinpath(utils_path, "eigensolver.jl"))


function solver_mixed!(
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
    (; is_pbk, is_bm) = PlasmaParameters.flags
    (; Tsz, Tsx, vtsz, vtsx, us0) = PlasmaParameters.thermo
    (; wps, wcs, rhocs, lambdaDs, ms, ns0) = PlasmaParameters.physics
    (; kappasz, kappasx, sgms) = PlasmaParameters.dist_params
    J = floor(Int, J_opt)
    lambdas = 0.5*kx^2*rhocs.^2  # for argument of the modified bessel function


    # PBK parameters
    S_pbk = count(is_pbk)
    Ns_pbk = Ns[is_pbk]
    kappasz_pbk = kappasz[is_pbk]
    kappasx_pbk = kappasx[is_pbk]
    vtsz_pbk = vtsz[is_pbk]
    vtsx_pbk = vtsx[is_pbk]
    # Tsz_pbk = Tsz[is_pbk]
    # Tsx_pbk = Tsx[is_pbk]
    wps_pbk = wps[is_pbk]
    wcs_pbk = wcs[is_pbk]
    us0_pbk = us0[is_pbk]
    # rhocs_pbk = rhocs[is_pbk]
    # lambdaDs_pbk = lambdaDs[is_pbk]
    lambdas_pbk = lambdas[is_pbk]
    sgms_pbk = sgms[is_pbk]


    # BM parameters
    S_bm = count(is_bm)
    Ns_bm = Ns[is_bm]
    # kappasx_bm = kappasx[is_bm]
    vtsz_bm = vtsz[is_bm]
    # vtsx_bm = vtsx[is_bm]
    Tsz_bm = Tsz[is_bm]
    Tsx_bm = Tsx[is_bm]
    wps_bm = wps[is_bm]
    wcs_bm = wcs[is_bm]
    us0_bm = us0[is_bm]
    # rhocs_bm = rhocs[is_bm]
    # lambdaDs_bm = lambdaDs[is_bm]
    lambdas_bm = lambdas[is_bm]
    sgms_bm = sgms[is_bm]


    # PBK functions
    @inline csn(s,n) = pbk_csn(s,n,kz,kappasz_pbk,vtsz_pbk,wcs_pbk,us0_pbk)
    @inline b1snl(s,n,l) = pbk_b1snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,sgms_pbk,wcs_pbk,lambdas_pbk,EPS0)
    @inline b2snl(s,n,l) = pbk_b2snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,vtsx_pbk,sgms_pbk,lambdas_pbk,EPS0)
    @inline b3snl(s,n,l) = pbk_b3snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,sgms_pbk,wcs_pbk,lambdas_pbk,EPS0)
    @inline b4snl(s,n,l) = pbk_b4snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,vtsx_pbk,sgms_pbk,lambdas_pbk,EPS0)
    @inline b5snl(s,n,l) = pbk_b5snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,sgms_pbk,wcs_pbk,lambdas_pbk,EPS0)
    @inline b6snl(s,n,l) = pbk_b6snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,vtsx_pbk,sgms_pbk,lambdas_pbk,EPS0)

    # x-component
    bx10_pbk = pbk_bx10(wps_pbk,phys)
    # bx30_pbk = pbk_bx30(theta,wps_pbk)
    @inline bx11snl(s,n,l) = pbk_bx11snl(s,n,l,b1snl,wps_pbk,phys)
    @inline bx12snl(s,n,l) = pbk_bx12snl(s,n,l,b2snl,wps_pbk,phys)
    @inline bx21snl(s,n,l) = pbk_bx21snl(s,n,l,b3snl,wps_pbk,phys)
    @inline bx22snl(s,n,l) = pbk_bx22snl(s,n,l,b4snl,wps_pbk,phys)
    @inline bx31snl(s,n,l) = pbk_bx31snl(s,n,l,theta,csn,b1snl,b2snl,wps_pbk,wcs_pbk,phys)
    @inline bx32snl(s,n,l) = pbk_bx32snl(s,n,l,theta,csn,b2snl,wps_pbk,wcs_pbk,phys)
    @inline bx33snl(s,n,l) = pbk_bx33snl(s,n,l,theta,b1snl,wps_pbk,wcs_pbk,phys)
    
    # y-component
    @inline by11snl(s,n,l) = -1.0*bx21snl(s,n,l)
    @inline by12snl(s,n,l) = -1.0*bx22snl(s,n,l)
    by20_pbk = 1im*phys.epsilon0*sum(wps_pbk.^2)
    @inline by21snl(s,n,l) = pbk_by21snl(s,n,l,b5snl,wps_pbk,phys)
    @inline by22snl(s,n,l) = pbk_by22snl(s,n,l,b6snl,wps_pbk,phys)
    @inline by31snl(s,n,l) = pbk_by31snl(s,n,l,theta,csn,b3snl,b4snl,wps_pbk,wcs_pbk,phys)
    @inline by32snl(s,n,l) = pbk_by32snl(s,n,l,theta,csn,b4snl,wps_pbk,wcs_pbk,phys)
    @inline by33snl(s,n,l) = pbk_by33snl(s,n,l,theta,b3snl,wps_pbk,wcs_pbk,phys)
    
    # z-component
    @inline bz11snl(s,n,l) = bx31snl(s,n,l)
    @inline bz12snl(s,n,l) = bx32snl(s,n,l)
    @inline bz13snl(s,n,l) = bx33snl(s,n,l)
    @inline bz21snl(s,n,l) = -1.0*by31snl(s,n,l)
    @inline bz22snl(s,n,l) = -1.0*by32snl(s,n,l)
    @inline bz23snl(s,n,l) = -1.0*by33snl(s,n,l)
    @inline bz31snl(s,n,l) = pbk_bz31snl(s,n,l,theta,csn,b1snl,b2snl,wps_pbk,wcs_pbk,phys)
    @inline bz32snl(s,n,l) = pbk_bz32snl(s,n,l,theta,csn,b2snl,wps_pbk,wcs_pbk,phys)
    @inline bz33snl(s,n,l) = pbk_bz33snl(s,n,l,theta,csn,b1snl,b2snl,wps_pbk,wcs_pbk,phys)
    @inline bz34snl(s,n,l) = pbk_bz34snl(s,n,l,theta,b1snl,wps_pbk,wcs_pbk,phys)

    # BM functions
    @inline csnj(s,n,jj) = maxwell_csnj(s,n,jj,kz,J_opt,vtsz_bm,wcs_bm,us0_bm)
    # bsnj = (s,n,jj) -> maxwell_bsnj(s,n,jj,kz,J_opt,vts_parallel_bm,Ts_parallel_bm,Ts_perp_bm,wcs_bm)
    @inline b12snj(s,n,jj) = maxwell_b12snj(s,n,jj,kz,J_opt,vtsz_bm,Tsz_bm,Tsx_bm,sgms_bm,wcs_bm,lambdas_bm,EPS0)
    @inline b34snj(s,n,jj) = maxwell_b34snj(s,n,jj,kz,J_opt,vtsz_bm,Tsz_bm,Tsx_bm,sgms_bm,wcs_bm,lambdas_bm,EPS0)
    @inline b56snj(s,n,jj) = maxwell_b56snj(s,n,jj,kz,J_opt,vtsz_bm,Tsz_bm,Tsx_bm,sgms_bm,wcs_bm,lambdas_bm,EPS0)

    # x-component
    bx1_bm = maxwell_bx1(S_bm,Ns_bm,J,b12snj,csnj,wps_bm,phys)
    @inline bx1snj(s,n,jj) = maxwell_bx1snj(s,n,jj,b12snj,csnj,wps_bm,phys)
    bx2_bm = maxwell_bx2(S_bm,Ns_bm,J,b34snj,csnj,wps_bm,phys)
    @inline bx2snj(s,n,jj) = maxwell_bx2snj(s,n,jj,b34snj,csnj,wps_bm,phys)
    bx3_bm = maxwell_bx3(S_bm,Ns_bm,J,theta,b12snj,csnj,wps_bm,phys)
    @inline bx3snj(s,n,jj) = maxwell_bx3snj(s,n,jj,theta,b12snj,csnj,wps_bm,wcs_bm,phys)
    
    # y-component
    by1_bm = maxwell_by1(S_bm,Ns_bm,J,b34snj,csnj,wps_bm,phys)
    @inline by1snj(s,n,jj) = maxwell_by1snj(s,n,jj,b34snj,csnj,wps_bm,phys)
    by2_bm = maxwell_by2(S_bm,Ns_bm,J,b56snj,csnj,wps_bm,phys)
    @inline by2snj(s,n,jj) = maxwell_by2snj(s,n,jj,b56snj,csnj,wps_bm,phys)
    by3_bm = maxwell_by3(S_bm,Ns_bm,J,theta,b34snj,csnj,wps_bm,phys)
    @inline by3snj(s,n,jj) = maxwell_by3snj(s,n,jj,theta,b34snj,csnj,wps_bm,wcs_bm,phys)
    
    # z-component
    bz1_bm = maxwell_bz1(S_bm,Ns_bm,J,theta,b12snj,csnj,wps_bm,phys)
    @inline bz1snj(s,n,jj) = maxwell_bz1snj(s,n,jj,theta,b12snj,csnj,wps_bm,wcs_bm,phys)
    bz2_bm = maxwell_bz2(S_bm,Ns_bm,J,theta,b34snj,csnj,wps_bm,phys)
    @inline bz2snj(s,n,jj) = maxwell_bz2snj(s,n,jj,theta,b34snj,csnj,wps_bm,wcs_bm,phys)
    bz3_bm = maxwell_bz3(S_bm,Ns_bm,J,theta,b12snj,csnj,wps_bm,phys)
    @inline bz3snj(s,n,jj) = maxwell_bz3snj(s,n,jj,theta,b12snj,csnj,wps_bm,wcs_bm,phys)
    
 
    # MEMORY ALLOCATION FOR THE GLOBAL MATRIX
    if !cache.is_initialized
        # 1. Calc PBK size
        len_M_pbk = sum((2n + 1) * (k + 4) * (k + 1) for (n, k) in zip(Ns_pbk, kappasz_pbk)) ÷ 2 + 1
        # 2. Calc BM size
        len_M_bm = J * (2 * sum(Ns_bm) + length(Ns_bm)) + 1
        # 3. Total structure
        total_rows = 3*len_M_pbk + 3*len_M_bm + 9  
        total_cols = total_rows

        cache.M = zeros(ComplexF64, total_rows, total_cols)
        cache.is_initialized = true
    end 

    M = cache.M
    fill!(M, 0.0 + 0.0im)  # Reset the matrix to zero before filling it with new values.


    # ASSEMBLE THE GLOBAL MATRIX
    #######################
    # ExyzNo = 8 for Jx1_bm
    # ExyzNo = 7 for Jy1_bm
    # ExyzNo = 6 for Jz1_bm
    # ExyzNo = 5 for Ex
    # ExyzNo = 4 for Ey
    # ExyzNo = 3 for Ez
    # BxyzNo = 2 for Bx
    # BxyzNo = 1 for By
    # BxyzNo = 0 for Bz
    #####################
    current_row = 1  # To keep track of the current row index for filling the matrix

    # Step 1: Mx PBK (MatrixNo=1)
    h_x_pbk = Mxy_pbk_mixed!(M, current_row, S_pbk, Ns_pbk, Ns_bm, kappasz_pbk, J, csn,
        bx11snl, bx21snl, bx31snl, bx12snl, bx22snl, bx32snl,
        bx33snl, bx10_pbk, 1, 5, 4, 3, 5)
    current_row += h_x_pbk
    
    # My PBK (MatrixNo=2)
    h_y_pbk = Mxy_pbk_mixed!(M, current_row, S_pbk, Ns_pbk, Ns_bm, kappasz_pbk, J, csn,
        by11snl, by21snl, by31snl, by12snl, by22snl, by32snl,
        by33snl, by20_pbk, 2, 5, 4, 3, 4)
    current_row += h_y_pbk

    # Mz PBK (MatrixNo=3)
    h_z_pbk = Mz_pbk_mixed!(M, current_row, S_pbk, Ns_pbk, Ns_bm, kappasz_pbk, J, csn,
        bz11snl, bz21snl, bz31snl, bz13snl, bz23snl, bz33snl,
        bz12snl, bz22snl, bz32snl, bz34snl, by20_pbk, 3, 5, 4, 3)
    current_row += h_z_pbk

    #Step 2:  bi-Maxwellian plasma matrix
    # Mx BM (MatrixNo=4)
    h_x_bm = M_maxwell_mixed!(M, current_row, S_bm, Ns_pbk, Ns_bm, kappasz_pbk, J, csnj, bx1snj, bx2snj, bx3snj, 4, 5, 4, 3)
    current_row += h_x_bm

    # My BM (MatrixNo=5)
    h_y_bm = M_maxwell_mixed!(M, current_row, S_bm, Ns_pbk, Ns_bm, kappasz_pbk, J, csnj, by1snj, by2snj, by3snj, 5, 5, 4, 3)
    current_row += h_y_bm

    # Mz BM (MatrixNo=6)
    h_z_bm = M_maxwell_mixed!(M, current_row, S_bm, Ns_pbk, Ns_bm, kappasz_pbk, J, csnj, bz1snj, bz2snj, bz3snj, 6, 5, 4, 3)
    current_row += h_z_bm
   

    # Step 3: The perturbed currents of first term of BM
    # Incorporation the first-term of BM (Jx1_bm, Jy1_bm, Jz1_bm).
    # BM perturbations (Jx1)
    M[end-8,end-5] += bx1_bm  # the first term of dEx
    M[end-8,end-4] += bx2_bm  # the second term of dEy
    M[end-8,end-3] += bx3_bm # the third term of dEz

    # BM perturbations (Jy1)
    M[end-7,end-5] += by1_bm # the first term of dEx
    M[end-7,end-4] += by2_bm # the second term of dEy
    M[end-7,end-3] += by3_bm # the third term of dEz

    # BM perturbations (Jz1)
    M[end-6,end-5] += bz1_bm  # the first term of dEx
    M[end-6,end-4] += bz2_bm  # the second term of dEy
    M[end-6,end-3] += bz3_bm  # the third term of dEz

    
    # Step 4.1: The perturbed currents of Maxwell's equations for the PBK plasmas
    # Maxwell's equations for Ex                   
    idx_Jx_pbk = getIndexOfBlkMatrix_mixed(Ns_pbk,Ns_bm,kappasz_pbk,J,1)
    M[end-5,idx_Jx_pbk[end]] -= 1im/phys.epsilon0  # for the second item of Ex
    # Maxwell's equations for Ey
    idx_Jy_pbk = getIndexOfBlkMatrix_mixed(Ns_pbk,Ns_bm,kappasz_pbk,J,2)
    M[end-4,idx_Jy_pbk[end]] -= 1im/phys.epsilon0  # for the third item of Ey
    # Maxwell's equations Ez
    idx_Jz_pbk = getIndexOfBlkMatrix_mixed(Ns_pbk,Ns_bm,kappasz_pbk,J,3)
    M[end-3,idx_Jz_pbk[end]] -= 1im/phys.epsilon0  # for the second item of Ez

    # Step 4.2: The perturbed currents of Maxwell's equations for the bi-Maxwellian plasmas
    # Maxwell's equations for Ex
    idx_Jx_bm = getIndexOfBlkMatrix_mixed(Ns_pbk,Ns_bm,kappasz_pbk,J,4)
    M[end-5,end-1] += phys.c2*kz # for the first item of Ex
    M[end-5,end-8] -= 1im/phys.epsilon0 # for dJx1_bm
    M[end-5,idx_Jx_bm[1:end-1]] .-= 1im/phys.epsilon0 # for dJx2_bm

    # Maxwell's equations for Ey
    idx_Jy_bm = getIndexOfBlkMatrix_mixed(Ns_pbk,Ns_bm,kappasz_pbk,J,5)
    M[end-4,end-2] -= phys.c2*kz  # for the first item of Ey
    M[end-4,end] += phys.c2*kx  # for the second item of Ey
    M[end-4,end-7] -= 1im/phys.epsilon0 # for dJy1_bm
    M[end-4,idx_Jy_bm[1:end-1]] .-= 1im/phys.epsilon0 # for dJy2_bm

    # Maxwell's equations for Ez
    idx_Jz_bm = getIndexOfBlkMatrix_mixed(Ns_pbk,Ns_bm,kappasz_pbk,J,6)
    M[end-3,end-1] -= phys.c2*kx    # for the first item of Ez
    M[end-3,end-6] -= 1im/phys.epsilon0  # for dJz1_bm
    M[end-3,idx_Jz_bm[1:end-1]] .-= 1im/phys.epsilon0  # for dJz2_bm

    # Step 5: Maxwell's equations for the perturbed quantities of Bx, By and Bz.
    M[end-2,end-4] -= kz   # for Bx
    M[end-1,end-5] += kz   # for the first item of By
    M[end-1,end-3] -= kx   # for the second item of By
    M[end,end-4] += kx     # for Bz

    
    ##Solver##
    eigenvalues, eigenvectors = eigensolver(M, select_eigsolver)
    # eigenvalues, eigenvectors = eigensolver(M, select_eigsolver=1, nev=10, sort_by=:real_desc)
    # eigenvalues, eigenvectors = eigensolver(M, select_eigsolver=1, nev=10)

    return (eigenvalues, eigenvectors)
end
