"""
Description: Top-level module entry for BOPBK.
Collects solvers and exposes a helper to select the appropriate solver
based on plasma parameter flags.
Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
Created: 2026-02-01
Changelog:
2026-02-08: Optimized select_solver, add SolverCache. 
"""

module BOPBK

# Core physical constants
include("PhysicalConstants.jl")
using .PhysicalConstants: PlasmaConfig

# Plasma parameters
include("PlasmaParameters.jl")


# Utils
const utils_path = abspath(joinpath(@__DIR__, "utils"))
include(joinpath(utils_path, "eigensolver.jl"))


# Solvers
const solvers_path = abspath(joinpath(@__DIR__, "solvers"))
include(joinpath(solvers_path, "solver_bm.jl"))
include(joinpath(solvers_path, "solver_pbk.jl"))
include(joinpath(solvers_path, "solver_mixed.jl"))

export solver_bm!, solver_pbk!, solver_mixed!, eigensolver, PlasmaParameters


"""
SolverCache
Holdds pre-allocated memory to reuse across solver iterations.
"is_initialized" flag allows for lazy allocation on the first call.
"""
mutable struct SolverCache
    M::Union{Matrix{ComplexF64}}
    is_initialized::Bool

    function SolverCache()
        #Initialize with an empty matrix, will be resized on first run 
        new(Matrix{ComplexF64}(undef, 0, 0), false)
    end
end



"""
select_solver(plasmaParams, select_eigsolver, EPS0, phys) -> (kx,kz,theta)->(eigs, vecs)
Return a solver closure based on `plasmaParams.flags`:
- PBK only: `solver_pbk`
- BM only:  `solver_bm`
- Mixed:    `solver_mixed`
"""
function select_solver(
    plasmaParams,
    select_eigsolver::Int,
    EPS0::Float64,
    phys::PlasmaConfig)

    (; is_pbk, is_bm) = plasmaParams.flags

    # Create a single SolverCache instance for this solver loop
    cache = SolverCache()

    if count(is_bm) == 0
        return (kx, kz, theta_rad) -> solver_pbk!(cache, kx, kz, theta_rad, plasmaParams, select_eigsolver, EPS0, phys)
    elseif count(is_pbk) == 0
        return (kx, kz, theta_rad) -> solver_bm!(cache, kx, kz, theta_rad, plasmaParams, select_eigsolver, EPS0, phys)
    else
        return (kx, kz, theta_rad) -> solver_mixed!(cache, kx, kz, theta_rad, plasmaParams, select_eigsolver, EPS0, phys)
    end

end # select_solver

end # module BOPBK