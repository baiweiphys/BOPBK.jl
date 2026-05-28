# solve_eigenproblem.jl
"""
Description: Solve eigenvalue problem and return eigenvalues and eigenvectors sorted by imaginary part
Parameters:
M: Input matrix (square matrix)
sp: 0 for dense (compute all), 1 for sparse (compute partial)
nev: Number of eigenvalues to compute (sparse case only)
sort_by: :imag_desc (default, descending by imaginary part), :imag_asc, :real_desc, :real_asc, :abs_desc, :abs_asc
Returns:
eigenvalues: Array of eigenvalues sorted by imaginary part (descending by default)
eigenvectors: Matrix of corresponding eigenvectors (columns are eigenvectors)
#
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2025-12-10
Last Modified: 2025-12-10
"""


# using MKL
using Arpack
using SparseArrays
using LinearAlgebra


function eigensolver(
    M::AbstractMatrix,
    select_eigsolver::Int=0,
    nev::Int=size(M, 2),
    sort_by::Symbol=:imag_desc)

    # Check if matrix is square
    n, m = size(M)
    if n != m
        error("Input matrix must be square, got size ($n, $m)")
    end

    if select_eigsolver == 0
        # Dense matrix: compute all eigenvalues and eigenvectors
        F = eigen(M)
        eigenvalues = F.values
        eigenvectors = F.vectors

    elseif select_eigsolver == 1
        # Sparse matrix: compute partial eigenvalues
        nev = min(max(nev, 1), n - 1)  # Arpack requires nev < n

        # Convert to sparse matrix
        M_sparse = sparse(M)
        # M_sparse = sparse(ComplexF64.(M))

        # Compute eigenvalues with largest imaginary parts
        # which=:LI gives eigenvalues with largest imaginary parts
        eigenvalues, eigenvectors = eigs(M_sparse, nev=nev - 2, which=:LI, maxiter=3000)
        # eigenvalues, eigenvectors = eigs(M_sparse, nev=nev - 2, which=:LI, maxiter=3000)
    else
        error("select_eigsolver parameter must be 0 (eigen) or 1 (eigs)")
    end

    # Sort eigenvalues and eigenvectors according to specified criterion
    if sort_by == :imag_desc
        # Sort by imaginary part in descending order (default)
        idx = sortperm(imag.(eigenvalues), rev=true)

    elseif sort_by == :imag_asc
        # Sort by imaginary part in ascending order
        idx = sortperm(imag.(eigenvalues))

    elseif sort_by == :real_desc
        # Sort by real part in descending order
        idx = sortperm(real.(eigenvalues), rev=true)

    elseif sort_by == :real_asc
        # Sort by real part in ascending order
        idx = sortperm(real.(eigenvalues))

    elseif sort_by == :abs_desc
        # Sort by magnitude in descending order
        idx = sortperm(abs.(eigenvalues), rev=true)

    elseif sort_by == :abs_asc
        # Sort by magnitude in ascending order
        idx = sortperm(abs.(eigenvalues))

    else
        # Default: no sorting
        idx = 1:length(eigenvalues)
        # @info "No sorting applied"
    end

    # Apply sorting
    eigenvalues_sorted = eigenvalues[idx]
    eigenvectors_sorted = eigenvectors[:, idx]

    # eigenvalues_sorted = ComplexF64.(eigenvalues[idx])
    # eigenvectors_sorted = ComplexF64.(eigenvectors[:, idx])

    # println("vals type: ", typeof(eigenvalues_sorted))
    # println("vals size: ", size(eigenvalues_sorted))
    # println("vecs type: ", typeof(eigenvectors_sorted))
    # println("vecs size: ", size(eigenvectors_sorted))

    return eigenvalues_sorted, eigenvectors_sorted
end
