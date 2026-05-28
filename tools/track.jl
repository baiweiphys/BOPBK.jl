# =============================================================================
# track.jl - Dispersion Branch Tracking and Mode Selection
# =============================================================================

using DataInterpolations
using LinearAlgebra

# -----------------------------------------------------------------------------
# 1. Struct Definitions (BranchPoint, SeedLocators)
# -----------------------------------------------------------------------------

abstract type AbstractSeedLocator end

struct SeedByReal <: AbstractSeedLocator end
struct SeedByImag <: AbstractSeedLocator end
struct SeedByAbs <: AbstractSeedLocator end

(::SeedByReal)(xs, x0) = argmin(x -> abs(real(x) - real(x0)), xs)
(::SeedByImag)(xs, x0) = argmin(x -> abs(imag(x) - imag(x0)), xs)
(::SeedByAbs)(xs, x0) = argmin(x -> abs(x - x0), xs)

# Select the best locator strategy based on the initial seed frequency
_default_seed_locator(ω) =
    iszero(imag(ω)) ? SeedByReal() : (iszero(real(ω)) ? SeedByImag() : SeedByAbs())

"""
    BranchPoint{T, W, L}

Initial point specification for 1D dispersion branch tracking (k-scan).
"""
struct BranchPoint{T,W,L}
    k::T
    ω::W
    locator::L
end

BranchPoint(k, ω) = BranchPoint(k, ω, _default_seed_locator(ω))

"""
    SurfaceBranchPoint{T, W, L}

Initial point specification for 2D dispersion surface tracking (k-theta scan).
"""
struct SurfaceBranchPoint{T,W,L}
    k::T
    θ::T
    ω::W
    locator::L
end

SurfaceBranchPoint(k, θ, ω) =
    SurfaceBranchPoint(k, θ, ω, _default_seed_locator(ω))


# -----------------------------------------------------------------------------
# 2. Interpolation Helper
# -----------------------------------------------------------------------------

"""
    interpolate_complex(z_prev, x_prev, x_new)

Predict the next complex frequency using independent PCHIP interpolation 
on real and imaginary parts. Enables extrapolation for branch prediction.
"""
function interpolate_complex(z_prev, x_prev, x_new)
    # Use simple linear extrapolation if only 1 or 2 points are available
    if length(z_prev) == 1
        return z_prev[1]
    elseif length(z_prev) == 2
        z1, z2 = z_prev
        x1, x2 = x_prev
        return z1 + (z2 - z1) * (x_new - x1) / (x2 - x1)
    else
        # PCHIP interpolation (Monotonicity preserving)
        # Note: Extrapolation is explicitly enabled via Extension

        # Real part interpolation
        itp_real = PCHIPInterpolation(real.(z_prev), x_prev; extrapolation=ExtrapolationType.Extension)
        xr = itp_real(x_new)

        # Imaginary part interpolation
        itp_imag = PCHIPInterpolation(imag.(z_prev), x_prev; extrapolation=ExtrapolationType.Extension)
        xi = itp_imag(x_new)

        return complex(xr, xi)
    end
end


# -----------------------------------------------------------------------------
# 3. Core Tracking Logic (track)
# -----------------------------------------------------------------------------

track(sol, point::NTuple{2,Any}) = track(sol, BranchPoint(point...))
track(sol, point::NTuple{3,Any}) = track(sol, SurfaceBranchPoint(point...))

"""
    track(solution, point)

Track a single dispersion branch across parameter space starting from `point`.
"""
track(sol, point::BranchPoint) = _track(sol.ks, sol.ωs, point.k, point.ω)

# Internal 1D tracking implementation
function _track(ks, ωs_all, k_start, ω_seed)
    nk = length(ks)

    # Find nearest k index to the starting point
    start_idx = argmin(abs.(ks .- k_start))

    # Initialize result array with NaNs
    ω_branch = Vector{ComplexF64}(undef, nk)
    fill!(ω_branch, NaN + NaN * im)

    # Find the closest match to the seed in the root pool at the start index
    current_roots = ωs_all[start_idx]
    if isempty(current_roots)
        @warn "No roots at start index $start_idx"
        return ks, ω_branch
    end
    best_idx = argmin(abs.(current_roots .- ω_seed))
    ω_branch[start_idx] = current_roots[best_idx]

    # Track Forward in k
    track_from_index!(ω_branch, ks, ωs_all, start_idx, (start_idx+1):nk)

    # Track Backward in k
    track_from_index!(ω_branch, ks, ωs_all, start_idx, (start_idx-1):-1:1)

    return ks, ω_branch
end

# Directional tracking helper that updates ω_branch in-place
function track_from_index!(ω_branch, ks, results, start_idx, indices)
    isempty(indices) && return

    is_forward = first(indices) > start_idx

    for i in indices
        # Get up to 3 previous points for interpolation/prediction
        prev_candidates = is_forward ? ((i-3):(i-1)) : ((i+1):(i+3))

        valid_prev = filter(idx -> idx >= 1 && idx <= length(ks) && !isnan(real(ω_branch[idx])), prev_candidates)

        if isempty(valid_prev)
            # Fallback to starting point if a gap is encountered
            valid_prev = [start_idx]
        end

        # Sort by x value for interpolation consistency
        sort!(valid_prev, by=idx -> ks[idx])

        x_prev = ks[valid_prev]
        y_prev = ω_branch[valid_prev]

        # Predict the next frequency value
        ω_pred = interpolate_complex(y_prev, x_prev, ks[i])

        # Match: Find the nearest eigenvalue in the pool at current index i
        pool = results[i]
        if isempty(pool)
            ω_branch[i] = NaN + NaN * im
        else
            best_match = pool[argmin(abs.(pool .- ω_pred))]
            ω_branch[i] = best_match
        end
    end
end

# 2D Tracking implementation (Surface tracking)
function track(sol, point::SurfaceBranchPoint)
    ks, θs = sol.ks, sol.θs
    nk, nθ = length(ks), length(θs)

    # Find starting coordinate indices
    idx_k0 = argmin(abs.(ks .- point.k))
    idx_θ0 = argmin(abs.(θs .- point.θ))

    ω_surface = Matrix{ComplexF64}(undef, nk, nθ)
    fill!(ω_surface, NaN + NaN * im)

    # Locate initial seed root
    roots_pool = sol.ωs[idx_k0, idx_θ0]
    if isempty(roots_pool)
        @warn "No roots at start indices ($idx_k0, $idx_θ0)"
        return ks, θs, ω_surface
    end
    ω_seed = point.locator(roots_pool, point.ω)

    # 1. Fill the initial Theta slice (Track along k for the fixed θ0)
    ωs_at_theta0 = sol.ωs[:, idx_θ0]
    _, ω_slice = _track(ks, ωs_at_theta0, point.k, ω_seed)
    ω_surface[:, idx_θ0] .= ω_slice

    # 2. Track in Theta direction for every k using the initial slice as seed

    # Forward along Theta
    for j in (idx_θ0+1):nθ
        for i in 1:nk
            ω_prev = ω_surface[i, j-1]
            isnan(real(ω_prev)) && continue

            # Simple zero-order prediction (assuming continuity in theta)
            ω_pred = ω_prev

            pool = sol.ωs[i, j]
            if !isempty(pool)
                ω_surface[i, j] = pool[argmin(abs.(pool .- ω_pred))]
            end
        end
    end

    # Backward along Theta
    for j in (idx_θ0-1):-1:1
        for i in 1:nk
            ω_prev = ω_surface[i, j+1]
            isnan(real(ω_prev)) && continue

            ω_pred = ω_prev
            pool = sol.ωs[i, j]
            if !isempty(pool)
                ω_surface[i, j] = pool[argmin(abs.(pool .- ω_pred))]
            end
        end
    end

    return ks, θs, ω_surface
end


# -----------------------------------------------------------------------------
# 4. Mode Selection and Auto-Tracking Interface
# -----------------------------------------------------------------------------

abstract type ModeSelector end

"""
    SelectMostUnstable(n::Int)

Strategy to find the top `n` modes with the largest imaginary part (highest growth rate).
"""
struct SelectMostUnstable <: ModeSelector
    n::Int
end

"""
    SelectLeastDamped(n::Int)

Strategy to find the top `n` modes with the smallest absolute imaginary part 
(closest to the real axis, i.e., least damping or slowest growth).
"""
struct SelectLeastDamped <: ModeSelector
    n::Int
end

"""
    SelectClosest(target::Number, n::Int=1)

Strategy to find the `n` modes geometrically closest to a target complex frequency.
"""
struct SelectClosest{T<:Number} <: ModeSelector
    target::T
    n::Int
end

"""
    SelectRange(r_min, r_max, i_min, i_max)

Strategy to select all modes within a rectangular box in the complex plane.
"""
struct SelectRange{T<:Real} <: ModeSelector
    r_min::T
    r_max::T
    i_min::T
    i_max::T
end

# --- track_modes Dispatch Wrappers ---

"""
    track_modes(sol; k_idx=1, theta_idx=1)

Default behavior: Track the single most unstable mode.
"""
function track_modes(sol; kwargs...)
    track_modes(sol, SelectMostUnstable(1); kwargs...)
end

"""
    track_modes(sol, selector::ModeSelector; k_idx=1, theta_idx=1)

Track branches identified by the given `selector` strategy starting from `k_idx`.
"""
function track_modes(sol, selector::ModeSelector; k_idx::Int=1, theta_idx::Int=1)
    # 1. Extract initial pool of roots at starting indices
    roots_pool = if sol.ωs isa AbstractMatrix # 2D case
        sol.ωs[k_idx, theta_idx]
    else # 1D case
        sol.ωs[k_idx]
    end

    if isempty(roots_pool)
        @warn "No roots found at initial index (k_idx=$k_idx)"
        return []
    end

    # 2. Select specific roots based on the strategy logic
    selected_roots = select_roots(selector, roots_pool)

    # 3. Track each selected root branch
    branches = []
    is_2d = sol.ωs isa AbstractMatrix

    for root in selected_roots
        if is_2d
            k_val = sol.ks[k_idx]
            th_val = sol.θs[theta_idx]
            bp = SurfaceBranchPoint(k_val, th_val, root)
            push!(branches, track(sol, bp))
        else
            k_val = sol.ks[k_idx]
            bp = BranchPoint(k_val, root)
            push!(branches, track(sol, bp))
        end
    end

    return branches
end

# --- Selection Logic Implementation ---

function select_roots(sel::SelectMostUnstable, roots::AbstractVector{<:Number})
    # Sort by imaginary part descending (Highest growth rate first)
    perm = sortperm(imag.(roots), rev=true)
    n_pick = min(length(roots), sel.n)
    return roots[perm[1:n_pick]]
end

function select_roots(sel::SelectLeastDamped, roots::AbstractVector{<:Number})
    # Sort by absolute value of imaginary part ascending (Closest to real axis first)
    perm = sortperm(abs.(imag.(roots)))
    n_pick = min(length(roots), sel.n)
    return roots[perm[1:n_pick]]
end

function select_roots(sel::SelectClosest, roots::AbstractVector{<:Number})
    # Sort by distance in the complex plane
    dists = abs.(roots .- sel.target)
    perm = sortperm(dists)
    n_pick = min(length(roots), sel.n)
    return roots[perm[1:n_pick]]
end

function select_roots(sel::SelectRange, roots::AbstractVector{<:Number})
    # Filter roots within the specified real and imaginary bounds
    filter(roots) do w
        r, i = real(w), imag(w)
        r >= sel.r_min && r <= sel.r_max && i >= sel.i_min && i <= sel.i_max
    end
end