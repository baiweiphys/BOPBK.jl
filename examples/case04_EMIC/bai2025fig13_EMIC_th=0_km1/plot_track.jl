using JLD2
using LaTeXStrings
using Plots
using Measures

const PROJECT_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "tools", "create_data_folder.jl"))
include(joinpath(PROJECT_ROOT, "tools", "export_plots.jl"))
include(joinpath(PROJECT_ROOT, "tools", "track.jl"))

# ==============================================================================
# 0. Helper Function: Global Search with Specific Rank Selection
# ==============================================================================
"""
    find_global_ranked_branches(sol; ranks=[1], strategy=:unstable)

1. Scans all k-space for candidates.
2. Sorts them by the chosen strategy.
3. Tracks and de-duplicates to identify unique physical branches.
4. Returns only the branches specified in the `ranks` list (e.g., [1, 4, 7]).
"""
function find_global_ranked_branches(sol; ranks=[1], strategy=:unstable)
    # 1. Collect all potential starting points (imag, k_idx, complex_w)
    candidates = []
    for (ik, roots) in enumerate(sol.ωs)
        for w in roots
            push!(candidates, (imag(w), ik, w))
        end
    end

    # 2. Sort candidates globally based on strategy
    if strategy == :unstable
        sort!(candidates, by=x -> x[1], rev=true) # Descending (Max Gamma first)
        println("Strategy: Most Unstable (Imag descending)")
    elseif strategy == :least_damped
        sort!(candidates, by=x -> abs(x[1]), rev=false) # Ascending (Closest to 0 first)
        println("Strategy: Least Damped (abs(Imag) ascending)")
    else
        error("Unknown strategy: $strategy")
    end

    unique_branches = []
    tracked_hashes = Set{UInt64}()
    max_rank_needed = maximum(ranks)

    println("Scanning and identifying unique branches...")

    # 3. Track candidates until we find enough unique branches
    for cand in candidates
        if length(unique_branches) >= max_rank_needed
            break
        end

        gamma_val, k_idx, w_target = cand

        # Use SelectClosest to lock onto the specific root at k_idx
        track_results = track_modes(sol, SelectClosest(w_target, 1); k_idx=k_idx)
        if isempty(track_results)
            continue
        end
        new_branch = track_results[1]

        # De-duplication: Use a hash of the real frequency at the middle of the branch
        # This prevents picking the same physical mode from different k-starting points
        mid_idx = length(new_branch[1]) ÷ 2
        id_hash = hash(round(real(new_branch[2][mid_idx]), digits=3))

        if !(id_hash in tracked_hashes)
            push!(unique_branches, (branch=new_branch, start_k=sol.ks[k_idx], start_w=w_target))
            push!(tracked_hashes, id_hash)
            current_rank = length(unique_branches)
            println("  -> Found Unique Branch #$current_rank at k=$(round(sol.ks[k_idx], digits=3))")
        end
    end

    # 4. Pick only the ranks requested by the user
    final_selection = []
    for r in ranks
        if r <= length(unique_branches)
            # Store the rank number along with the branch for labeling
            push!(final_selection, (rank=r, data=unique_branches[r]))
        else
            @warn "Rank $r requested, but only $(length(unique_branches)) unique branches were found."
        end
    end

    return final_selection
end

# ==============================================================================
# 1. Main Script
# ==============================================================================

# Load Data
@load "output/bopbk_data.jld2" theta_deg kk0 ww wcs

# Normalized Solution object
norm_val = abs(wcs[1])
ww_vec = [ww[i, :] ./ norm_val for i in 1:size(ww, 1)]
sol = (ks=kk0, θs=[theta_deg], ωs=ww_vec)

# --- USER CONFIGURATION ---
# target_strategy: :unstable or :least_damped
# target_ranks: e.g., [1] for only the best, [1, 3, 5] for specific ones
target_strategy = :unstable
# target_strategy = :least_damped
target_ranks = [1,3] # Example: top 3 branches

branches_to_plot = find_global_ranked_branches(sol; ranks=target_ranks, strategy=target_strategy)

println("Tracking complete. Plotting...")

# ==============================================================================
# 2. Plotting
# ==============================================================================
gr()

title_str = target_strategy == :unstable ? "Most Unstable" : "Least Damped"

p_real = plot(ylabel=L"\omega_r/\omega_{pe}",
    grid=true,
    ylims=(0, 1),
    legend=:topleft,
    title="$title_str Selection (Ranks: $target_ranks)",)

p_imag = plot(xlabel=L"ck/\omega_{pe}",
    ylabel=L"\gamma/\omega_{pe}",
    grid=true,
    ylims=(-0.1, 0.05),
    legend=:false,
    # legend=:outerright
    )

# 2.1 Background (All roots in grey)
scatter!(p_real, kk0, real.(ww) ./ norm_val,
    color=:grey,
    alpha=0.1,
    markersize=1.5,
    label="",
    markerstrokewidth=0)

scatter!(p_imag, kk0, imag.(ww) ./ norm_val,
    color=:grey,
    alpha=0.1,
    markersize=1.5,
    label="",
    markerstrokewidth=0)

# 2.2 Highlighted Ranked Branches
colors = [:red, :blue, :green, :purple, :orange, :magenta, :cyan]

for (i, item) in enumerate(branches_to_plot)
    rank_num = item.rank
    branch_data = item.data.branch

    ks_track = branch_data[1]
    ws_track = branch_data[2]

    # wr = real.(ws_track)
    wr = abs.(real.(ws_track))
    wi = imag.(ws_track)

    col = colors[mod1(i, length(colors))]

    # Plot curves
    plot!(p_real, ks_track, wr, lw=2, color=col, label="Rank $rank_num")
    plot!(p_imag, ks_track, wi, lw=2, color=col, label="Rank $rank_num")

    # Mark the optimal point for this branch
    if target_strategy == :unstable
        val, idx = findmax(wi)
    else
        val_abs, idx = findmin(abs.(wi))
        val = wi[idx]
    end
    scatter!(p_imag, [ks_track[idx]], [val], marker=:star5, markersize=5, color=col, label="")

end

# Combine
plt = plot(p_real, p_imag, layout=(2, 1), size=(800, 800), margin=5Plots.mm)
display(plt)


# Save in multiple formats for different requirements
prefix::String = "ranked_$(target_strategy)_trace"
plot_dir = create_data_folder(pwd(), "fig")
export_plots(plt, prefix; dir=plot_dir, formats=[".pdf"])
# export_plots(plt, prefix, dir=plot_dir, formats=[".pdf", ".png", ".svg"], dpi=300)
