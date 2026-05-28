"""
Description: Plots for plasma dispersion relation: frequency and growth rate plots.
Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
Date: 2026-01-31
Last Modified: 2026-02-01
"""

using JLD2  # for .jld2 files
using LaTeXStrings
using Plots
using Measures

const PROJECT_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "tools", "create_data_folder.jl"))
include(joinpath(PROJECT_ROOT, "tools", "export_plots.jl"))

# Define the path to the generated data file
# dat_file = "output/bopbk_data.jld2"

# Load the entire dataset into a Dictionary
# res = load(dat_file)

# Extract specific variables from the Dictionary
# theta_deg = res["theta_deg"]
# ww = res["ww"]
# kk0 = res["kk0"]
# wcs = res["wcs"]
# wn = res["wn"]

@load "output/bopbk_data.jld2" kk0 ww wn




# -------------------------------
# Prepare data for contour plots
wr = real(ww)
wi = imag(ww)

# Create and display the plot directly
firstNo = 1
lastNo = 500
rootsNo = firstNo:lastNo
plt = plot(
    # First subplot: Real frequency (ω_r/ω_pe)
    scatter(kk0, wr[:, rootsNo] / wn,
        color=:black, markersize=2,
        xlabel=L"k\lambda_{De}", ylabel=L"\omega_r/\omega_{pe}",
        xlims=(0, 0.37), 
        ylims=(0, 3),
        grid=true, legend=false),

    # Second subplot: Growth rate (γ/ω_pe)
    scatter(kk0, wi[:, rootsNo] / wn,
        color=:black, markersize=2,
        xlabel=L"k\lambda_{De}", ylabel=L"\gamma/\omega_{pe}",
        xlims=(0, 0.37), ylims=(-1.0, 0.25),
        grid=true, legend=false),
        #
        layout=(2, 1),     # 2 rows, 1 column layout
        size=(600, 800)    # Figure size: width=600px, height=800px
)

# Display the plot
display(plt)

# Save in multiple formats for different requirements
prefix::String = "dispersion_relation"
plot_dir = create_data_folder(pwd(), "fig")
export_plots(plt, prefix; dir=plot_dir, formats=[".pdf"])
# export_plots(plt, prefix, dir=plot_dir, formats=[".pdf", ".png", ".svg"], dpi=300)
