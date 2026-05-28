"""
Description: # Plots for 3D plasma dispersion relation: frequency and growth rate contour plots
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
# kk = res["kk"]
# ww = res["ww"]
# wcs = res["wcs"]
# kk0 = res["kk0"]

@load "output/bopbk_data.jld2" theta_deg kk0 ww wcs


# -------------------------------
# Prepare data for contour plots
ww_contour = copy(ww)

idx1 = imag.(ww_contour ./ wcs[2]) .< 3.0e-3
idx2 = imag.(ww_contour ./ wcs[2]) .> 2.0e1

# Combine masks using element-wise logical OR
idx = idx1 .| idx2

# Assign NaN to filtered indices (NaNs are ignored by the plotting engine)
ww_contour[idx] .= NaN + im * NaN

# Extract real (frequency) and imaginary (growth rate) parts
wr_contour = real.(ww_contour)
wi_contour = imag.(ww_contour)


# Define shared plotting attributes
# plot_attr = (
#     fill=true,
#     linecolor=:transparent,
#     c=:jet,
#     xlabel=L"k \rho_{cp}",
#     ylabel=L"\theta^\circ",
#     lw=1,   # contour line width
#     colorbar_title="",
#     grid=false,
#     tickfont=font(12, "Arial"),
#     guidefont=font(14, "Arial")
# )

plot_attr = (
    fill=true,
    linecolor=:transparent,
    c=:jet,
    xlabel=L"k d_{i}",
    ylabel=L"\theta\ (\mathrm{^\circ})",
    lw=0.5,   # contour line width
    colorbar_title="",
    grid=false,
    tickfont=font(12, "Arial"),
    guidefont=font(14, "Arial"),
    colorbar_tickfont=font(11, "Arial"),
    framestyle=:box,
    right_margin = 12mm, 
    left_margin = 5mm
)

# --- Generate Subplots ---
num_start = 1
num_end = 200
#
# Subplot 1: Real Frequency (Normalized by cyclotron frequency)
p1 = contourf(kk0, theta_deg, wr_contour[:, :, num_start]' ./ abs(wcs[2]);
    title=L"\omega_r / \Omega_{cp}",
    levels=20,
    # Additional publication enhancements:
    tick_direction=:out,  # Outward ticks are cleaner
    minorticks=true,  # Minor ticks for better precision
    colorbar_position=:right,
    aspect_ratio=:auto,  # Maintain aspect ratio
    plot_attr...)

# Subplot 2: Growth Rate (Normalized by cyclotron frequency)
# p2 = contourf(kk0, theta_deg, wi_contour[:, :, num_start]' ./ abs(wcs[2]);
p2 = contourf(kk0, theta_deg, wi_contour[:, :, num_start]' ./ abs(wcs[2]);
    title=L"\omega_i / \Omega_{cp}",
    levels=30,
    # Additional publication enhancements:
    tick_direction=:out,  # Outward ticks are cleaner
    minorticks=true,  # Minor ticks for better precision
    colorbar_position=:right,
    aspect_ratio=:auto,  # Maintain aspect ratio
    plot_attr...)

# Combine subplots into a single figure
# plt = plot(p1, p2, layout=(2, 1), size=(600, 800))
plt = plot(p2, layout=(1, 1), size=(600, 400))

# Display the plot
display(plt)

# Save in multiple formats for different requirements
prefix::String = "plasmaDispersion_3D"
plot_dir = create_data_folder(pwd(), "fig")
export_plots(plt, prefix; dir=plot_dir, formats=[".pdf"])
# export_plots(plt, prefix, dir=plot_dir, formats=[".pdf", ".png", ".svg"], dpi=300)
