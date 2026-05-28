using Plots
using Dates

"""
    export_plots(plt, base_name; dir="fig", formats=[".pdf", ".png", ".svg"], dpi=300)
"""
function export_plots(plt, base_name; dir="fig", formats=[".pdf", ".png", ".svg"], dpi=300)
    # Ensure the output directory exists
    mkpath(dir)
    pure_name = splitext(base_name)[1]

    println("--- Saving Plots ---")
    for ext in formats
        full_path = joinpath(dir, pure_name * ext)
        # default: dpi=300
        if ext == ".png"
            savefig(plt, full_path)
            # savefig(plt, full_path; dpi=dpi)
        else
            savefig(plt, full_path)
        end
        println("[✓] Saved: $full_path")
    end
end
