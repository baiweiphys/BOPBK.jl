# tools/CaseRunner.jl
module CaseRunner

using Revise

export run_case

"""
    run_case(name, path; main_file, plot_file, run_main, run_plot)

Arguments:
- `name`: Unique ID for the module (e.g., "Case01").
- `path`: Directory containing the scripts.
- `main_file`: Name of the calculation script.
- `plot_file`: Name of the plotting script.
- `run_main`: Boolean, set to false to skip the calculation.
- `run_plot`: Boolean, set to false to skip the plotting.
"""

"""
1. Example: skip plot_bopbk.jl
CaseRunner.run_case("Case01", "path/to/case", run_plot=false)

2. Example: run only plot_bopbk.jl
CaseRunner.run_case("Case01", "path/to/case", run_main=false)

3. Example: run both scripts (default)
CaseRunner.run_case("Case01", "path/to/case")
or,
CaseRunner.run_case("Case05", path_case, 
         main_file="main_bopbk.jl", 
         plot_file="plot_bopbk.jl");
"""
function run_case(name::String, path::String;
    main_file::String="main_bopbk.jl",
    plot_file::String="plot_bopbk.jl",
    run_main::Bool=true,   # New flag: toggle for main script
    run_plot::Bool=true    # New flag: toggle for plot script
)

    # Create the unique module name
    mod_name = Symbol("Module_", name)

    # Define the module in the Main (Jupyter global) scope
    @eval module $mod_name
    using Revise

    # This allows scripts to either use it or attempt to redefine it 
    # (though redefining const will still trigger a warning, not an error).
    # if !@isdefined(PROJECT_ROOT)
    #     PROJECT_ROOT = Main.PROJECT_ROOT
    # end

    # Inject local variables for this specific run
    target_path = $path
    main_script = $main_file
    plot_script = $plot_file
    do_main = $run_main
    do_plot = $run_plot

    println(">>> STARTING CASE: ", $name)

    cd(target_path) do
        # 1. Logic for Calculation Script
        if do_main
            if isfile(main_script)
                @info "Running $main_script ..."
                includet(main_script)
            else
                @error "Main script not found: $main_script"
            end
        else
            @info "Skipping main script as requested."
        end

        # 2. Logic for Plotting Script
        if do_plot
            if isfile(plot_script)
                @info "Running $plot_script ..."
                includet(plot_script)
            else
                @error "Plot script not found: $plot_script"
            end
        else
            @info "Skipping plot script as requested."
        end
    end

    println(">>> FINISHED CASE: ", $name, "\n")
    end

    # Return nothing here so Jupyter doesn't display the Module object
    return nothing
end

end # module CaseRunner