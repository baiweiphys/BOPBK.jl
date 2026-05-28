using HDF5, JLD2, MAT

"""
    save_results(prefix::String, dat_dir::String; formats=[:h5, :jld2, :mat], data...)

Save simulation results in multiple scientific data formats.
# Arguments
- `prefix`: String prefix for the filename.
- `dat_dir`: Destination directory path.
- `formats`: An array or tuple containing desired formats, e.g., `[:h5, :jld2, :mat]`.
- `data...`: Keyword arguments representing the physical variables to be stored.
"""
function save_results(prefix::String, dat_dir::String; formats=[:h5, :jld2, :mat], data...)
    # Ensure the output directory exists
    mkpath(dat_dir)

    # Preprocess data into a Dictionary (required for HDF5 and MAT formats).
    # Convert AbstractRange to Array (via collect) to ensure compatibility with HDF5.
    data_dict = Dict(string(k) => (v isa AbstractRange ? collect(v) : v) for (k, v) in data)

    # --- Save in HDF5 format: Standard for cross-platform scientific data ---
    if :h5 in formats
        h5_file = joinpath(dat_dir, "$(prefix)_data.h5")
        h5open(h5_file, "w") do file
            for (k, v) in data_dict
                write(file, k, v)
            end
        end
        println("[✓] HDF5 file saved: $h5_file")
    end

    # --- Save in JLD2 format: Native Julia format with compression ---
    if :jld2 in formats
        jld_file = joinpath(dat_dir, "$(prefix)_data.jld2")
        # Direct keyword argument passing for JLD2 serialization
        jldsave(jld_file; compress=true, data...)
        println("[✓] JLD2 file saved: $jld_file")
    end

    # --- Save in MAT format: Compatible with MATLAB environment ---
    if :mat in formats
        mat_file = joinpath(dat_dir, "$(prefix)_data.mat")
        matwrite(mat_file, data_dict)
        println("[✓] MAT file saved: $mat_file")
    end
end