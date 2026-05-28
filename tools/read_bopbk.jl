"""
    read_bopbk(path::AbstractString; delim::Char=' ')

read_bopbk(path::AbstractString; delim=' ') -> (df, par_data)
Read a bopbk.in file and return a DataFrame and its column-ordered matrix.
The file is expected to have a header row with column names.

Usage:
    include(joinpath(PROJECT_ROOT, "tools", "read_bopbk.jl"))
    df, par_data = read_bopbk("./bopbk.in")
    df.Tsz
    df.kappasz
"""

using CSV
using DataFrames

function read_bopbk(path::AbstractString; delim::Char=' ')
    df = CSV.read(path, DataFrame;
        delim=delim,
        header=true,
        ignorerepeated=true,
        types=Float64)
    par_data = Matrix(df)  # preserves column order
    return df, par_data
end
