using Dates

"""
    display_runtime(runtime::Real)
    
Prints the total elapsed time.
Converts raw seconds into a structured duration (hours, minutes, seconds, milliseconds).
% Author: Bai Wei (baiwei12@mail.ustc.edu.cn,baiweiphys@gmail.com)
% Date: 2026.01.30
"""
function display_runtime(runtime::Real)
    # Convert seconds to milliseconds and create a Dates.Period object
    # We round to Int because Millisecond expects an integer value
    ms = Dates.Millisecond(round(Int, runtime * 1000))

    # canonicalize() automatically breaks down large units into smaller, 
    # readable ones (e.g., 3661000ms becomes 1 hour, 1 minute, 1 second)
    dur = Dates.canonicalize(ms)

    println("Elapsed time: $dur")
end
