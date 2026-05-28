"""
    create_data_folder(path, name)
    
Description: Create a new data folder if it doesn't exist.
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2023-10-01
LastEditors: Bai Wei
Last Modified: 2025-12-10
"""
function create_data_folder(path, name)
    dir = joinpath(path, name)
    isdir(dir) || mkdir(dir)
    println(isdir(dir) ? "Folder ready: $dir" : "Created: $dir")
    return dir
end