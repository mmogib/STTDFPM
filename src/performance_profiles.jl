include("dependencies.jl")
include("included_files.jl")

pyplot()
# input_file = "./results/paper/data/profiles_2023_05_23_15_06.csv"
input_file = "./results/paper/data/profiles_2023_06_07_12_00.csv"
plts = producePaperProfileImages(input_file)
