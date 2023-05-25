include("dependencies.jl")
include("included_files.jl")
algorithms = [
    (BBSTTDFPAlgorithmNoIN, false),
    (BBSTTDFPAlgorithm, false),
    (MOPCGAlgorithm, false),
    (CGDFPAlgorithm, false),
    (AHDFPMAlgorithm, false)
]

input_file = "./results/profiles_2023_05_23_15_06.csv"
dfs = CSV.File(input_file) |> DataFrame
dfs = transform(dfs, :algorithm => ByRow(strip) => :algorithm)
algorithm_names = unique(dfs[!, :algorithm])
bbsttd_df = filter(row -> row[:algorithm] == "BBSTTDFPAlgorithm(old_params)", dfs)
bbsttd_no_int_df = filter(row -> row[:algorithm] == "BBSTTDFPNoINAlgorithm(old_params)", dfs)
mopcg_df = filter(row -> row[:algorithm] == "MOPCGAlgorithm", dfs)
cgdfp_df = filter(row -> row[:algorithm] == "CGDFPAlgorithm", dfs)
ahdfpm_df = filter(row -> row[:algorithm] == "AHDFPMAlgorithm", dfs)

pyplot()
plt1, plt2, plt3 = plotResultsProfiles(algorithms, [bbsttd_no_int_df, bbsttd_df, mopcg_df, cgdfp_df, ahdfpm_df])
plt2
out_folder = outputfolder("./results/profiles")
zdt = now(tz"Asia/Riyadh")
saving_file_name = "$out_folder/profiles_$(Dates.format(zdt,"yyyy_mm_dd_HH_MM"))"
map([(plt1, "Iterations"), (plt2, "Evalauations"), (plt3, "CPUTime")]) do (plt, label)
    savefig(plt, "$(saving_file_name)_$label.png")
    savefig(plt, "$(saving_file_name)_$label.svg")
    savefig(plt, "$(saving_file_name)_$label.pdf")
    savefig(plt, "$(saving_file_name)_$label.eps")
end
numericalResultDF2CSV(dfs, "$saving_file_name.csv")

