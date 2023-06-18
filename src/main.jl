include("dependencies.jl")
include("included_files.jl")

# [0] settings
# rng = MersenneTwister(202305150300)
# dims = [1000, 10_000, 100_000]
dims = [1000, 10_000, 100_000]
options = AlgorithmOptions(2000, 1e-11, Inf, (fx0, fx1) -> norm(fx0 - fx1) / norm(fx0) <= 1e-7)
saveit = true
algorithms = [
    (BBSTTDFPAlgorithm, true),
    # (BBSTTDFPAlgorithm, false),
    (BBSTTDFPAlgorithmNoIN, true),
    # (BBSTTDFPAlgorithmNoIN, false),
    (MOPCGAlgorithm, false),
    (CGDFPAlgorithm, false),
    (AHDFPMAlgorithm, false)
]






dfs = runExperiment(algorithms, dims, options; saveit=saveit)
# algorithmsname = [
#     "BBSTTDFPAlgorithm(old_params)",
#     "BBSTTDFPNoINAlgorithm(old_params)",
#     "CGDFPAlgorithm",
#     "MOPCGAlgorithm",
#     "AHDFPMAlgorithm"
# ]
# ndfs = map(algorithmsname) do algo
#     filter(df -> all(strip.(df[!, :algorithm]) .== "$algo"), dfs)
# end
# nndfs = vcat(ndfs...)
pyplot()
plts = plotResultsProfiles(algorithms, dfs)
p1 = plot(plts..., layout=(3, 1), size=(900, 800))
output_folder = outputfolder("./results/profiles")
map(["png", "svg", "pdf", "eps"]) do ext
    savefig(p1, "$output_folder/profile.$ext")
end
numericalResultDF2CSV(vcat(dfs...), "$output_folder/profile.csv")
# pgfplotsx()
# p2 = plot(plts..., layout=(3, 1), size=(800, 900), tex_output_standalone=true)
# savefig(p, "./results/latex/profiles.tex")
