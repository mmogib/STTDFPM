include("dependencies.jl")
include("included_files.jl")


options = AlgorithmOptions(2000, 1e-5)
"""
n: length of original signal
m: number of observations to take
r: number of nonzero spikes in the signal
"""
n = 2^5
m, r = fld(n, 4), fld(n, 8)
σ = 0.01
rng = MersenneTwister(20230515)
ATA, initial_point, c, original_signal, observed_signal = createCSData(m, n, r, σ)

flder_to_save = "./results/stored_data"
saveDataToFile(flder_to_save, ATA, initial_point, c, original_signal, observed_signal)

ATA, initial_point, c, original_signal, observed_signal = getSavedData(flder_to_save, m, n, r)
testProblems = createCSProblems(ATA, vec(initial_point), vec(c), n)

# [1] Pick parameters
# dim = 10_000
szs = 1000
itrs = 2000
# BBSTTDFPAlgorithmNoIN
params_dfs = runNumericalExperiments(
    BBSTTDFPParams, testProblems, itrs, szs;
    rng=rng,
    newparams=true
)











# [2] Smoothen Parameters

# testProblems = createTestProlems(10_000)
params_string = [
    # "t=0.48642251451762664, β=0.05197735523036173,  σ=0.8944837832002654, γ=1.0619048801778004, αmin=0.9205894528139231, αmax=1.45012474683393, r=0.983911897366891, ψ=0.45610249242743195, η=0.42930287193237837, ξ=0.6046090888298925",
    # "t=0.6803973835990211, β=0.9936779876442459,  σ=0.1043358242207244, γ=1.0619048801778004, αmin=0.9205894528139231, αmax=1.45012474683393, r=0.983911897366891, ψ=0.45610249242743195, η=0.08913779697247537, ξ=0.7262906361015746"
    # "t=0.48642251451762664, β=0.9687683369158444,  σ=0.5274431734206948, γ=1.0619048801778004, αmin=0.9205894528139231, αmax=1.45012474683393, r=0.983911897366891, ψ=0.45610249242743195, η=0.1670918140075346, ξ=0.8786033839504519"
    "t=0.39037961468610294, β=0.919013101031336,  σ=0.6305683958101485, γ=1.1397194144500726, αmin=0.007183651395935575, αmax=0.6550674547390956, r=0.4533115757253876, ψ=0.02310981356882147, η=0.3783809432766798, ξ=0.7512792447314791"
]
randomizing_args = (save=true,
    experiments_count=100,
    paramstr=params_string,
    randomized=[:σ],
    options=options,
    rng=rng,
    newparams=true
)
solving_args = (save=false,
    experiments_count=1,
    paramstr=params_string,
    randomized=nothing,
    options=options,
    rng=rng,
    newparams=true
)
smothened_params_dfs = runNumericalExperiments(BBSTTDFPParams, testProblems; solving_args...)
vscodedisplay.(smothened_params_dfs)