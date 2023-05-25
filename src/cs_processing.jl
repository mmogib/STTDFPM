include("dependencies.jl")
include("included_files.jl")


# cameraman = testimage("cameraman.tif")
# barbara = testimage("barbara_gray_512.bmp")
# lena = testimage("lena_gray_512.tif")
# chart = testimage("resolution_test_512.tif")
# x = vec(cameraman);
# rng = MersenneTwister(20230511)
# # n = length(x)
# # m = 100
# # r = 64
# # A = rand(rng, m, n)
# # w = rand(rng, m)
# # b = A * x + w |> y->Float64.(y)
# # ATA = A'*A
# m, n, r = 512, 2048, 64
# x = zeros(n)
# q = randperm(n)
# A = randn(rng, m, n)
# B = qr(A).Q * Matrix(I, size(A)...)
# b = A * x + vec(0.001 * rand(rng, m, 1)) #|> y->Float64.(x)
# plot(b)
# denoise(b) |> plot

# # img = load("./imgs/lenna.png")
# # ndims(img)
# # x = permutedims(Float64.(lena))
# # L = 2
# # xts = wplotim(x, L, wavelet(WT.db3, WT.Filter))
# # Gray.(xts)

"""
Solve 
        min ||Ax-b||² + λ||x||₁  or   min ||Ax-b||² + τ||x||₁ 
"""

# rng = MersenneTwister(20230515)
# params1 = BBSTTDFPParams(
#     0.19091055694541392,
#     0.8248500773070584,
#     0.6059149417443106,
#     1.0160210808836148,
#     0.968177563496758,
#     1.8227259069163115,
#     0.9507650425301077,
#     0.8808475746290094,
#     0.030826250174403214,
#     0.3719729951071617
# )
# params0 = BBSTTDFPParams(
#     0.98,
#     0.8,
#     0.0001,
#     1.8,
#     1e-30,
#     Inf,
#     0.1,
#     0.5,
#     nothing,#0.001,
#     nothing#0.6
# )
do_expriment(;
    k=7, σ=0.01, inertial=true,
    options=AlgorithmOptions(5000, 1e-12),
    params=BBSTTDFPParams(
        0.98,
        0.9,
        0.0001,
        1.8,
        1e-30,
        Inf,
        0.1,
        0.5,
        nothing,#0.001,
        nothing#0.6
    )) = begin

    n = 2^k
    m, r = fld(n, 4), fld(n, 8)
    # flder_to_save = "./results/stored_data"
    # ATA, initial_point, c, original_signal, observed_signal = getSavedData(flder_to_save, m, n, r)
    ATA, initial_point, c, original_signal, observed_signal = createCSData(m, n, r, σ)
    problems = createCSProblems(ATA, vec(initial_point), vec(c), n)
    problem = problems[1]
    sol_bbs = BBSTTDFP(problem; options, params, inertial)
    z_bbs = sol_bbs.result.solution
    reconstructed_bbs = z_bbs[1:n] - z_bbs[n+1:end]
    sol_ahd = AHDFPM(problem; options=options, params=AHDFPMParameters(0.6, 1.75, 0.0001, 0.001, 0.4, 2.4, 0.8, 0.1, 1))
    z_ahd = sol_ahd.result.solution
    reconstructed_ahd = z_ahd[1:n] - z_ahd[n+1:end]

    sol_cgd = CGDFP(problem; options=options, params=CGDFPParameters(0.6, 1, 0.001, 0.7, 0.3))
    z_cgd = sol_cgd.result.solution
    reconstructed_cgd = z_cgd[1:n] - z_cgd[n+1:end]

    sol_mopcg = MOPCG(problem; options=options, params=MOPCGParameters(1, 0.6, 0.1, 0.2))

    z_mopcg = sol_mopcg.result.solution
    reconstructed_mopcg = z_mopcg[1:n] - z_mopcg[n+1:end]

    solution = (Dict(
            :name => "DATA",
            :dims => (n, m, r),
            :ATA => ATA,
            :c => c,
            :problem => problem,
            :initial => initial_point,
            :original => original_signal,
            :observed => observed_signal,
        ),
        Dict(
            :name => "BBSTTDFP",
            :reconstructed => reconstructed_bbs,
            :fullsolution => sol_bbs
        ),
        Dict(
            :name => "AHDFPM",
            :reconstructed => reconstructed_ahd,
            :fullsolution => sol_ahd
        ),
        Dict(
            :name => "CGDFP",
            :reconstructed => reconstructed_cgd,
            :fullsolution => sol_cgd
        ),
        Dict(
            :name => "MOPCG",
            :reconstructed => reconstructed_mopcg,
            :fullsolution => sol_mopcg
        )
    )
    solution
end
no_of_exprs = 100
solutions = map(1:no_of_exprs) do _
    # nparams = getNewParams(BBSTTDFPParams(
    #         0.9,
    #         0.85,
    #         0.9,
    #         1.75,
    #         1e-30,
    #         Inf,
    #         0.94,
    #         0.1,
    #         nothing,#0.001,
    #         nothing#0.6
    #     ), [:σ, :r, :ψ])
    nparams = BBSTTDFPParams(
        0.9,
        0.85,
        0.9,
        1.75,
        1e-30,
        Inf,
        0.94,
        0.1,
        nothing,#0.001,
        nothing#0.6
    )
    do_expriment(k=11,
        inertial=true,
        params=nparams,
        options=AlgorithmOptions(2_000, 1e-7, Inf, (fx0, fx1) -> (norm(fx0 - fx1) / norm(fx0)) <= 1e-7),
    )
end
msis = map(1:no_of_exprs) do i
    solution = solutions[i]
    data, bbsttdfp, ahdfpm, cgdfp, mopcg = solution
    original_signal, observed_signal = data[:original], data[:observed]
    bbsttdfp_reconstructed = bbsttdfp[:reconstructed]
    ahdfpm_reconstructed = ahdfpm[:reconstructed]
    cgdfp_reconstructed = cgdfp[:reconstructed]
    mopcg_reconstructed = mopcg[:reconstructed]
    (i,
        MSI(original_signal, bbsttdfp_reconstructed),
        MSI(original_signal, ahdfpm_reconstructed),
        MSI(original_signal, cgdfp_reconstructed),
        MSI(original_signal, mopcg_reconstructed)
    )
end
df = DataFrame(
    :experiment => map(x -> x[1], msis),
    :bbsttdfp => map(x -> x[2], msis),
    :ahdfpm => map(x -> x[3], msis),
    :cgdfp => map(x -> x[4], msis),
    :mopcg => map(x -> x[5], msis),
)
# msis, sol = msis_full[1:end][1:3], msis_full[1:end][4]
flder_to_save = "./results/stored_data"
writedlm("$flder_to_save/cs_proccessing_2023_05_22_11_00.csv", Iterators.flatten(([names(df)], eachrow(df))), ',')
sorted_msis = sort(msis, lt=(x, y) -> x[2] < y[2])
println(sorted_msis)
# p1, p2, p3, p4 = createCSPlots(original_signal, observed_signal, reconstructed)
# p4
# PSNR(original_signal, reconstructed)
solutions
