include("dependencies.jl")
include("included_files.jl")
include("algorithms_parameters.jl")


"""
Solve 
        min ||Ax-b||² + λ||x||₁  or   min ||Ax-b||² + τ||x||₁ 
"""

function getConstructedAndMSE(sol::NumericalResult, original::Vector{Float64}, n::Int64, name::String, experiment::Int64)
    z_bbs = sol.result.solution
    reconstructed = z_bbs[1:n] - z_bbs[n+1:end]
    mse_val = MSI(original, reconstructed)
    df = numericalResult2DataFrame(sol, algorithm=name)
    hcat(DataFrame(experiment_no=[experiment], mse=[mse_val]), df)
end

do_expriment(;
    name::Int64=0,
    k::Int64=7,
    σ::Float64=0.01,
    options::AlgorithmOptions=AlgorithmOptions(5000, 1e-12),
    params::BBSTTDFPParams=BBSTTDFPParams(),
    d::Distribution=Normal(0, 0.01)
) = begin

    n = 2^k
    m, r = fld(n, 4), fld(n, 8)
    # flder_to_save = "./results/stored_data"
    # ATA, initial_point, c, original_signal, observed_signal = getSavedData(flder_to_save, m, n, r)
    ATA, initial_point, c, original_signal, observed_signal = createCSData(m, n, r, σ; d=d)
    printstyled("Creating data for experiment $name \n"; color=:bold)
    problems = createCSProblems(ATA, vec(initial_point), vec(c), n)
    problem = problems[1]

    printstyled("Solving with BBSTTDFP for experiment $name \n"; color=:cyan)
    sol_bbs = BBSTTDFP(problem; options, params, inertial=false)

    printstyled("Solving with IBBSTTDFP for experiment $name \n"; color=:light_magenta)
    sol_ibbs = BBSTTDFP(problem; options, params, inertial=true)

    printstyled("Solving with AHDFPM for experiment $name \n"; color=:light_green)
    sol_ahd = AHDFPM(problem; options=options, params=AHDFPMParameters(0.6, 1.75, 0.0001, 0.001, 0.4, 2.4, 0.8, 0.1, 1))

    printstyled("Solving with CGDFP for experiment $name \n"; color=:magenta)
    sol_cgd = CGDFP(problem; options=options, params=CGDFPParameters(0.6, 1, 0.001, 0.7, 0.3))

    printstyled("Solving with MOPCG for experiment $name \n"; color=:light_red)
    sol_mopcg = MOPCG(problem; options=options, params=MOPCGParameters(1, 0.6, 0.1, 0.2))

    solution = (Dict(
        :name => "DATA",
        :dims => (n, m, r),
        :ATA => ATA,
        :c => c,
        :problem => problem,
        :initial => initial_point,
        :original => original_signal,
        :observed => observed_signal,
        :df => vcat(
            getConstructedAndMSE(sol_ibbs, original_signal, n, "IBBSTTDFP", name),
            getConstructedAndMSE(sol_bbs, original_signal, n, "BBSTTDFP", name),
            getConstructedAndMSE(sol_ahd, original_signal, n, "AHDFPM", name),
            getConstructedAndMSE(sol_cgd, original_signal, n, "CGDFP", name),
            getConstructedAndMSE(sol_mopcg, original_signal, n, "MOPCG", name)
        )
    )
    )
    solution
end
no_of_exprs = 100
algo_names = [
    "IBBSTTDFP",
    "BBSTTDFP",
    "AHDFPM",
    "CGDFP",
    "MOPCG",
]
params = getAlgorithmParams("cs_expr", 1e-7; newparams=true)
Random.seed!(123)
d = Normal(0, 0.01)
sols = map(1:no_of_exprs) do i
    printstyled("Starting experiment $i \n"; color=:green)
    do_expriment(name=i, k=11, σ=0.01, d=d,
        params=params,
        options=AlgorithmOptions(2_000, 1e-7, Inf, (fx0, fx1) -> (norm(fx0 - fx1) / norm(fx0)) <= 1e-7)
    )
end
colors = [:red, :black, :green, :blue, :purple]
solutions = vcat(map(sol -> sol[:df], sols)...)
data = readPaperData(solutions, algo_names)
T1 = dataFrames2Matrix(data, :mse)
T2 = dataFrames2Matrix(data, :time)

p1 = plotPerformanceProfile(T1, "MSE", algo_names; logscale=false)
p2 = plotPerformanceProfile(T2, "CPU TIME", algo_names; logscale=false)
folder = "./results/paper/data/cs"
map([(p1, "MSE"), (p2, "CPU TIME")]) do (p, title)
    pyplot()
    file_name = replace(title, " " => "_") |> lowercase
    printstyled("saving plots for $title to file \n"; color=:reverse)
    png_file = outputfilename(file_name, "png"; root=folder, suffix="cs_pp")
    savefig(p, png_file)
    svg_file = outputfilename(file_name, "svg"; root=folder, suffix="cs_pp")
    savefig(p, svg_file)
    pdf_file = outputfilename(file_name, "pdf"; root=folder, suffix="cs_pp")
    savefig(p, pdf_file)
    # eps_file = outputfilename(file_name, "eps"; root=folder, suffix="cs_pp")
    # savefig(p, eps_file)
end
file_name = outputfilename("cs", "csv"; root=folder, suffix="df")
numericalResultDF2CSV(solutions, file_name)
