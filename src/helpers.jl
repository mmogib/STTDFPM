function cleanVectorFromEmpty(x, v::Float64)
    x .|> d -> isempty(d) ? v : d
end
function createTestProlem(name::String, F::Function, P::Function, Pcheck::Function, x0::Vector{<:Real}, x0label::String)
    TestProblem(
        name,
        F,
        P,
        Pcheck,
        x0,
        x0label
    )
end

function createTestProlem(fns::Vector{Tuple{String,Function,Function,Function}}, x::Tuple{String,Vector{Float64}})
    map(fns) do fn
        name, F, P, Pcheck = fn
        x0label, x0 = x
        createTestProlem(name, F, P, Pcheck, x0, x0label)
    end
end

createParams(::Type{T}) where {T<:GenericProblem} = createParams(T, 1000)
function createParams(::Type{BBSTTDFPProblem}, sz::Int64, rng::MersenneTwister)
    ϵ = eps(1.0)
    ηs = rand(rng, ϵ:ϵ:0.5, sz, 1)
    ξs = ηs .+ rand(rng, ϵ:ϵ:(0.5-ϵ), sz, 1)
    hcat(ηs, ξs)
end
function createParams(::Type{BBSTTDFPProblem},
    sz::Int64;
    rng::MersenneTwister=MersenneTwister(parse(Int, Dates.format(now(), "ddHHMM"))),
    newparams::Bool=false
)
    rand_init = rand(rng, sz, 8)
    rand_init[:, 4] = rand_init[:, 4] + rand(rng, sz)
    rand_init[:, 6] = rand_init[:, 5] + rand(rng, sz)
    rand_init = newparams ? hcat(rand_init, createParams(BBSTTDFPProblem, sz, rng)) : rand_init
    params = map(p -> BBSTTDFPParams(p...), eachrow(rand_init))
    params
end

function numericalResult2DataFrame(result::NumericalResult; algorithm::String="")
    algorithm = algorithm
    problem_name = result.problem.name
    x0_label = result.problem.x0label
    dim = result.dim
    if (isa(result.result, SuccessfullResult))
        iters = result.result.iterations
        evals = result.result.functionEvalauations
        fNorm = @sprintf("%.5e", result.result.Fnorm)
        exe_time = result.result.exec_time
        message = isnothing(result.result.flag) ? "✔ problem solved" : "🤔 maximum iterations reached"
        return DataFrame(
            algorithm=[algorithm],
            problem_name=[problem_name],
            x0=[x0_label],
            dim=[dim],
            iters=[iters],
            evals=[evals],
            time=isnothing(exe_time) ? ["--"] : [exe_time],
            fnorm=[result.result.Fnorm],
            fNorm=[fNorm],
            message=[message],
            options=[result.options],
            params=[result.params]
        )

    else
        return DataFrame(
            algorithm=[algorithm],
            problem_name=[problem_name],
            x0=[x0_label],
            dim=[dim],
            iters=[""],
            evals=[""],
            time=[""],
            fnorm=[""],
            fNorm=[""],
            message=["❌😒" * result.result.message],
            options=[result.options],
            params=[result.params]
        )
    end
end


function numericalResult2DataFrame(results::Vector{NumericalResult}, algorithm::String="")
    df = DataFrame()
    dfs = results .|> result -> numericalResult2DataFrame(result, algorithm=algorithm)
    vcat(df, dfs...)
end

function numericalResultDF2CSV(df::DataFrame, filename::String)
    df |> CSV.write(filename)
end



function optimizeParameters(p::TestProblem)
    optimizeParameters([p])
end
function optimizeParameters(ps::Vector{TestProblem})
    function goodParams(vparams::Vector{Float64}, ps::Vector{TestProblem})
        params = BBSTTDFPParams(vparams...)
        options = AlgorithmOptions(50_000, 1e-5)

        results = ps .|> p -> BBSTTDFP(p; params=params, options=options)
        if all(map(result -> isa(result.result, SuccessfullResult), results))
            return max(map(result -> result.result.iterations, results)...)
        else
            return options.maxiters
        end
    end
    names_of_all = ps .|> (p -> p.name) |> names -> join(names, ",")
    println("Working on problems $(names_of_all), please wait..")
    ga = GA(populationSize=100,
        selection=susinv,
        crossover=DC,
        mutation=PLM()
    )
    lower = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0002, 0.0001, 0.0001]
    upper = [0.9999, 0.9999, 0.9999, 1.9999, 0.9999, 0.9999, 0.9999, 0.9999]
    result = Evolutionary.optimize(
        x -> goodParams(x, ps),
        BoxConstraints(lower, upper),
        ones(8),
        ga
    )
    root_dir = outputfolder()
    params = BBSTTDFPParams(Evolutionary.minimizer(result)...)
    options = AlgorithmOptions(50000, 1e-5)
    for p in ps
        println("Solving problem $p.name with optimized parameters...")
        result_p = BBSTTDFP(p; params=params, options=options)
        dfs = numericalResult2DataFrame(result_p)
        println("Saving to file ... ")
        numericalResultDF2CSV(dfs, "$root_dir/optimized_params_$(p.name).csv")
    end
end


function getNewParams(oldparams::BBSTTDFPParams, paramsToRenew::Vector{Symbol})
    t, β, σ, γ, αmin = oldparams.t, oldparams.β, oldparams.σ, oldparams.γ, oldparams.αmin
    αmax, r, ψ, η, ξ = oldparams.αmax, oldparams.r, oldparams.ψ, oldparams.η, oldparams.ξ
    # t::Float64 ∈ (0,1) ~ 0.11
    t = (:t in paramsToRenew) ? rand() : t
    # β::Float64 ∈ (0,1) ~ 0.5
    β = (:β in paramsToRenew) ? rand() : β
    # σ::Float64 ∈ (0,1) ~ 0.01
    σ = (:σ in paramsToRenew) ? rand() : σ
    # γ::Float64 ∈ (0,2) ~ 1.8 or 1.72
    γ = (:γ in paramsToRenew) ? rand([0, 1]) + rand() : γ
    # αmin::Float64 0 < αmin < αmax < 1 ~ 0
    αmin = (:αmin in paramsToRenew) ? rand() : αmin
    # αmax::Float64 0 < αmin < αmax < 1 ~ Inf64
    αmax = (:αmax in paramsToRenew) ? αmin + rand() : αmax
    # r::Float64  ∈ ( 0,1) ~ 0.1
    r = (:r in paramsToRenew) ? rand() : r
    # ψ::Float64 ∈ (0 ,1) ~ 0.5 
    ψ = (:ψ in paramsToRenew) ? rand() : ψ
    # η::Union{Nothing,Float64} optional 0 < η <= ξ  ~ 0.001 == μ in a2023LiuWuShaoZhangCao
    η = (:η in paramsToRenew) ? rand() : η
    # ξ::Union{Nothing,Float64} optional 0 < η <= ξ  ~ 0.6 == ν in a2023LiuWuShaoZhangCao
    ξ = (:ξ in paramsToRenew) ? rand() : ξ


    BBSTTDFPParams(t, β, σ, γ, αmin, αmax, r, ψ, η, ξ)
end

function str2params(strparams::String, rng::MersenneTwister;
    randomized::Union{Nothing,Vector{Symbol}}=nothing,
    newparams::Bool=false
)

    parts = Dict(split(strparams, ",") .|> x -> split(x, "=") .|> x -> strip.(x))

    namedtuples = (; zip(Symbol.(keys(parts)), parse.(Float64, (values(parts))))...)


    getnewparam(s::Symbol, passedparams::NamedTuple, randomized::Union{Nothing,Vector{Symbol}}) = begin
        if isnothing(randomized)
            get(passedparams, s, rand(rng))
        else
            s in randomized ? rand(rng) : get(passedparams, s, rand(rng))
        end
    end
    pms = BBSTTDFPParams(
        getnewparam(:t, namedtuples, randomized),
        getnewparam(:β, namedtuples, randomized),
        getnewparam(:σ, namedtuples, randomized),
        getnewparam(:γ, namedtuples, randomized),
        getnewparam(:αmin, namedtuples, randomized),
        getnewparam(:αmax, namedtuples, randomized),
        getnewparam(:r, namedtuples, randomized),
        getnewparam(:ψ, namedtuples, randomized),
    )
    npms = if newparams
        ϵ = eps(1.0)
        η = rand(rng, ϵ:ϵ:0.2)
        ξ = rand(rng, (η+0.3):ϵ:(1-ϵ))
        BBSTTDFPParams(pms.t,
            pms.β,
            pms.σ,
            pms.γ,
            pms.αmin,
            pms.αmax,
            pms.r,
            pms.ψ,
            η,
            ξ
        )
    else
        pms
    end
    npms
end

function outputfolder(root::String="./results")
    t = now()
    dayfolder = Dates.format(t, "yyyy_mm_dd")
    hourfolder = Dates.format(t, "HH_MM")
    root_dir = mkpath("$root/$dayfolder/$hourfolder")
    return root_dir
end

function outputfilename(name::String, extension::String;
    root::String="./results",
    suffix::Union{Nothing,String}=nothing)
    root_dir = outputfolder(root)
    filename = if isnothing(suffix)
        "$root_dir/$name.$extension"
    else
        "$root_dir/$name_$suffix.$extension"
    end
    filename
end

function runNumericalExperiments(::Type{Ding2017Problem}, ps::Vector{TestProblem}; save::Bool=false)
    println("Running Ding 2017 Experiments for comparing.")
    options = AlgorithmOptions(500_000, 1e-5)
    params = Ding2017Parameters(0.74, 1e-4, 0.1, 0.5, 1)

    dfs = ps .|> problem -> begin
        println("Solving $(problem.name)...")
        ding2017(problem; options=options, params=params)
    end |> numericalResult2DataFrame


    df = vcat(dfs...)
    if save
        filename = outputfilename("ding2017", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end

function runNumericalExperiments(::Type{BBSTTDFPProblem}, ps::Vector{TestProblem},
    params::BBSTTDFPParams,
    options::AlgorithmOptions;
    save::Bool=false,
    experiment_name::Union{Nothing,String}=nothing
)

    println("Running Experiments $(isnothing(experiment_name) ? "" : "($experiment_name)") testing new algorithm.")
    dfs = ps .|> problem -> begin
        println("Solving ... $(problem.name)")
        results = BBSTTDFP(problem; params, options=options)
        numericalResult2DataFrame(results, algorithm="BBSTTDFPAlgorithm")
    end
    df = vcat(dfs...)
    if (save)
        println("Saving to file ... ")
        filename = isnothing(experiment_name) ? "" : replace(experiment_name, " " => "_")
        file_name = outputfilename(filename, "csv")
        numericalResultDF2CSV(df, file_name)
        println("DONE!")
    end
    df
end
function runNumericalExperiments(
    ::Type{BBSTTDFPParams},
    problems::Vector{TestProblem};
    save::Bool=false,
    experiments_count::Int64=5,
    paramstr::Vector{String}=["t=0.13824781669963615, β=0.5909641129216927,  σ=0.25514946641464453, γ=1.4117390335497881, αmin=0.912323896610999, αmax=1.693242503065831, r=0.8025450374268133, ψ=0.6642550140797701"],
    randomized::Union{Nothing,Vector{Symbol}}=nothing,
    options::AlgorithmOptions=AlgorithmOptions(100_000, 1e-5),
    rng::MersenneTwister=MersenneTwister(parse(Int, Dates.format(now(), "ddHHMM"))),
    newparams::Bool=false
)
    params = Iterators.flatten(str2params.(paramstr, rng; randomized=randomized, newparams=newparams) for _ in 1:(experiments_count)) |> collect


    dfs = enumerate(params) .|> p -> begin
        i, param = p
        runNumericalExperiments(BBSTTDFPProblem, problems, param, options; experiment_name="experiment $i")
    end
    df = vcat(dfs...)
    if save
        file_name = outputfilename("testing_bbsttdp", "csv")
        numericalResultDF2CSV(df, file_name)
    end
    dfs
end

function runNumericalExperiments(
    ::Type{BBSTTDFPParams},
    ps::Vector{TestProblem},
    itrs::Int64,
    szs::Int64;
    rng::MersenneTwister=MersenneTwister(parse(Int, Dates.format(now(), "ddHHMM"))),
    newparams::Bool=false
)
    println("Running Experiments to pick best parameters.")

    params = createParams(BBSTTDFPProblem, szs; rng=rng, newparams=newparams)
    options = AlgorithmOptions(itrs, 1e-5)
    dfs = enumerate(ps) .|> p -> begin
        i, problem = p
        println("Solving ... $(problem.name)")
        results = BBSTTDFP(problem, params; options=options)

        numericalResult2DataFrame(results)
    end
    df = vcat(dfs...)

    println("Saving to file ... ")
    file_name = outputfilename("paramters", "csv")
    numericalResultDF2CSV(df, file_name)
    println("DONE!")
    return df
end



function runNumericalExperiments(::Type{MOPCGAlgorithm},
    ps::Vector{TestProblem};
    save::Bool=false,
    options::AlgorithmOptions=AlgorithmOptions(),
    params::MOPCGParameters=MOPCGParameters()
)
    println("Running MOPCG - SabiuShahStanimirovicIvanovWaziri (2023) Experiments for comparing.")

    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("MOPCG($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        MOPCG(problem; options=options, params=params)
    end |> res -> numericalResult2DataFrame(res, algorithm="MOPCGAlgorithm")


    df = vcat(dfs...)
    if save
        filename = outputfilename("MOPCG", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end


function runExperiment(::Type{MOPCGAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false)

    problems = createTestProlems(dim)
    params = MOPCGParameters(0.1, 0.9, 0.0001, 0.1)
    dfs = runNumericalExperiments(
        MOPCGAlgorithm,
        problems;
        save=save,
        options=options,
        params=params
    )
    dfs

end

function runNumericalExperiments(::Type{FISAPAlgorithm},
    ps::Vector{TestProblem};
    save::Bool=false,
    options::AlgorithmOptions=AlgorithmOptions(),
    params::FISAPParameters=FISAPParameters()
)
    println("Running FISAP Algorithm - ZhangLiuZhangLu (2023) Experiments for comparing.")

    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("ISAP($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        FISAP(problem; options=options, params=params)
    end |> res -> numericalResult2DataFrame(res, algorithm="FISAPAlgorithm ")


    df = vcat(dfs...)
    if save
        filename = outputfilename("FISAP", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end
function runExperiment(::Type{FISAPAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false)

    problems = createTestProlems(dim)
    params = FISAPParameters(0.01, 0.5, 0.24, 1.7, 0.001, 1 / 12)
    dfs = runNumericalExperiments(
        FISAPAlgorithm,
        problems;
        save=save,
        options=options,
        params=params
    )
    dfs

end

function runNumericalExperiments(::Type{CGDFPAlgorithm},
    ps::Vector{TestProblem};
    save::Bool=false,
    options::AlgorithmOptions=AlgorithmOptions(),
    params::CGDFPParameters=CGDFPParameters()
)
    println("Running CGDFP Algorithm - 2020ZhengYangLiang (2020) Experiments for comparing.")

    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("CGDFP($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        CGDFP(problem; options=options, params=params)
    end |> res -> numericalResult2DataFrame(res, algorithm="CGDFPAlgorithm ")


    df = vcat(dfs...)
    if save
        filename = outputfilename("CGDFP", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end
function runExperiment(::Type{CGDFPAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false)

    problems = createTestProlems(dim)
    params = CGDFPParameters(0.6, 1, 0.001, 0.7, 0.3)
    dfs = runNumericalExperiments(CGDFPAlgorithm, problems;
        save=save,
        options=options,
        params=params
    )
    dfs

end



function runNumericalExperiments(::Type{AHDFPMAlgorithm},
    ps::Vector{TestProblem};
    save::Bool=false,
    options::AlgorithmOptions=AlgorithmOptions(),
    params::AHDFPMParameters=AHDFPMParameters()
)
    println("Running AHDFPM Algorithm Algorithm - 2023LiuWuShaoZhangCao (2023) Experiments for comparing.")
    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("AHDFPM($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        AHDFPM(problem; options=options, params=params)
    end |> res -> numericalResult2DataFrame(res, algorithm="AHDFPMAlgorithm ")


    df = vcat(dfs...)
    if save
        filename = outputfilename("AHDFPM", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end

function runExperiment(::Type{AHDFPMAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false)

    problems = createTestProlems(dim)
    params = AHDFPMParameters(0.6, 1.72, 0.001, 0.001, 0.6, 2.3, 0.7, 0.01, 1)
    dfs = runNumericalExperiments(AHDFPMAlgorithm, problems;
        save=save,
        options=options,
        params=params
    )
    dfs

end

# 

# function runExperiment(::Type{BBSTTDFPAlgorithm}, dim::Int64, options::AlgorithmOptions)

#     problems = createTestProlems(dim)
#     paramstr = [
#         # "t=0.11408250586841473, β=0.6570240538036174,  σ=0.008575295412900363, γ=1.8928223633049506, αmin=0.10912531544765902, αmax=0.8312657796763663, r=0.12207595744012179, ψ=0.20929112731122723",
#         # "t=0.08117923071701227, β=0.9863891970568359,  σ=0.012299072746424944, γ=1.0745780136703662, αmin=0.7792661839876198, αmax=0.898216880691173, r=0.4521518638573503, ψ=0.9914360613166049",
#         # "t=0.11408250586841473, β=0.6570240538036174,  σ=0.008575295412900363, γ=1.8928223633049506, αmin=0.10912531544765902, αmax=0.8312657796763663, r=0.1452713517728983, ψ=0.20929112731122723",
#         # "t=0.9104508426215066, β=0.6432565761660647,  σ=0.09303982544261569, γ=1.7165457478753454, αmin=0.6963744294272964, αmax=1.5275295393983448, r=0.7460784575151529, ψ=0.18304771537151976",
#         # "t=0.6745380844386832, β=0.4188321351360458,  σ=0.09303982544261569, γ=1.7165457478753454, αmin=0.6963744294272964, αmax=1.5275295393983448, r=0.7460784575151529, ψ=0.18304771537151976",
#         "t=0.11408250586841473, β=0.6570240538036174,  σ=0.008575295412900363, γ=1.8928223633049506, αmin=0.10912531544765902, αmax=0.8312657796763663, r=0.12207595744012179, ψ=0.20929112731122723"
#     ]
#     rng = MersenneTwister(230505151000)
#     dfs = runNumericalExperiments(BBSTTDFPAlgorithm,
#         problems;
#         save=false, experiments_count=1,
#         # randomized=[:t],
#         paramstr=paramstr,
#         options=AlgorithmOptions(2000, 1e-11),
#         rng=rng
#     )
#     dfs

# end
function runExperiment(::Type{BBSTTDFPAlgorithmNoIN}, dim::Int64, options::AlgorithmOptions; save::Bool=false, newparams::Bool=false)
    problems = createTestProlems(dim)
    # t=0.11408250586841473, β=0.6570240538036174,  
    # σ=0.008575295412900363, γ=1.8928223633049506, 
    # αmin=0.10912531544765902, αmax=0.8312657796763663, 
    # r=0.12207595744012179, ψ=0.20929112731122723
    params0 = BBSTTDFPParams(
        0.11,
        0.5,
        0.01,
        1.8,
        1e-30,
        1e+30,
        0.1,
        0.5,
        newparams ? 0.001 : nothing,
        newparams ? 0.6 : nothing
    )
    params1 = BBSTTDFPParams(
        0.11408250586841473,
        0.6570240538036174,
        0.008575295412900363,
        1.8928223633049506,
        0.10912531544765902,
        0.8312657796763663,
        0.12207595744012179,
        0.20929112731122723,
        newparams ? 0.001 : nothing,
        newparams ? 0.6 : nothing
    )
    params2 = BBSTTDFPParams(
        0.19091055694541392,
        0.8248500773070584,
        0.6059149417443106,
        1.0160210808836148,
        0.968177563496758,
        1.8227259069163115,
        0.9507650425301077,
        0.8808475746290094,
        newparams ? 0.030826250174403214 : nothing,
        newparams ? 0.3719729951071617 : nothing
    )

    dfs = runNumericalExperiments(
        BBSTTDFPAlgorithmNoIN,
        problems;
        save=save,
        options=options,
        params=params0
    )
    dfs
end

function runNumericalExperiments(::Type{BBSTTDFPAlgorithmNoIN},
    ps::Vector{TestProblem};
    save::Bool=false,
    options::AlgorithmOptions=AlgorithmOptions(1000, 1e-6),
    params::BBSTTDFPParams=BBSTTDFPParams(
        0.11408250586841473,
        0.6570240538036174,
        0.008575295412900363,
        1.8928223633049506,
        0.10912531544765902,
        0.8312657796763663,
        0.12207595744012179,
        0.20929112731122723,
        0.001,
        0.6
    )
)
    # t=0.11408250586841473, β=0.6570240538036174,  σ=0.008575295412900363, γ=1.8928223633049506, αmin=0.10912531544765902, αmax=0.8312657796763663, r=0.12207595744012179, ψ=0.20929112731122723
    println("Running our BBSTTDFP_NO_INERTIAL Algorithm")
    counter = 1
    total = length(ps)
    title_appendage = isnothing(params.η) ? "old_params" : "new_params"
    dfs = ps .|> problem -> begin
        println("BBSTTDFP_NO_INERTIAL($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        BBSTTDFP(problem, params=params, options=options; inertial=false)

    end |> res -> numericalResult2DataFrame(res, algorithm="BBSTTDFPNoINAlgorithm($title_appendage)")


    df = vcat(dfs...)
    if save
        filename = outputfilename("BBSTTDFP_NO_INERTIAL_$title_appendage", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end

function runNumericalExperiments(::Type{BBSTTDFPAlgorithm},
    ps::Vector{TestProblem};
    save::Bool=false,
    options::AlgorithmOptions=AlgorithmOptions(1000, 1e-6),
    params::BBSTTDFPParams=BBSTTDFPParams(
        0.11408250586841473,
        0.6570240538036174,
        0.008575295412900363,
        1.8928223633049506,
        0.10912531544765902,
        0.8312657796763663,
        0.12207595744012179,
        0.20929112731122723,
        0.001,
        0.6
    )
)
    # t=0.11408250586841473, β=0.6570240538036174,  σ=0.008575295412900363, γ=1.8928223633049506, αmin=0.10912531544765902, αmax=0.8312657796763663, r=0.12207595744012179, ψ=0.20929112731122723
    println("Running our BBSTTDF Algorithm")
    counter = 1
    total = length(ps)
    title_appendage = isnothing(params.η) ? "old_params" : "new_params"
    dfs = ps .|> problem -> begin
        println("BBSTTDFP($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        BBSTTDFP(problem, params=params, options=options)

    end |> res -> numericalResult2DataFrame(res, algorithm="BBSTTDFPAlgorithm($title_appendage)")


    df = vcat(dfs...)
    if save
        filename = outputfilename("BBSTTDFP_$title_appendage", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end

function runExperiment(::Type{BBSTTDFPAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false, newparams::Bool=false)
    problems = createTestProlems(dim)
    # t=0.11408250586841473, β=0.6570240538036174,  
    # σ=0.008575295412900363, γ=1.8928223633049506, 
    # αmin=0.10912531544765902, αmax=0.8312657796763663, 
    # r=0.12207595744012179, ψ=0.20929112731122723
    params0 = BBSTTDFPParams(
        0.11, # 0.98
        0.5,
        0.01,
        1.8,
        1e-30,
        1e+30,
        0.1,
        0.5,
        newparams ? 0.001 : nothing,
        newparams ? 0.6 : nothing
    )
    params1 = BBSTTDFPParams(
        0.11408250586841473,
        0.6570240538036174,
        0.008575295412900363,
        1.8928223633049506,
        0.10912531544765902,
        0.8312657796763663,
        0.12207595744012179,
        0.20929112731122723,
        newparams ? 0.001 : nothing,
        newparams ? 0.6 : nothing
    )
    params2 = BBSTTDFPParams(
        0.19091055694541392,
        0.8248500773070584,
        0.6059149417443106,
        1.0160210808836148,
        0.968177563496758,
        1.8227259069163115,
        0.9507650425301077,
        0.8808475746290094,
        newparams ? 0.030826250174403214 : nothing,
        newparams ? 0.3719729951071617 : nothing
    )

    dfs = runNumericalExperiments(
        BBSTTDFPAlgorithm,
        problems;
        save=save,
        options=options,
        params=params0
    )
    dfs
end


function runExperiment(algorithms::Vector{Tuple{T,Bool}} where {T<:DataType}, dims::Vector{Int64}, options::AlgorithmOptions;
    saveit::Bool=false)
    ## [*] Run our experiments  
    alldfs = map(dims) do dim
        df = map(algorithms) do (algo, flag)
            if flag # new params 
                runExperiment(algo, dim, options, save=saveit, newparams=true)
            else
                runExperiment(algo, dim, options, save=saveit)
            end
        end
        df
    end
    # dfs = 

    vcat(alldfs...)

end
function plotResultsProfiles(algorithms::Vector{Tuple{T,Bool}} where {T<:DataType}, dfs::Vector{DataFrame}; k=11)

    labels = algorithms .|> t -> replace("$(t[1])", "Algorithm" => "")
    labels = collect(1:length(algorithms)) .|> i -> algorithms[i][2] ? "$(labels[i])_NEW_PARAMS" : labels[i]
    # xlabel = join(labels, " vs ")
    xlabel = L"\tau"
    legendfontsize, leg = 8, :bottomright
    plts = map([
        (:iters, "Iterations"),
        (:evals, "Evalauations"),
        (:time, "CPU Time")]) do (item, title)
        T = map(dfs) do df
            values = df[!, item]
            val = Float64(maximum(values[findall(!isempty, values)]))
            values |> x -> cleanVectorFromEmpty(x, 2 * val) .|> Float64
        end
        theme(:default)
        performance_profile(
            PlotsBackend(),
            hcat(T...),
            labels,
            title=title,
            xlabel=xlabel,
            ylabel=L"\rho(\tau)",
            legendfontsize=legendfontsize,
            leg=leg,
            palette=:Dark2_5
        )
    end
    plts

end
function createTestProlems(dim::Int64)
    # PExponetialIII(x0=x0),
    # problem1 = createTestProlems(
    #     "ExponetialIII",
    #     x -> ExponetialIII(x),
    #     x -> projectOnBox(x; bounds=(0.0, Inf64)),
    #     x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
    #     dim
    # )
    problem2 = createTestProlems(
        "PolynomialI",
        x -> PolynomialI(x),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        dim
    )
    problem3 = createTestProlems(
        "SmoothSine",
        x -> SmoothSine.(x),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        dim
    )

    # problem4 = createTestProlems(
    #     "PolynomialSineCosine",
    #     x -> PolynomialSineCosine.(x),
    #     x -> projectOnBox(x; bounds=(0.0, Inf64)),
    #     x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
    #     dim
    # )
    problem5 = createTestProlems(
        "ExponetialI",
        x -> ExponetialI.(x),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        dim
    )
    problem6 = createTestProlems(
        "NonsmoothSine",
        x -> NonsmoothSine.(x),
        x -> projectOnTriangle(x; lower=0),
        x -> projectOnTriangleCheck(x; lower=0),
        dim
    )
    problem7 = createTestProlems(
        "ModifiedNonsmoothSine",
        x -> ModifiedNonsmoothSine.(x),
        x -> projectOnTriangle(x; lower=-1),
        x -> projectOnTriangleCheck(x; lower=-1),
        dim
    )
    problem8 = createTestProlems(
        "ModifiedNonsmoothSine2",
        x -> ModifiedNonsmoothSine2.(x),
        x -> projectOnTriangle(x; lower=-1),
        x -> projectOnTriangleCheck(x; lower=-1),
        dim
    )
    problem9 = createTestProlems(
        "ExponetialSineCosine",
        x -> ExponetialSineCosine.(x),
        x -> projectOnTriangle(x; lower=-1),
        x -> projectOnTriangleCheck(x; lower=-1),
        dim
    )
    problem10 = createTestProlems(
        "ModifiedTrigI",
        x -> ModifiedTrigI(x),
        x -> projectOnBox(x; bounds=(-3.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(-3.0, Inf64)),
        dim
    )
    # problem11 = createTestProlems(
    #     "ModifiedTrigII",
    #     x -> ModifiedTrigII(x),
    #     x -> projectOnBox(x; bounds=(-2.0, Inf64)),
    #     dim
    # )
    problem12 = createTestProlems(
        "ModifiedTridiagonal",
        x -> ModifiedTridiagonal(x),
        x -> projectOnTriangle(x; lower=0, β=1),
        x -> projectOnTriangleCheck(x; lower=0, β=1),
        dim
    )
    problem13 = createTestProlems(
        "Logarithmic",
        x -> Logarithmic(x),
        x -> projectOnBox(x; bounds=(-1.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(-1.0, Inf64)),
        dim
    )
    problem14 = createTestProlems(
        "NonmoothLogarithmic",
        x -> NonmoothLogarithmic(x),
        x -> projectOnBox(x; bounds=(0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0, Inf64)),
        dim
    )

    # problem15 = createTestProlems(
    #     "ARWHEADGrad",
    #     x -> ARWHEADGrad(x),
    #     x -> projectOnTriangle(x; lower=0),
    #     x -> projectOnTriangleCheck(x; lower=0),
    #     dim
    # )

    problem16 = createTestProlems(
        "ENGVAL1Grad",
        x -> ENGVAL1Grad(x),
        x -> x,
        x -> true,
        dim
    )

    return vcat(
        problem2,
        problem3,
        # problem4,
        problem5,
        problem6,
        problem7,
        problem8,
        problem9,
        problem10,
        problem12,
        problem13,
        problem14,
        problem16
    )

end
function createTestProlems(name::String, F::Function, P::Function, Pcheck::Function, dim::Int64)
    onez = ones(dim)
    indxs = 1:dim
    startin_points_1 = [
        ("x1: zeros", zeros(dim)),
        ("x2: 0.2", 0.2 * onez),
        ("x3: 0.4", 0.4 * onez),
        ("x4: 0.5", 0.5 * onez),
        ("x5: 0.6", 0.6 * onez),
        ("x6: 0.8", 0.8 * onez),
        # ("x6: -1", -1 * onez),
        ("x7: 1", onez),
        ("x8: 1.1", 1.1 * onez),
    ]
    startin_points_2 = [
        # ("x9: 1/2ᵏ", indxs .|> t -> 1 / 2^t),
        ("x10: 1 - 1/n", indxs .|> t -> 1 - (1 / dim)),
        ("x11: 1/k", indxs .|> t -> 1 / t),
        ("x12: (k-1)/n", indxs .|> t -> (t - 1) / dim),
        ("x13: 1/n", indxs .|> t -> 1 / dim),
        ("x14: 1/3ᵏ", indxs .|> t -> 1 / BigFloat(3^t)),
        ("x15: k/n", indxs .|> t -> t / dim)
    ]
    startin_points = [
        startin_points_1..., startin_points_2...
    ]
    # startin_points = [
    #     ("x1: 1", onez),
    #     ("x2: 0.2", 0.2 * onez),
    #     ("x3: 1/2ᵏ", indxs .|> t -> 1 / 2^t),
    #     ("x4: (k-1)/n", indxs .|> t -> (t - 1) / dim),
    #     ("x5: 1/k", indxs .|> t -> 1 / t),
    #     ("x6: 1/n", indxs .|> t -> 1 / dim),
    #     ("x7: 1 - 1/n", indxs .|> t -> 1 - (1 / dim)),
    #     ("x8: 1.1", 1.1 * onez),
    # ]
    map(startin_points) do (label, x0)
        createTestProlem(name, F, P, Pcheck, x0, label)
    end
end


"""
CS:Compressed Sensing helpers
"""

function createCSData(m, n, r, σ; d::Distribution=Normal())
    x = zeros(n)
    q = randperm(n)
    x[q[1:r]] = rand(d, r)
    A = rand(d, m, n)
    A = Matrix(qr(A').Q)'
    noise = σ * rand(d, m)
    b = A * x + noise
    x0 = A' * b
    τ = 0.01 * norm(x0, Inf)
    c = τ * ones(2n) + vcat(-x0, x0)
    z0 = vcat(max.(x0, 0), max.(-x0, 0))
    A' * A, z0, c, x, b
end

function createCSPlots(original, observed, reconstructed)
    p1 = bar(original, ylims=(-1.1, 1.1), bar_width=0.001, title="Original Signal", label=nothing)
    p2 = plot(observed, title="Observed Signal", label=nothing)
    p3 = bar(reconstructed, ylims=(-1.1, 1.1), bar_width=0.001, title="Reconstrcucted Signal", label=nothing)
    p4 = plot(p1, p2, p3, layout=(3, 1), size=(900, 300))
    p1, p2, p3, p4
end

# l₁ Regularized Least Squares 
function L1LS(z, ATA, c, n)
    u = z[1:n]
    v = z[n+1:2n]
    # z = vcat(u, v)
    Bu = ATA * (u - v)
    Hz = vcat(Bu, -Bu)
    min.(z, Hz + c)
end
function createCSProblems(ATA, z0, c, n)

    problem1 = TestProblem(
        "SignalProcessing - no projection",
        x -> L1LS(x, ATA, c, n),
        x -> x,
        x -> true,
        z0,
        "A'b"
    )
    problem2 = TestProblem(
        "SignalProcessing - with projection",
        x -> L1LS(x, ATA, c, n),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        z0,
        "A'b"
    )
    [
        problem1,
        problem2
    ]
end
function saveDataToFile(flder_to_save::String,
    ATA, initial_point, c, original_signal, observed_signal)

    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_ATA.txt", "w") do io
        writedlm(io, ATA)
    end
    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_initial_point.txt", "w") do io
        writedlm(io, initial_point)
    end
    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_c.txt", "w") do io
        writedlm(io, c)
    end
    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_original_signal.txt", "w") do io
        writedlm(io, original_signal)
    end
    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_observed_signal.txt", "w") do io
        writedlm(io, observed_signal)
    end

end

function getSavedData(flder_where_saved::String, m, n, r)
    ATA = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_ATA.txt")
    initial_point = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_initial_point.txt")
    c = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_c.txt")
    original_signal = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_original_signal.txt")
    observed_signal = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_observed_signal.txt")
    ATA, initial_point, c, original_signal, observed_signal
end

MSI(original, reconstructed) = (1 / length(original)) * norm(reconstructed - original)^2


PSNR(original, reconstructed) = 10 * log10(norm(reconstructed, Inf)^2 / MSI(original, reconstructed))