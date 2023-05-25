abstract type GenericAlogrithm end
abstract type Ding2017Algorithm <: GenericAlogrithm end
abstract type BBSTTDFPAlgorithm <: GenericAlogrithm end
abstract type BBSTTDFPAlgorithmNoIN <: GenericAlogrithm end
abstract type MOPCGAlgorithm <: GenericAlogrithm end
abstract type FISAPAlgorithm <: GenericAlogrithm end
abstract type CGDFPAlgorithm <: GenericAlogrithm end
abstract type AHDFPMAlgorithm <: GenericAlogrithm end


abstract type GenericProblem end
abstract type Ding2017Problem <: GenericProblem end
abstract type BBSTTDFPProblem <: GenericProblem end
abstract type MOPCGProblem <: GenericProblem end

struct Parameters
    hh::Float64
    r::Float64
    ψ::Float64
    θₘᵢₙ::Float64
    θₘₐₓ::Float64
    α::Float64
    σ::Float64
end

Base.show(io::IO, p::Parameters) = print(io, "hh=$(p.hh), r=$(p.r), ψ=$(p.ψ), θₘᵢₙ=$(p.θₘᵢₙ), θₘₐₓ=$(p.θₘₐₓ), α=$(p.α), σ=$(p.σ)")
t, ϵ, β, σ, γ = 0.5, 1e-5, 0.5, 0.5, 1.5
αmin, αmax = 0.5, 0.8
r, ψ = 0.5, 0.5

struct BBSTTDFPParams
    t::Float64
    β::Float64
    σ::Float64
    γ::Float64
    αmin::Float64
    αmax::Float64
    r::Float64
    ψ::Float64
    η::Union{Nothing,Float64}
    ξ::Union{Nothing,Float64}
end
BBSTTDFPParams(p::BBSTTDFPParams, σ::Float64) = BBSTTDFPParams(p.t, p.β, σ, p.γ, p.αmin, p.αmax, p.r, p.ψ, nothing, nothing)
BBSTTDFPParams(t::Float64, β::Float64, σ::Float64, γ::Float64, αmin::Float64, αmax::Float64, r::Float64, ψ::Float64) = BBSTTDFPParams(t, β, σ, γ, αmin, αmax, r, ψ, nothing, nothing)
Base.show(io::IO, p::BBSTTDFPParams) =
    if isnothing(p.η)
        print(io, "t=$(p.t), β=$(p.β),  σ=$(p.σ), γ=$(p.γ), αmin=$(p.αmin), αmax=$(p.αmax), r=$(p.r), ψ=$(p.ψ)")
    else
        print(io, "t=$(p.t), β=$(p.β),  σ=$(p.σ), γ=$(p.γ), αmin=$(p.αmin), αmax=$(p.αmax), r=$(p.r), ψ=$(p.ψ), η=$(p.η), ξ=$(p.ξ)")
    end

struct MOPCGParameters
    ρ::Float64
    η::Float64
    ξ::Float64
    λ::Float64
end

struct FISAPParameters
    σ::Float64
    r::Float64
    ρ::Float64
    ρk::Float64
    ν::Float64
    τ::Float64
end
struct CGDFPParameters
    ρ::Float64
    β::Float64
    σ::Float64
    σ1::Float64
    σ2::Float64
end

struct AHDFPMParameters
    ρ::Float64
    γ::Float64
    σ::Float64
    μ::Float64
    ν::Float64
    η1::Float64
    η2::Float64
    η3::Float64
    ξ::Float64
end



struct AlgorithmOptions
    maxiters::Int64
    tol::Float64
    maxevals::Union{Int64,Float64}
    stopping::Function
end
AlgorithmOptions(maxiters::Int64, tol::Float64) = AlgorithmOptions(maxiters, tol, Inf, (args...) -> false)
Base.show(io::IO, o::AlgorithmOptions) = print(io, "maxiters=$(o.maxiters), tol=$(o.tol), maxevals=$(o.maxevals)")

struct Ding2017Parameters
    ρ::Float64
    σ::Float64
    θ::Float64
    η::Float64
    ξ::Float64
end

Base.show(io::IO, p::Ding2017Parameters) = print(io, "ρ=$(p.ρ), σ=$(p.σ), θ=$(p.θ), η=$(p.η), ξ=$(p.ξ)")

struct TestProblem
    name::String
    f::Function
    P::Function
    Pcheck::Function
    x0::Vector{Float64}
    x0label::String
end
# TestProblem(name::String, f::Function, P::Function, x0::Vector{Float64}) = TestProblem(name, f, P, x0, (-Inf64, Inf64))
# TestProblem(name::String, f::Function, P::Function, x0::Vector{Float64}, l::Float64) = TestProblem(name, f, P, x0, (l, Inf64))
# TestProblem(name::String, f::Function, P::Function, x0::Vector{Float64}, bounds::Tuple{Float64,Float64}) = TestProblem(name, f, P, x0, bounds)

struct SuccessfullResult
    iterations::Int64
    functionEvalauations::Int64
    solution::Vector{Float64}
    Fnorm::Float64
    exec_time::Union{Nothing,Float64}
    flag::Union{Nothing,Symbol}
end
SuccessfullResult(iter::Int64, funevals::Int64, sols::Vector{Float64}, fnorm::Float64) =
    SuccessfullResult(iter, funevals, sols, fnorm, nothing, nothing)
SuccessfullResult(iter::Int64, funevals::Int64, sols::Vector{Float64}, fnorm::Float64, flag::Symbol) =
    SuccessfullResult(iter, funevals, sols, fnorm, nothing, flag)
SuccessfullResult(result::SuccessfullResult, exe_time::Float64) =
    SuccessfullResult(result.iterations, result.functionEvalauations, result.solution, result.Fnorm, exe_time, result.flag)

struct FailedResult
    message::String
end

struct NumericalResult
    problem::TestProblem
    dim::Int64
    options::AlgorithmOptions
    params::Union{Ding2017Parameters,Parameters,BBSTTDFPParams,MOPCGParameters,FISAPParameters,CGDFPParameters,AHDFPMParameters}
    result::Union{SuccessfullResult,FailedResult}
end
# NumericalResult(problem::TestProblem, dim::Int64, message::String) = NumericalResult(problem, dim, message, nothing, nothing, nothing, nothing)
