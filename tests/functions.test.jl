using LinearAlgebra, Random, BenchmarkTools, Dates, Printf, Test
using CUTEst, NLPModels
using Flux
include.(["../src/types.jl",
    "../src/functions.jl", "../src/test_proplems.jl"])


fARWHEADgrad(x) = gradient(ARWHEADFun, x)[1]
fPENALTY1grad(x) = gradient(PENALTY1Fun, x)[1]
fDIXON3DQgrad(x) = gradient(DIXON3DQFun, x)[1]
fGENHUMPSgrad(x) = gradient(GENHUMPSFun, x)[1]
fENGVAL1grad(x) = gradient(ENGVAL1Fun, x)[1]
fDIXMAANHgrad(x) = gradient(DIXMAANHFun, x)[1]
fDIXMAANIgrad(x) = gradient(DIXMAANIFun, x)[1]

rng = MersenneTwister(parse(Int, Dates.format(now(), "ddHHmm")))

@testset verbose = true "Testing `ARWHEAD` Function" begin
    @testset verbose = true "Compared with CUTEst" begin
        nlp = CUTEstModel("ARWHEAD")
        x0 = nlp.meta.x0
        @testset "Testing the function" begin
            @test obj(nlp, x0) ≈ ARWHEADFun(x0)
        end
        @testset "Testing the gradient" begin
            @test norm(grad(nlp, x0) - ARWHEADGrad(x0)) ≈ 0
        end
        finalize(nlp)
    end

    @testset "Compared with FlUX AD" begin

        x0 = ones(500)
        @test norm(ARWHEADGrad(x0) - fARWHEADgrad(x0)) ≈ 0
        for i in 1:10
            x0 = rand(rng, 500)
            @test norm(ARWHEADGrad(x0) - fARWHEADgrad(x0)) <= 1e-11
        end
    end

    @testset "CUTEst vs FlUX AD" begin
        nlp = CUTEstModel("ARWHEAD")
        x0 = nlp.meta.x0
        @test norm(grad(nlp, x0) - fARWHEADgrad(x0)) <= 1e-10
        finalize(nlp)
    end

end


@testset verbose = true "Testing `PENALTY1` Function" begin
    # test 1 with CUTEst
    @testset verbose = true "Compared with CUTEst" begin
        nlp = CUTEstModel("PENALTY1")
        x0 = nlp.meta.x0
        @testset "Testing the function" begin
            @test obj(nlp, x0) ≈ PENALTY1Fun(x0)
        end
        @testset "Testing the gradient" begin
            @test norm(grad(nlp, x0) - PENALTY1Grad(x0)) ≈ 0
        end
        finalize(nlp)
    end

    @testset "Compared with FlUX AD" begin
        x0 = ones(500)
        @test norm(PENALTY1Grad(x0) - fPENALTY1grad(x0)) ≈ 0
        for i in 1:10
            x0 = rand(rng, 500)
            @test norm(PENALTY1Grad(x0) - fPENALTY1grad(x0)) <= 1e-10
        end
    end

    @testset "CUTEst vs FlUX AD" begin
        nlp = CUTEstModel("PENALTY1")
        x0 = nlp.meta.x0
        @test norm(grad(nlp, x0) - fPENALTY1grad(x0)) <= 1e-10
        finalize(nlp)
    end
end

@testset verbose = true "Testing `DIXON3DQ` Function" begin
    # test 1 with CUTEst
    @testset verbose = true "Compared with CUTEst" begin
        nlp = CUTEstModel("DIXON3DQ")
        x0 = nlp.meta.x0
        @testset "Testing the function" begin
            @test obj(nlp, x0) ≈ DIXON3DQFun(x0)
        end

        @testset "Testing the gradient" begin
            @test norm(grad(nlp, x0) - DIXON3DQGrad(x0)) ≈ 0
        end
        finalize(nlp)
    end

    @testset "Compared with FlUX AD" begin
        x0 = -ones(500)
        @test norm(DIXON3DQGrad(x0) - fDIXON3DQgrad(x0)) ≈ 0
        for i in 1:10
            x0 = rand(rng, 500)
            @test norm(DIXON3DQGrad(x0) - fDIXON3DQgrad(x0)) <= 1e-11
        end
    end
    @testset "CUTEst vs FlUX AD" begin
        nlp = CUTEstModel("DIXON3DQ")
        x0 = nlp.meta.x0
        @test norm(grad(nlp, x0) - fDIXON3DQgrad(x0)) <= 1e-10
        finalize(nlp)
    end

end


@testset verbose = true "Testing `GENHUMPS` Function" begin
    # test 1 with CUTEst
    @testset verbose = true "Compared with CUTEst" begin
        nlp = CUTEstModel("GENHUMPS")
        x0 = nlp.meta.x0

        @testset "Testing the function" begin
            @test obj(nlp, x0) ≈ GENHUMPSFun(x0)
        end
        @testset "Testing the gradient" begin
            @test norm(grad(nlp, x0) - GENHUMPSGrad(x0)) <= 1e-11
        end
        finalize(nlp)
    end

    @testset "Compared with FlUX AD" begin
        x0 = [-506.0, -506.2 * ones(99)...]
        @test norm(GENHUMPSGrad(x0) - fGENHUMPSgrad(x0)) ≈ 0
        for i in 1:10
            x0 = rand(rng, 500)
            @test norm(GENHUMPSGrad(x0) - fGENHUMPSgrad(x0)) <= 1e-11
        end
    end
    @testset verbose = true "CUTEst vs FlUX AD" begin
        nlp = CUTEstModel("GENHUMPS")
        @testset "CUTEst vs FlUX AD with meta x0" begin
            x0 = nlp.meta.x0
            @test norm(grad(nlp, x0) - fGENHUMPSgrad(x0)) <= 1e-11
        end
        @testset "CUTEst vs FlUX AD with starting point" begin
            x0 = -[506.0, 506.2 * ones(length(nlp.meta.x0) - 1)...]
            @test norm(grad(nlp, x0) - fGENHUMPSgrad(x0)) <= 1e-11
        end
        finalize(nlp)
    end

end




@testset verbose = true "Testing `ENGVAL1` Function" begin
    # test 1 with CUTEst
    @testset verbose = true "Compared with CUTEst" begin
        nlp = CUTEstModel("ENGVAL1")
        x0 = nlp.meta.x0

        @testset "Testing the function" begin
            @test obj(nlp, x0) ≈ ENGVAL1Fun(x0)
        end
        @testset "Testing the gradient" begin
            @test norm(grad(nlp, x0) - ENGVAL1Grad(x0)) <= 1e-11
        end
        finalize(nlp)
    end

    @testset "Compared with FlUX AD" begin
        x0 = 2.0 * ones(100)
        @test norm(ENGVAL1Grad(x0) - fENGVAL1grad(x0)) ≈ 0
        for i in 1:10
            x0 = rand(rng, 500)
            @test norm(ENGVAL1Grad(x0) - fENGVAL1grad(x0)) <= 1e-11
        end
    end
    @testset verbose = true "CUTEst vs FlUX AD" begin
        nlp = CUTEstModel("ENGVAL1")
        @testset "CUTEst vs FlUX AD with meta x0" begin
            x0 = nlp.meta.x0
            @test norm(grad(nlp, x0) - fENGVAL1grad(x0)) <= 1e-11
        end
        finalize(nlp)
    end
end


@testset verbose = true "Testing `DIXMAANH` Function" begin
    # test 1 with CUTEst
    @testset verbose = true "Compared with CUTEst" begin
        nlp = CUTEstModel("DIXMAANH")
        x0 = nlp.meta.x0

        @testset "Testing the function" begin
            @test obj(nlp, x0) ≈ DIXMAANHFun(x0)
        end
        @testset "Testing the gradient" begin
            @test norm(grad(nlp, x0) - DIXMAANHGrad(x0)) <= 1e-7
        end
        finalize(nlp)
    end

    @testset "Compared with FlUX AD" begin
        x0 = 2.0 * ones(600)
        @test norm(DIXMAANHGrad(x0) - fDIXMAANHgrad(x0)) <= 1e-11
        for i in 1:10
            x0 = rand(rng, 600)
            @test norm(DIXMAANHGrad(x0) - fDIXMAANHgrad(x0)) <= 1e-11
        end
    end
    @testset verbose = true "CUTEst vs FlUX AD" begin
        nlp = CUTEstModel("DIXMAANH")
        @testset "CUTEst vs FlUX AD with meta x0" begin
            x0 = nlp.meta.x0
            @test norm(grad(nlp, x0) - fDIXMAANHgrad(x0)) <= 1e-7
        end
        finalize(nlp)
    end
end



@testset verbose = true "Testing `DIXMAANI` Function" begin
    # test 1 with CUTEst
    @testset verbose = true "Compared with CUTEst" begin
        nlp = CUTEstModel("DIXMAANI")
        x0 = nlp.meta.x0

        @testset "Testing the function" begin
            @test obj(nlp, x0) ≈ DIXMAANIFun(x0)
        end
        @testset "Testing the gradient" begin
            @test norm(grad(nlp, x0) - DIXMAANIGrad(x0)) <= 1e-7
        end
        finalize(nlp)
    end

    @testset "Compared with FlUX AD" begin
        x0 = 2.0 * ones(600)
        @test norm(DIXMAANIGrad(x0) - fDIXMAANIgrad(x0)) <= 1e-11
        for i in 1:10
            x0 = rand(rng, 600)
            @test norm(DIXMAANIGrad(x0) - fDIXMAANIgrad(x0)) <= 1e-11
        end
    end
    @testset verbose = true "CUTEst vs FlUX AD" begin
        nlp = CUTEstModel("DIXMAANI")
        @testset "CUTEst vs FlUX AD with meta x0" begin
            x0 = nlp.meta.x0
            @test norm(grad(nlp, x0) - fDIXMAANIgrad(x0)) <= 1e-7
        end
        finalize(nlp)
    end
end
