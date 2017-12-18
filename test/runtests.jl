
using Base.Test
using StatsBase
using Distributions
using DependentBootstrap

#Quick and dirty function for simulating AR(1) data with deterministic start point for random number generation
function temp_ar(seedInt::Int)
    srand(seedInt)
    N = 100
    e = randn(N)
    x = Array{Float64}(N)
    x[1] = 0.0
    for n = 2:N
        x[n] = 0.8 * x[n-1] + e[n]
    end
    return(x)
end

#Get the AR(1) data
x = temp_ar(1234);

#Bootstrap methods
bmVec = Symbol[:iid, :stationary, :moving, :circular, :nooverlap];

#Expected results for each bootstrap method given seed 1234
correctBlockLength = Float64[8.015752150100226, 7.002404488495542, 8.015752150100226, 8.015752150100226, 8.015752150100226];
correctLevel2 = Float64[0.021527103583357056, 0.07168871559296867, 0.07330809341723817, 0.07069130008694534, 0.07327722639804433];

#Check optimal block length and bootstrapped variance of mean for every bootstrap method
@testset "Univariate bootstrap block length and method tests" begin
    for k in eachindex(bmVec)
        println("Testing $(bmVec[k]) (univariate)")
        bi = BootInput(x, bootmethod=bmVec[k]);
        @test isapprox(bi.blocklength, correctBlockLength[k])
        setblocklength!(bi, 5.0);
        srand(1234);
        @test isapprox(dboot(x, bi), correctLevel2[k])
    end
end

#Check that multivariate bootstrapping works too (only bother for stationary)
x = Vector{Float64}[ temp_ar(n) for n = 1:4 ];
@testset "Multivariate bootstrap block length and method tests" begin
    println("Testing stationary (multivariate)")
    bi = BootInput(x, bootmethod=:stationary);
    @test isapprox(bi.blocklength, 9.375963199884655)
    setblocklength!(bi, 5.0);
    fl1 = (x -> maximum(Float64[ mean(x[k]) for k = 1:length(x) ])); #Maximum of sample means
    fl2 = (x -> quantile(x, [0.025, 0.975])); #5% confidence interval
    setflevel1!(bi, fl1);
    setflevel2!(bi, fl2);
    srand(1234)
    y1 = dboot(x, bi)
    @test isapprox(y1[1], -0.047408575390325614)
    @test isapprox(y1[2], 1.0913695579833216)
    xMat = Array{Float64}(length(x[1]), length(x))
    for k = 1:length(x)
        xMat[:, k] = x[k]
    end
    srand(1234)
    y2 = dboot(x, bi)
    @test isapprox(y1[1], y2[1])
    @test isapprox(y1[2], y2[2])
end
