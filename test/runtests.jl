
using Random, Test
using StatsBase
using Distributions
using DependentBootstrap
using DataFrames, TimeSeries

#Quick and dirty function for simulating AR(1) data with deterministic start point for random number generation
function temp_ar(seedint::Int)
    #srand(seedInt)
    Random.seed!(seedint)
    N = 100
    e = randn(N)
    x = NaN*ones(Float64,N)
    x[1] = 0.0
    for n = 2:N
        x[n] = 0.8 * x[n-1] + e[n]
    end
    return x
end

#Get the AR(1) data
x = temp_ar(1234);
xmat = hcat(x, temp_ar(5678));
xvv = [ xmat[:, k] for k = 1:size(xmat,2) ]

#bootstrap methods
bootmethodvec = [:iid, :stationary, :moving, :circular, :nooverlap];
bootmethodtypevec = [DependentBootstrap.BootIID(), DependentBootstrap.BootStationary(), DependentBootstrap.BootMoving(),
                     DependentBootstrap.BootCircular(), DependentBootstrap.BootNoOverlap()]
blocklengthvec = [0.0, 5.0];
blocklengthmethodvec = [:ppw2009]
blocklengthmethodtypevec = [DependentBootstrap.BLPPW2009()]

#Test constructor
@testset "BootInput constructor test" begin
    for kbm = 1:length(bootmethodvec)
        for kbl = 1:length(blocklengthvec)
            for kblm = 1:length(blocklengthmethodvec)
                bi = BootInput(x, 100, blocklengthvec[kbl], 200, bootmethodtypevec[kbm],
                               blocklengthmethodtypevec[kblm], var, std, 300, mean)
                @test bi.numresample == 200
                @test bi.bootmethod == bootmethodtypevec[kbm]
                @test bi.blocklengthmethod == blocklengthmethodtypevec[kblm]
                @test bi.flevel1 == var
                @test bi.flevel2 == std
                @test bi.numobsperresample == 300
                @test bi.fblocklengthcombine == mean
                bi = BootInput(x, numresample=200, bootmethod=bootmethodvec[kbm], blocklength=blocklengthvec[kbl],
                               blocklengthmethod=blocklengthmethodvec[kblm], flevel1=var, flevel2=std,
                               numobsperresample=300, fblocklengthcombine=mean)
                @test bi.numresample == 200
                @test bi.bootmethod == bootmethodtypevec[kbm]
                @test bi.blocklengthmethod == blocklengthmethodtypevec[kblm]
                @test bi.flevel1 == var
                @test bi.flevel2 == std
                @test bi.numobsperresample == 300
                @test bi.fblocklengthcombine == mean
            end
        end
    end
end

bootmethodvec = Symbol[:iid, :stationary, :moving]
correctblocklength1 = Vector{Float64}[Float64[1.0, 7.002404488495543, 8.015752150100226]];
correctblocklength2 = Vector{Float64}[Float64[1.0, 6.837001835299588, 7.826413377230708]];

#Test block length selection procedure
@testset "Block length selection test" begin
    for kbm = 1:length(bootmethodvec)
        for kblm = 1:length(blocklengthmethodvec)
            bi = BootInput(xmat[:, 1], numresample=200, bootmethod=bootmethodvec[kbm], blocklengthmethod=blocklengthmethodvec[kblm])
            @test isapprox(bi.blocklength, correctblocklength1[kblm][kbm])
            bi = BootInput(xmat[:, 2], numresample=200, bootmethod=bootmethodvec[kbm], blocklengthmethod=blocklengthmethodvec[kblm])
            @test isapprox(bi.blocklength, correctblocklength2[kblm][kbm])
            bi = BootInput(xmat, numresample=200, bootmethod=bootmethodvec[kbm], blocklengthmethod=blocklengthmethodvec[kblm], fblocklengthcombine=mean)
            @test isapprox(bi.blocklength, mean([correctblocklength1[kblm][kbm],correctblocklength2[kblm][kbm]]))
        end
    end
end

#Test univariate bootstrap method
bootmethodvec = [:iid, :stationary, :moving, :circular, :nooverlap];
correctbootunivariate = Float64[0.02231785457420185,0.0805280736777265,0.07755558515710691,0.07248357346968781,0.0718412998084052]
correctbootunivariatebl1 = 0.021374650187840242*ones(Float64, length(bootmethodvec))
@testset "Univariate bootstrap test" begin
    for kbm = 1:length(bootmethodvec)
        Random.seed!(1234)
        y = dboot(x, numresample=1000, bootmethod=bootmethodvec[kbm], blocklength=5, flevel1=mean, flevel2=var)
        @test isapprox(y, correctbootunivariate[kbm])
    end
    for kbm = 1:length(bootmethodvec)
        Random.seed!(1234)
        y = dboot(x, numresample=500, bootmethod=bootmethodvec[kbm], blocklength=1, flevel1=mean, flevel2=var)
        @test isapprox(y, correctbootunivariatebl1[kbm])
    end
end

#Test multivariate bootstrap
bootmethodvec = [:iid, :stationary, :moving, :circular, :nooverlap];
correctbootmultmatrix = Float64[0.19194008990823777,0.40101052326568176,0.5186276897440075,0.44038974430819566,0.3960765922525029]
correctbootmultvv = Float64[0.1461957546031049,0.3546585480604615,0.316914521186327,0.34363235435327766,0.361937416083102]
@testset "Multivariate bootstrap test" begin
    for kbm = 1:length(bootmethodvec)
        Random.seed!(1234)
        y = dbootvar(xmat, numresample=1000, bootmethod=bootmethodvec[kbm], blocklength=5, flevel1=x->minimum(x[:,1].*x[:,2]))
        @test isapprox(y, correctbootmultmatrix[kbm])
        Random.seed!(1234)
        y = dbootconf(xvv, numresample=1000, bootmethod=bootmethodvec[kbm], blocklength=5, flevel1=flevel1=x->mean([mean(x[1]),mean(x[2])]))
        @test isapprox(y[2], correctbootmultvv[kbm])
    end
end

#Test dbootdata and dbootinds
@testset "Test dbootdata and dbootinds" begin
    a = dbootinds(x, numresample=20, blocklength=5, bootmethod=:stationary)
    @test length(a) == 20
    a = dbootdata(x, numresample=20, blocklength=5, bootmethod=:stationary)
    @test length(a) == 20
end

#Test exotic dataset types
@testset "Exotic data types test" begin
    Random.seed!(1234)
    xdf = DataFrame(xmat)
    #y = dboot(xdf, bootmethod=:stationary, blocklength=5, numresample=1000, flevel1=x->mean(DataFrames.columns(x)[1]))
    y = dboot(xdf, bootmethod=:stationary, blocklength=5, numresample=1000, flevel1=(x->mean(x[:,1])))
    @test isapprox(y, 0.0805280736777265)
    dtvec = [ Date(2000)+Day(n) for n = 1:size(xmat,1) ]
    xta1 = TimeSeries.TimeArray(dtvec, x)
    xta2 = TimeSeries.TimeArray(dtvec, xmat)
    Random.seed!(1234)
    y1 = dboot(xta1, bootmethod=:stationary, blocklength=5, numresample=1000, flevel1=x->mean(values(x)))
    Random.seed!(1234)
    y2 = dboot(xta2, bootmethod=:stationary, blocklength=5, numresample=1000, flevel1=x->mean(values(x)))
    @test isapprox(y1, 0.0805280736777265)
    @test isapprox(y2, 0.05077003035163576)
end
