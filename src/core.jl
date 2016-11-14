
"""
    dbootdata{T<:Number}(x::AbstractVector{T}, bi::BootInput)::Vector{Vector{T}}
    dbootdata{T<:Number}(x::Vector{Vector{T}}, bi::BootInput)::Vector{Vector{Vector{T}}}
    dbootdata{T<:Number}(x::Matrix{T}, bi::BootInput)::Vector{Matrix{T}}

Get dependent bootstrap resampled data for the input dataset x using the bootstrap methodology defined in BootInput.
Note, this function should always have output Vector{typeof(x)} as output.
A keyword variant is also available.
"""
dbootdata{T<:Number}(x::AbstractVector{T}, inds::Vector{Vector{Int}})::Vector{Vector{T}} = Vector{T}[ x[inds[n]] for n = 1:length(inds) ]
dbootdata{T<:Number}(x::AbstractVector{T}, bi::BootInput, bm::BootMethod)::Vector{Vector{T}} = dbootdata(x, dbootinds(bi))
function dbootdata{T<:Number}(x::AbstractVector{T}, bi::BootInput, bm::BootTapered)::Vector{Vector{Float64}}
    error("Tapered block bootstrap method is still under construction. If you would like to contribute, please visit the packages github page.")
    error("Tapered block bootstrap currently assumes that input data is of type Float64")
    xOut = dbootdata(apply_influence_function(x), dbootinds(bi))
    dboot_weight!(xOut, bi)
    return(xOut)
end
dbootdata{T<:Number}(x::AbstractVector{T}, bi::BootInput)::Vector{Vector{T}} = dbootdata(x, bi, bi.bootmethod)
#Multivariate methods
function dbootdata{T<:Number}(x::Vector{Vector{T}}, bi::BootInput, bm::BootMethod)::Vector{Vector{Vector{T}}}
    inds = dbootinds(bi)
    return(Vector{Vector{T}}[ dbootdata(x[k], inds) for k = 1:length(x) ])
end
function dbootdata{T<:Number}(x::Matrix{T}, bi::BootInput, bm::BootMethod)::Vector{Matrix{T}}
    inds = dbootinds(bi)
    return(Matrix{T}[ vvtomat(dbootdata(x[:, k], inds)) for k = 1:size(x, 2) ])
end
function dbootdata{T<:Number}(x::Vector{Vector{T}}, bi::BootInput, bm::BootTapered)::Vector{Vector{Vector{T}}}
    error("Tapered block bootstrap method is still under construction. If you would like to contribute, please visit the packages github page.")
end
function dbootdata{T<:Number}(x::Matrix{T}, bi::BootInput, bm::BootTapered)::Vector{Matrix{T}}
    error("Tapered block bootstrap method is still under construction. If you would like to contribute, please visit the packages github page.")
end
dbootdata{T<:Number}(x::Vector{Vector{T}}, bi::BootInput)::Vector{Vector{Vector{T}}} = dbootdata(x, bi, bi.bootmethod)
#Keyword method
function dbootdata(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))
    return(dbootdata(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
end
#Local functions for converting Vector{Vector} to Matrix and back
vvtomat{T}(x::Vector{Vector{T}})::Matrix{T} = T[ x[m][n] for n = 1:data_length(x), m = 1:length(x) ]
mattovv{T}(x::Matrix{T})::Vector{Vector{T}} = Vector{T}[ x[:, m] for m = 1:size(x, 2) ]

"""
    dbootlevel1(data, bi::BootInput)

Get the level 1 bootstrapped statistics associated with dataset in data, and bootstrap methodology in BootInput.
Note, the return type of the output will be determined by bi.flevel1, which must be a function that accepts
the whatever type data is as input. For example, if data is a Vector{T<:Number} then bi.flevel1 might be the
function mean. A more complicated example: if data is Matrix{T<:Number} then bi.flevel1 might be an anonymous
function that implements w'cov(data)w, i.e. a quadratic form over the covariance matrix of the data.
"""
function dbootlevel1(data, bi::BootInput)
    rdata = dbootdata(data, bi)
    try
        level1Out = [ bi.flevel1(rdata[s]) for s = 1:length(rdata) ]
        return(level1Out)
    catch
        error("The function $(bi.flevel1) stored in flevel1 in BootInput does not accept $(typeof(rdata[1])) as input. Please specify an appropriate level1 transformation function.")
    end
    error("Logic fail in dbootlevel1. Please file an issue")
end
function dbootlevel1(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))
    return(dbootlevel1(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
end

"""
    dboot(data, bi::BootInput)

Get the level 2 bootstrapped statistics associated with dataset in data, and bootstrap methodology in BootInput.
Note, the return type of the output will be determined by bi.flevel2, which must be a function that accepts
the output type of dbootlevel1. Note, this will be a vector of the output type of flevel1.

For example, if data is a Vector{T<:Number} and bi.flevel1 is mean, then bi.flevel2 could be var
(the variance function), since the output of mean is likely to be Float64, and var will accept Vector{Float64} as input.

Alternatively, bi.flevel2 could be an anonymous function corresponding to quantile(dbootlevel1_output, [0.025, 0.975]),
in which case the output will be a 95% confidence interval for the test statistic in flevel1.

Note, a keyword argument variant is also available.
"""
function dboot(data, bi::BootInput)
    level1data = dbootlevel1(data, bi)
    try
        level2Out = bi.flevel2(level1data)
        return(level2Out)
    catch
        error("The function $(bi.flevel2) stored in flevel2 in BootInput does not accept $(typeof(level1data)) as input. Please specify an appropriate level2 transformation function.")
    end
    error("Logic fail in dboot. Please file an issue")
end
function dboot(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))
    return(dboot(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
end

"""
    dbootlevel2(data, bi::BootInput)

Core function for obtaining a bootstrapped statistic. The function name dbootlevel2 is included for consistency
with dbootlevel1, but under the hood this function just wraps the main module function dboot. Please see dboot
documentation for more detail.
"""
dbootlevel2(data, bi::BootInput)= dboot(data, bi)
function dbootlevel2(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))
    return(dbootlevel2(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
end



"""
    dbootvar(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                     blmethod::Symbol=:dummy, flevel1::Function=mean, numobsperresample::Number=numobs)

Shortcut method to bootstrap the variance of test statistic implied by flevel1. For example:
dbootvar(data, flevel1=mean) will bootstrap the variance of the sample mean of the dataset in data.
"""
function dbootvar(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, numobsperresample::Number=data_length(data))
    return(dbootvar(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=var, numobsperresample=numobsperresample)))
end

"""
    dbootconf(data ; alpha::Float64=0.05, blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                     blmethod::Symbol=:dummy, flevel1::Function=mean, numobsperresample::Number=numobs)

Shortcut method to bootstrap the 95% confidence interval of the test statistic implied by flevel1. Note, keyword argument
alpha can be used to change the confidence interval, e.g. alpha=0.05 -> 95% confidence interval, alpha = 0.1 -> 90% confidence
interval, etc.

For example:
dbootconf(data, alpha=0.99, flevel1=mean) will bootstrap the a 1% confidence interval for the sample mean of the dataset in data.
"""
function dbootconf(data ; alpha::Float64=0.05, blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, numobsperresample::Number=data_length(data))
    !(0.0 < alpha < 0.5) && error("Invalid alpha for confidence interval. It must lie on the interval (0.0, 0.5)")
    confFunc = (x -> quantile(data, [alpha / 2.0, 1.0 - (alpha / 2.0)]))
    return(dbootvar(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=confFunc, numobsperresample=numobsperresample)))
end
