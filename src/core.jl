
"""
    dbootdata(data, bi::BootInput) \n
    dbootdata(data; kwargs...) \n

Get dependent bootstrap resampled data for the input dataset using the bootstrap methodology defined in BootInput. \n
There is also a keyword variant that calls the keyword constructor for BootInput. \n
Note, this function should always have output type Vector{typeof(x)}.
"""
(dbootdata(x::AbstractVector{T}, inds_vv::Vector{Vector{Int}})::Vector{Vector{T}}) where {T<:Number} = [ x[inds] for inds in inds_vv ]
(dbootdata(x::AbstractVector{T}, bi::BootInput, bm::BootMethod)::Vector{Vector{T}}) where {T<:Number} = dbootdata(x, dbootinds(bi))
function dbootdata(x::AbstractVector{T}, bi::BootInput, bm::BootTapered)::Vector{Vector{Float64}} where {T<:Number}
    error("Tapered block bootstrap method is still under construction. If you would like to contribute, please visit the packages github page.")
    error("Tapered block bootstrap currently assumes that input data is of type Float64")
    xOut = dbootdata(apply_influence_function(x), dbootinds(bi))
    dboot_weight!(xOut, bi)
    return(xOut)
end
#Multivariate methods
function dbootdata(x::Vector{Vector{T}}, bi::BootInput, bm::BootMethod)::Vector{Vector{Vector{T}}} where {T<:Number}
    inds = dbootinds(bi)
    return Vector{Vector{T}}[ Vector{T}[ x[j][indsk] for j = 1:length(x) ] for indsk in inds ]
end
function dbootdata(x::Matrix{T}, bi::BootInput, bm::BootMethod)::Vector{Matrix{T}} where {T<:Number}
    inds = dbootinds(bi)
    return(Matrix{T}[ x[indsk, :] for indsk in inds ])
end
function dbootdata(x::Vector{Vector{T}}, bi::BootInput, bm::BootTapered)::Vector{Vector{Vector{T}}} where {T<:Number}
    error("Tapered block bootstrap method is still under construction. If you would like to contribute, please visit the packages github page.")
end
function dbootdata(x::Matrix{T}, bi::BootInput, bm::BootTapered)::Vector{Matrix{T}} where {T<:Number}
    error("Tapered block bootstrap method is still under construction. If you would like to contribute, please visit the packages github page.")
end
#Core method
dbootdata(x, bi::BootInput) = dbootdata(x, bi, bi.bootmethod)
#Keyword method
function dbootdata(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))
    return (dbootdata(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
end
#Local functions for converting Vector{Vector} to Matrix and back
(vvtomat(x::Vector{Vector{T}})::Matrix{T}) where {T} = T[ x[m][n] for n = 1:data_length(x), m = 1:length(x) ]
(mattovv(x::Matrix{T})::Vector{Vector{T}}) where {T} = Vector{T}[ x[:, m] for m = 1:size(x, 2) ]

"""
    dbootlevel1(data, bi::BootInput) \n
    dbootlevel1(data; kwargs...) \n

Get the level 1 bootstrapped statistics associated with dataset in data, and bootstrap methodology in BootInput. \n

A keyword method that calls the BootInput keyword constructor is also accepted. \n

Note, the return type is determined by bi.flevel1, which must be a function that accepts typeof(data)
as input. \n

For example, if data is a Vector{T<:Number} then bi.flevel1 might be the function mean \n

A more complicated example: if data is Matrix{T<:Number} then bi.flevel1 might be an anonymous
function that implements w'cov(data)w, i.e. a quadratic form over the covariance matrix of the data
"""
function dbootlevel1(data, bi::BootInput)
    rdata = dbootdata(data, bi)
    try
        lvl1out = [ bi.flevel1(y) for y in rdata ]
        return lvl1out
    catch
        error("The function $(bi.flevel1) stored in flevel1 in BootInput does not accept $(typeof(rdata[1])) as input. Please specify an appropriate level1 transformation function.")
    end
    error("Logic fail in dbootlevel1. Please file an issue")
end
function dbootlevel1(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))
    return (dbootlevel1(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
end

"""
    dboot(data, bi::BootInput) \n
    dboot(data ; kwargs...) \n

Get the level 2 bootstrapped statistics associated with dataset in data, and bootstrap methodology in BootInput. \n
Note, the return type of the output will be determined by bi.flevel2, which must be a function that accepts
Vector{T}, where T is the output type of bi.flevel1.

For example, if data is a Vector{T<:Number} and bi.flevel1 is mean, then bi.flevel2 could be var
(the variance function), since the output of mean is likely to be Float64, and var will accept Vector{Float64} as input.
Alternatively, bi.flevel2 could be the anonymous function (x -> quantile(x, [0.025, 0.975])), in which case
the output will be a bootstrapped 95% confidence interval for the level1 statistic of the input dataset \n

A keyword method that calls the keyword constructor of BootInput is also available.
"""
function dboot(data, bi::BootInput)
    lvl1data = dbootlevel1(data, bi)
    try
        lvl2out = bi.flevel2(lvl1data)
        return lvl2out
    catch
        error("The function $(bi.flevel2) stored in flevel2 in BootInput does not accept $(typeof(lvl1data)) as input. Please specify an appropriate level2 transformation function.")
    end
    error("Logic fail in dboot. Please file an issue")
end
function dboot(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))
    return(dboot(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
end

"dbootlevel2 <- Identical to the dboot function. This function is only included for naming consistency with dbootlevel1"
dbootlevel2(data, bi::BootInput)= dboot(data, bi)
function dbootlevel2(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))
    return (dbootlevel2(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
end

"dbootvar <- Identical to dboot but with the level 2 statistic set to variance"
function dbootvar(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, numobsperresample::Number=data_length(data))
    return(dbootvar(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=var, numobsperresample=numobsperresample)))
end

"dbootconf <- Identical to dboot but with the level 2 statistic set to a confidence interval with width determined by keyword alpha. Default alpha=0.05 corresponds to a 95% confidence interval."
function dbootconf(data ; alpha::Float64=0.05, blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, numobsperresample::Number=data_length(data))
    !(0.0 < alpha < 0.5) && error("Invalid alpha for confidence interval. It must lie on the interval (0.0, 0.5)")
    confFunc = (x -> quantile(data, [alpha / 2.0, 1.0 - (alpha / 2.0)]))
    return(dbootvar(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=confFunc, numobsperresample=numobsperresample)))
end
