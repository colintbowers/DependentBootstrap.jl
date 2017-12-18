
"""
    apply_inds_to_data(data::T, inds::Vector{Int})::T
    apply_inds_to_data(data::T, inds::Vector{Vector{Int}})::Vector{T}

Apply the output from a call to dbootinds_one to data to obtain a single,
resampled dataset. \n

NOTE 1: If adding a new type for data to this module, this function must be
extended to the new type. Current accepted data types are: \n
    - Vector{<:Number} \n
    - Matrix{<:Number} \n
    - Vector{Vector{<:Number}} \n

NOTE 2: I deliberately do not build views into data, since the level 1 bootstrap
function is typically able to take advantage of BLAS routines (eg mean, var, etc),
and in these cases it will be more efficient to build the new, resampled data
sequentially in memory.
"""
(apply_inds_to_data(data::Vector{T}, inds::Vector{Int})::Vector{T}) where {T<:Number} = data[inds]
(apply_inds_to_data(data::Matrix{T}, inds::Vector{Int})::Matrix{T}) where {T<:Number} = data[inds, :]
(apply_inds_to_data(data::Vector{Vector{T}}, inds::Vector{Int})::Vector{Vector{T}}) where {T<:Number} = [ y[inds] for y in data ]

"""
    dbootdata_one(data::T, bi::BootInput)::T \n
    dbootdata_one(data::T; kwargs...)::T \n

Get a single resampled dataset of the input data using the dependent boostrap
methodology defined in BootInput. \n
Note, the output type will always be the same as the type of the input data.
"""
(dbootdata_one(data::TD, bi::BootInput, bm::TM)::TD) where {TD, TM<:BootMethod} = apply_inds_to_data(data, dbootinds_one(bi))
function dbootdata_one_infl(data_infl::TD, bi::BootInput, bm::BootTapered)::TD where {TD}
    error("The tapered block bootstrap is currently not available in this package. Users interested in contributing should check the package github page.")
    dataout = apply_inds_to_data(data_infl, dbootinds_one(bi))
    dboot_weight!(dataout, bi)
    return dataout
end
(dbootdata_one(data::T, bi::BootInput, bm::BootTapered)::T) where {T} = dbootdata_one_infl(apply_influence_function(data), bi, bm)
(dbootdata_one(data::T, bi::BootInput)::T) where {T} = dbootdata_one(data, bi, bi.bootmethod)

"""
dbootdata(data::T, bi::BootInput)::Vector{T} \n
dbootdata(data::T; kwargs...)::Vector{T} \n

Get the resampled datasets of the input data using the dependent bootstrap
methodology defined in BootInput. \n
There is also a keyword variant that calls the keyword constructor for BootInput. \n
Note, this function should always have output type Vector{T}.
"""
(dbootdata(data::TD, bi::BootInput, bm::TM)::Vector{T}) where {TD, TM<:BootMethod} = [ apply_inds_to_data(data, dbootinds_one(bi)) for j = 1:bi.numresample ]
function dbootdata(data::T, bi::BootInput, bm::BootTapered)::Vector{T} where {T}
    data_infl = apply_influence_function(data)
    return [ dbootdata_one_infl(data_infl, bi, bm) for j = 1:bi.numresample ]
end
(dbootdata(data::T, bi::BootInput)::Vector{T}) where {T} = dbootdata(data, bi, bi.bootmethod)
function dbootdata(data::T ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                 blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))::T where {T}
    return dbootdata(data, BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample))
end


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
dbootlevel1(data, bi::BootInput) = [ bi.flevel1(dbootdata_one(data, bi)) for j = 1:bi.numresample ]
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
dboot(data, bi::BootInput) = bi.flevel2(dbootlevel1(data, bi))
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
