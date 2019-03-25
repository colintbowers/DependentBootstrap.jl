

"num_obs <- Internal function used to determine the number of observations in the input dataset"
(num_obs(data::AbstractVector{T})::Int) where {T} = length(data)
(num_obs(data::AbstractMatrix{T})::Int) where {T} = size(data, 1)
(num_obs(data::Vector{Vector{T}})::Int) where {T} = (isempty(data) || any(length.(data) .!= length(data[1]))) ? error("Input dataset is empty, or inner vectors of input dataset do not have matching length: $(length.(data))") : length(data[1])
(num_obs(data::DataFrames.DataFrame)::Int) = size(data, 1)
(num_obs(data::TimeSeries.TimeArray{T,1})::Int) where {T} = size(data, 1)
(num_obs(data::TimeSeries.TimeArray{T,2})::Int) where {T} = size(data, 1)

"num_var <- Internal function used to determine the number of variables in the input dataset"
(num_var(data::AbstractVector{T})::Int) where {T} = 1
(num_var(data::AbstractMatrix{T})::Int) where {T} = size(data, 2)
(num_var(data::Vector{Vector{T}})::Int) where {T} = isempty(data) ? error("Input dataset is empty") : length(data)
(num_var(data::DataFrames.DataFrame)::Int) = size(data, 2)
(num_var(data::TimeSeries.TimeArray{T,1})::Int) where {T} = 1
(num_var(data::TimeSeries.TimeArray{T,2})::Int) where {T} = size(data, 2)

"local_get_var <- Internal function used to get the ith variable in dataset"
(local_get_var(data::AbstractVector{T}, i::Int)::Vector{T}) where {T} = i == 1 ? data[:] : error("Invalid index $(i) given data $(typeof(data))")
(local_get_var(data::AbstractMatrix{T}, i::Int)::Vector{T}) where {T} = (1 <= i <= size(data,2)) ? data[:, i] : error("Invalid index $(i) given data $(typeof(data)) with number of columns $(size(data, 2))")
(local_get_var(data::Vector{Vector{T}}, i::Int)::Vector{T}) where {T} = (1 <= i <= length(data)) ? data[i] : error("Invalid index $(i) given data $(typeof(data)) with outer length: $(length(data))")
(local_get_var(data::DataFrames.DataFrame, i::Int)) = (1 <= i <= num_var(data)) ? data[:, i] : error("Invalid index $(i) given data $(typeof(data)) with number of columns $(num_var(data))")
(local_get_var(data::TimeSeries.TimeArray{T,1}, i::Int)::Vector{T}) where {T} = i == 1 ? values(data) : error("Invalid index $(i) given data $(typeof(data)) with number of columns $(num_var(data))")
(local_get_var(data::TimeSeries.TimeArray{T,2}, i::Int)::Vector{T}) where {T} = (1 <= i <= num_var(data)) ? values(data)[:,i] : error("Invalid index $(i) given data $(typeof(data)) with number of columns $(num_var(data))")

"local_get_index <- Internal function used to resample the dataset data using the input resampling index inds"
(local_get_index(data::AbstractVector{T}, inds::Vector{Int})::Vector{T}) where {T} = data[inds]
(local_get_index(data::AbstractMatrix{T}, inds::Vector{Int})::Matrix{T}) where {T} = data[inds, :]
(local_get_index(data::Vector{Vector{T}}, inds::Vector{Int})::Vector{Vector{T}}) where {T} = [ y[inds] for y in data ]
(local_get_index(data::DataFrames.DataFrame, inds::Vector{Int})::DataFrames.DataFrame) = data[inds, :]
(local_get_index(data::TimeSeries.TimeArray{T,1}, inds::Vector{Int})::TimeSeries.TimeArray{T,1}) where {T} = TimeSeries.TimeArray(TimeSeries.timestamp(data), values(data)[inds], TimeSeries.colnames(data) ; unchecked=true)
(local_get_index(data::TimeSeries.TimeArray{T,2}, inds::Vector{Int})::TimeSeries.TimeArray{T,2}) where {T} = TimeSeries.TimeArray(TimeSeries.timestamp(data), values(data)[inds, :], TimeSeries.colnames(data) ; unchecked=true)


"""
    dbootdata_one(data::T, bi::BootInput)::T
    dbootdata_one(data::T; kwargs...)::T

Get a single resampled dataset of the input data using the dependent boostrap
methodology defined in BootInput.

A keyword method that calls the keyword constructor for BootInput is also provided. Please
use ?BootInput at the REPL for more detail on feasible keywords.

Note, the output type will always be the same as the type of the input data.
"""
(dbootdata_one(data::Td, bi::BootInput)::Td) where {Td} = local_get_index(data, dbootinds_one(bi))
(dbootdata_one(data::Td ; kwargs...)::Td) where {Td} = dbootdata_one(data, BootInput(data ; kwargs...))
# function dbootdata_one_infl(data_infl::Td, bi::BootInput{BootTapered})::TD where {Td}
#     dataout = apply_inds_to_data(data_infl, dbootinds_one(bi))
#     dboot_weight!(dataout, bi)
#     return dataout
# end
# function (dbootdata_one(data::Td, bi::BootInput{BootTapered})::Td) where {Td}
#     return dbootdata_one_infl(apply_influence_function(data), bi, bm)
# end

"""
    dbootdata(data::T , bi::BootInput)::Vector{T}
    dbootdata(data::T ; kwargs...)::Vector{T}

Get the resampled datasets of the input data using the dependent bootstrap
methodology defined in BootInput.

A keyword method that calls the keyword constructor for BootInput is also provided. Please
use ?BootInput at the REPL for more detail on feasible keywords.

Note, this function should always have output type Vector{T}.
"""
(dbootdata(data::Td, bi::BootInput)::Vector{Td}) where {Td} = [ local_get_index(data, dbootinds_one(bi)) for j = 1:bi.numresample ]
(dbootdata(data::Td ; kwargs...)::Vector{Td}) where {Td} = dbootdata(data, BootInput(data ; kwargs...))
# function dbootdata(data::Td, bi::BootInput{BootTapered})::Vector{Td} where {Td}
#     data_infl = apply_influence_function(data)
#     return [ dbootdata_one_infl(data_infl, bi, bm) for j = 1:bi.numresample ]
# end

"""
    dbootlevel1(data::T1, bi::BootInput)
    dbootlevel1(data::T1; kwargs...)

Get the level 1 bootstrapped statistics associated with dataset in data, and bootstrap methodology in BootInput.

A keyword method that calls the keyword constructor for BootInput is also provided. Please
use ?BootInput at the REPL for more detail on feasible keywords.

Note, the return type is determined by bi.flevel1, which must be a function that accepts T1,
ie typeof(data), as input. It may return any output type T2, as long as bi.flevel2 will
accept Vector{T2} as input.

For example, if data is a Vector{<:Number} then bi.flevel1 might be the function mean,
which in this case will return Float64, so bi.flevel2 must be some function that can
accept Vector{Float64} as input.

A more complicated example: if data is Matrix{<:Number} then bi.flevel1 might be the anonymous
function x->mean(x,dims=1), which in this case will return a single row Matrix{Float64}, and
so bi.flevel2 must be some function that can accept Vector{Matrix{Float64}} as input.
"""
dbootlevel1(data, bi::BootInput) = [ bi.flevel1(dbootdata_one(data, bi)) for j = 1:bi.numresample ]
dbootlevel1(data ; kwargs...) = dbootlevel1(data, BootInput(data ; kwargs...))

"""
    dboot(data, bi::BootInput)
    dboot(data ; kwargs...)

Get the level 2 bootstrapped statistics associated with dataset in data, and bootstrap methodology in BootInput.

A keyword method that calls the keyword constructor for BootInput is also provided. Please
use ?BootInput at the REPL for more detail on feasible keywords.

Note, the return type of the output will be determined by bi.flevel2, which must be a function that accepts
Vector{T}, where T is the output type of bi.flevel1.

For example, if data is a Vector{<:Number} and bi.flevel1 is mean, then in this case, bi.flevel1 will return
Float64, and so bi.flevel2 must be some function that accepts Vector{Float64} as input (and can have any output
type.)

Alternatively, bi.flevel2 could be the anonymous function (x -> quantile(x, [0.025, 0.975])), in which case
the input should be Vector{Float64}, and so bi.flevel1 should return Float64. Note, the output of bi.flevel2
in this case will be a 2-element Vector{Float64} with elements corresponding bootstrapped 95% confidence interval
for the level1 statistic of the input dataset
"""
dboot(data, bi::BootInput) = bi.flevel2(dbootlevel1(data, bi))
dboot(data ; kwargs...) = dboot(data, BootInput(data ; kwargs...))

"dbootlevel2 <- Identical to the dboot function. This function is only included for naming consistency with dbootlevel1"
dbootlevel2(data, bi::BootInput)= dboot(data, bi)
dbootlevel2(data ; kwargs...)= dbootlevel2(data, BootInput(data ; kwargs...))

"dbootvar <- Identical to dboot but with the level 2 statistic set to variance"
dbootvar(data ; kwargs...) = dboot(data, BootInput(data ; flevel2=var, kwargs...))

"dbootconf <- Identical to dboot but with the level 2 statistic set to a confidence interval with width determined by keyword alpha. Default alpha=0.05 corresponds to a 95% confidence interval."
dbootconf(data ; alpha::Float64=0.05, kwargs...) = (0.0 < alpha < 0.5) ? dboot(data, BootInput(data ; flevel2=(x -> quantile(x, [alpha/2, 1-(alpha/2)])), kwargs...)) : error("Invalid alpha of $(alpha) for confidence interval. alpha must lie on the interval (0.0, 0.5)")
