__precompile__()

"""
Module for dependent bootstrap procedures, by Colin T Bowers

Implemented bootstrap methods: \n
	- IID
	- Stationary
	- Moving Block
	- Circular Block
	- NoOverlapBlock

Implemented block length selection procedures: \n
	- Patton, Politis, and White (2009) Correction to Automatic Block Length Selection For The Dependent Bootstrap

Accepted input dataset types: \n
	- Vector{<:Number}
	- Matrix{<:Number} (where rows are observations and columns are variables)
	- Vector{Vector{<:Number}} (where elements of inner vectors are observations and outer vectors are variables)
	- DataFrame
	- TimeSeries.TimeArray{T,N} (only for N = 1 and N = 2)

Additional input dataset types are easily added. Please open an issue at https://github.com/colintbowers/DependentBootstrap.jl

The module has only a single exported type: \n
	- BootInput  <-- Core input type accepted by all exported functions. Typically constructed via keyword method. See ?BootInput for more detail.

All exported functions exhibit the following keyword signatures: \n
	- exported_func(data, bootinput::BootInput)
	- exported_func(data ; kwargs...)

Most users will be content to use the keyword argument method. In practice, this method wraps a
keyword argument BootInput constructor, which is then input to the exported function BootInput
method. For more detail on accepted keywords, see ?BootInput. All exported functions then use
the input dataset and bootstrap methodology described in BootInput in order to return the
appropriate statistics. A list of exported functions follows: \n
	- optblocklength <-- Estimate the optimal block length for the input dataset
	- dbootinds      <-- Get a vector of bootstrap resampling indices
	- dbootdata      <-- Get a vector of resampled datasets
	- dbootlevel1    <-- Get a vector of level 1 resampled bootstrap statistics
	- dboot          <-- Get the level 2 bootstrap statistic
	- dbootlevel2    <-- Identical to dboot. Included for completeness
	- dbootvar       <-- Wrapper on dboot that sets the level 2 statistic as the variance
	- dbootconf      <-- Wrapper on dboot that sets the level 2 statistic as a confidence interval

I use the phrases level 1 and level 2 statistics in this package in the same manner discussed in
Chapter 1 of Lahiri's textbook Resampling Methods for Dependent Data.

This package has an MIT license. Please see associated LICENSE.md file for more detail.
"""
module DependentBootstrap

using 	StatsBase, Distributions, DataFrames, TimeSeries

import 	Base: 	show

export 	BootInput,
		optblocklength,
		dbootinds,
		dbootdata,
		dbootlevel1,
		dbootlevel2,
		dboot,
		dbootvar,
		dbootconf

const NUM_RESAMPLE = 1000::Int #Default number of bootstrap resamples

#-----------------------------------------------------------------------------------------
# OLD CODE THAT USED Requires.jl TO LOAD DataFrames.jl and TimeSeries.jl
#-----------------------------------------------------------------------------------------
#In order to accommodate different types of datasets, we use the lazy loading
#features offered by the Requires package here, so that the modules in which these
#datasets are defined are not loaded unless they are actually needed. All functions
#that reference these modules are placed here inside __init__() (which is a
#requirement when using Requires package on Julia v0.7+). If you want to add a new
#dataset type to DependentBootstrap, it should be as simple as adding appropriate
#methods here. This module is structured so that nothing else should need to be done.
#The methods that need to be added are:
#    num_obs(data::T)::Int <- Number of observations in dataset
#	 num_var(data::T)::Int <- Number of variables in dataset
#    local_get_var(data::T, i::Int)::Vector <- Get the data associated with the ith variable and output as a Vector
#	 local_get_index(data::T, inds::Vector{Int})::T <- Resample data using resampling indices in inds
#For local_get_index, note that the output type should always match the input type of the dataset.
# function __init__()
# 	@require TimeSeries="9e3dc215-6440-5c97-bce1-76c03772f85e" begin
# 		(num_obs(data::TimeSeries.TimeArray{T,1})::Int) where {T} = size(data, 1)
# 		(num_obs(data::TimeSeries.TimeArray{T,2})::Int) where {T} = size(data, 1)
# 		(num_var(data::TimeSeries.TimeArray{T,1})::Int) where {T} = 1
# 		(num_var(data::TimeSeries.TimeArray{T,2})::Int) where {T} = size(data, 2)
# 		(local_get_var(data::TimeSeries.TimeArray{T,1}, i::Int)::Vector{T}) where {T} = i == 1 ? data.values[:] : error("Invalid index $(i) given data $(typeof(data)) with number of columns $(num_var(data))")
# 		(local_get_var(data::TimeSeries.TimeArray{T,2}, i::Int)::Vector{T}) where {T} = (1 <= i <= num_var(data)) ? data.values[:, i] : error("Invalid index $(i) given data $(typeof(data)) with number of columns $(num_var(data))")
# 		(local_get_index(data::TimeSeries.TimeArray{T,1}, inds::Vector{Int})::TimeSeries.TimeArray{T,1}) where {T} = TimeSeries.TimeArray(TimeSeries.timestamp(data), data.values[inds], TimeSeries.colnames(data) ; unchecked=true)
# 		(local_get_index(data::TimeSeries.TimeArray{T,2}, inds::Vector{Int})::TimeSeries.TimeArray{T,2}) where {T} = TimeSeries.TimeArray(TimeSeries.timestamp(data), data.values[inds, :], TimeSeries.colnames(data) ; unchecked=true)
# 	end
# 	@require DataFrames="a93c6f00-e57d-5684-b7b6-d8193f3e46c0" begin
# 		(num_obs(data::DataFrames.DataFrame)::Int) = size(data, 1)
# 		(num_var(data::DataFrames.DataFrame)::Int) = size(data, 2)
# 		(local_get_var(data::DataFrames.DataFrame, i::Int)) = (1 <= i <= num_var(data)) ? data[:, i] : error("Invalid index $(i) given data $(typeof(data)) with number of columns $(num_var(data))")
# 		(local_get_index(data::DataFrames.DataFrame, inds::Vector{Int})::DataFrames.DataFrame) = data[inds, :]
# 	end
# end
#-----------------------------------------------------------------------------------------

include("types.jl")
include("blocklength.jl")
include("bootinds.jl")
#include("tapered.jl")
include("core.jl")

end # module
