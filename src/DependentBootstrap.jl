module DependentBootstrap
#-----------------------------------------------------------
#PURPOSE
#	Colin T. Bowers module for dependent bootstraps
#NOTES
#	Valid input datasets are currently:
#		-AbstractVector{T<:Number}
#		-Vector{Vector{T<:Number}}
#		-AbstractMatrix{T<:Number}
#	To extend the package to a new type of dataset, you will need to add new methods to the following functions:
#		data_length (in types.jl)
#		dbootdata (in core.jl) (where output must be of type Vector{typeof(inputdata)})
#		optblocklength (in blocklength.jl)
#		several functions in tapered.jl IF you want to use tapered block bootstrap (which is actually not operational yet)
#-----------------------------------------------------------

using 	StatsBase, Distributions

import 	Base: 	show

export 	BootMethod,
		BLMethod,
		BootStationary,
		BootMoving,
		BootCircular,
		BootNoOverlap,
		BootTapered,
		BLPPW2009,
		BLPP2002,
		BootInput,
		setblocklength!,
		setnumresample!,
		setflevel1!,
		setflevel2!,
		optblocklength,
		dbootinds,
		dbootdata,
		dbootlevel1,
		dbootlevel2,
		dboot,
		dbootvar,
		dbootconf

const NUM_RESAMPLE = 1000::Int #Default value
const BOOT_METHOD_TYPES = Symbol[:BootIID, :BootStationary, :BootMoving, :BootNoOverlap, :BootCircular, :BootTapered]::Vector{Symbol}
const BOOT_METHODS = Symbol[:iid, :stationary, :moving, :nooverlap, :circular, :tapered]::Vector{Symbol}
const BLOCK_LENGTH_METHOD_TYPES = Symbol[:BLPP2002, :BLPPW2009]::Vector{Symbol}
const BLOCK_LENGTH_METHODS = Symbol[:pp2002, :ppw2009]::Vector{Symbol}
const BANDWIDTH_METHODS = Symbol[:p2003]::Vector{Symbol}

include("types.jl")
include("common.jl")
include("blocklength.jl")
include("bootinds.jl")
include("tapered.jl")
include("core.jl")

end # module
