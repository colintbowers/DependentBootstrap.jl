"""
Module for dependent bootstrap procedures, by Colin T Bowers \n
\n
Implemented bootstrap methods: \n
	- IID \n
	- Stationary \n
	- Moving Block \n
	- Circular Block \n
	- NoOverlapBlock \n
\n
Implemented block length selection procedures: \n
	- Patton, Politis, and White (2009) Correction to Automatic Block Length Selection For The Dependent Bootstrap \n
	- Paparoditis and Politis (2002) The Tapered Block Bootstrap For General Statistics From Stationary Sequences \n
\n
Accepted input dataset types: \n
	- Vector{<:Number}
	- Matrix{<:Number}
	- Vector{Vector{<:Number}}
\n
Additional input dataset types are easily added. Please open an issue on at https://github.com/colintbowers/DependentBootstrap.jl
\n
This package has an MIT license. Please see associated LICENSE.md file for more detail.
"""
module DependentBootstrap

using 	Requires
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
