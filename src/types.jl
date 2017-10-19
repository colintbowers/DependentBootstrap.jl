
#Two abstract types to nest all bootstrap method types and block length selection method types
abstract type BootMethod end
abstract type BLMethod end

#Bootstrap method type definitions
mutable struct BootIID <: BootMethod ; end
mutable struct BootStationary <: BootMethod ; end
mutable struct BootMoving <: BootMethod ; end
mutable struct BootNoOverlap <: BootMethod ; end
mutable struct BootCircular <: BootMethod ; end
mutable struct BootTapered <: BootMethod
	kernelfunction::Symbol
	function BootTapered(kernelfunction::Symbol)
		!(kernelfunction == :trap || kernelfunction == :smooth) && error("Invalid kernel function symbol of $(kernelfunction) for tapered block bootstrap.")
		new(kernelfunction)
	end
end

#Block length method type definitions
mutable struct BLPP2002 <: BLMethod
	bootmethod::Symbol
	bandwidthmethod::Symbol
	kernelfunction::Symbol
	function BLPP2002(bootmethod::Symbol, bandwidthmethod::Symbol, kernelfunction::Symbol)
		!(any(bootmethod .== BOOT_METHODS)) && error("Invalid bootstrap method of $(bootmethod)")
		!(any(bandwidthmethod .== BANDWIDTH_METHODS)) && error("Invalid bandwidth method of $(bandwidthmethod)")
		!(kernelfunction == :trap || kernelfunction == :smooth) && error("Invalid kernel function symbol of $(kernelfunction) for tapered block bootstrap.")
		new(bootmethod, bandwidthmethod, kernelfunction)
	end
end
mutable struct BLPPW2009 <: BLMethod
	bootmethod::Symbol
	bandwidthmethod::Symbol
	function BLPPW2009( bootmethod::Symbol, bandwidthmethod::Symbol)
		!(any(bootmethod .== BOOT_METHODS)) && error("Invalid bootstrap method of $(bootmethod)")
		!(any(bandwidthmethod .== BANDWIDTH_METHODS)) && error("Invalid bandwidth method of $(bandwidthmethod)")
		new(bootmethod, bandwidthmethod)
	end
end

#Shortened symbol representation of each type
#Note, these should be consistent with values in module constants BOOT_METHODS and BLOCK_LENGTH_METHODS
shortsymbol(bm::BootIID) = :iid
shortsymbol(bm::BootStationary) = :stationary
shortsymbol(bm::BootMoving) = :moving
shortsymbol(bm::BootNoOverlap) = :nooverlap
shortsymbol(bm::BootCircular) = :circular
shortsymbol(bm::BootTapered) = :tapered
shortsymbol(blm::BLPP2002) = :pp2002
shortsymbol(blm::BLPPW2009) = :ppw2009

#Bootstrap method type constructors
BootTapered() = BootTapered(:trap)

#Block length type constructors
BLPP2002() = BLPP2002(:stationary, :p2003, :trap)
BLPP2002(bootmethod::Symbol) = BLPP2002(bootmethod, :p2003, :trap)
BLPPW2009() = BLPPW2009(:stationary, :p2003)
BLPPW2009(bootmethod::Symbol) = BLPPW2009(bootmethod, :p2003)

#Get type based on symbol input
function symboltobootmethod(s::Symbol)::BootMethod
    s == :iid && return (BootIID())
    s == :stationary && return (BootStationary())
    s == :moving && return (BootMoving())
    s == :nooverlap && return (BootNoOverlap())
    s == :circular && return (BootCircular())
    s == :tapered && return (BootTapered())
    error("Invalid input symbol of $(s)")
end
function symboltoblmethod(s::Symbol)::BLMethod
    s == :pp2002 && return (BLPP2002())
    s == :ppw2009 && return (BLPPW2009())
    error("Invalid input symbol of $(s)")
end

#Get appropriate block length method given the bootstrap method
bootmethodtoblmethod(bm::BootMethod)::BLMethod = BLPPW2009(shortsymbol(bm), :p2003)
bootmethodtoblmethod(bm::BootTapered)::BLMethod = BLPP2002(shortsymbol(bm), :p2003, :trap)
bootmethodtoblmethod(s::Symbol)::BLMethod = bootmethodtoblmethod(symboltobootmethod(s))


#BootInput type definition. Note, almost all core functions accept this as input
mutable struct BootInput
	numobs::Int #Number of observations in the dataset that is to be bootstrapped
    blocklength::Float64 #Block length for bootstrapping
    numresample::Int #Number of resamples
	bootmethod::BootMethod #Bootstrap method
	blmethod::BLMethod #Block length detection method
	flevel1::Function #The estimator of interest
	flevel2::Function #The distributional parameter of flevel1 that we want to bootstrap
    numobsperresample::Int #Number of observations to be drawn (with replacement) per resample (almost always equal to numobs)
	function BootInput(numobs::Int, blocklength::Float64, numresample::Int, bootmethod::BootMethod, blmethod::BLMethod,
                       flevel1::Function, flevel2::Function, numobsperresample::Int)
		blocklength <= 0.0 && error("Block length must be strictly positive")
		new(numobs, blocklength, numresample, bootmethod, blmethod, flevel1, flevel2, numobsperresample)
	end
end

#BootInput constructors
BootInput() = BootInput(0, 1.0, 0, NUM_RESAMPLE, BootStationary(), BLPPW2009(), mean, var, 0)
function BootInput(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                    blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))::BootInput
	blmethod == :dummy ? (blmethod = bootmethodtoblmethod(bootmethod)) : (blmethod = symboltoblmethod(blmethod)) #auto-adjust blmethod if not specified
	blocklength <= 0.0 && (blocklength = optblocklength(data, blmethod))
    return BootInput(data_length(data), Float64(blocklength), Int(numresample), symboltobootmethod(bootmethod), blmethod, flevel1, flevel2, Int(numobsperresample))
end

#Local function for getting the number of observations in a dataset
(data_length(x::AbstractVector{T})::Int) where {T} = length(x)
(data_length(x::AbstractMatrix{T})::Int) where {T} = size(x, 1)
function data_length(x::Vector{Vector{T}})::Int where {T}
	length(x) == 0 && error("Empty input dataset")
    any(length(x[1]) .!= [ length(y) for y in x ]) && error("Inner vector length mismatch in input data")
	return length(x[1])
end
