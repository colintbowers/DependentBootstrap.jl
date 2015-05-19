module DependentBootstrap
#-----------------------------------------------------------
#PURPOSE
#	Colin T. Bowers module for dependent bootstraps
#NOTES
#	The following functions will probably need new methods if you add a new bootstrap method
#		dBootstrapIndex (probably)
#		replaceBlockLength!
#	Throughout this module, automatic block length estimation (and hence the type BlockLengthMethod) is only ever used if the current block length is <= 0 (note, default block lengths are -1)
#LICENSE
#	MIT License (see github repository for more detail: https://github.com/colintbowers/DependentBootstrap.jl.git)
#-----------------------------------------------------------


#Load any entire modules that are needed (use import ModuleName1, ModuleName2, etc)
using 	Distributions,
		StatsBase,
		KernelStat

#Load any specific variables/functions that are needed (use import ModuleName1.FunctionName1, ModuleName2.FunctionName2, etc)
import 	Base.string,
		Base.show,
		Base.copy

#Specify the variables/functions to export (use export FunctionName1, FunctionName2, etc)
export 	BootstrapMethod,
		BlockLengthMethod,
		BootstrapIID,
		BootstrapStationary,
		BootstrapMovingBlock,
		BootstrapCircularBlock,
		BootstrapNonoverlappingBlock,
		BootstrapTaperedBlock,
		BlockLengthPPW2009,
		BlockLengthPP2002,
		BootstrapParam,
		getBlockLength,
		replaceNumObsData!,
		replaceNumObsResample!,
		replaceNumResample!,
		replaceBlockLength!,
		dBootstrapBlockLength,
		dBootstrapBlockLength!,
		dBootstrapIndex,
		dBootstrapIndex!,
		dBootstrapData,
		dBootstrapData!,
		dBootstrapStatistic,
		dBootstrapStatistic!,
		dBootstrap,
		dBootstrap!,
		dBootstrapVar,
		dBootstrapStd,
		dBootstrapConf

#******************************************************************************

#----------------------------------------------------------
#SET CONSTANTS FOR MODULE
#----------------------------------------------------------
const defaultNumResample = 1000::Int #Sets the default number of re-samples used in BootstrapParam constructors
const defaultBlockLength = -1::Int #A default block length that is useful for constructors that initialize a BootstrapParam without estimating the block length
const methodsThatUseWeighting = ["taperedBlock"]::Vector{ASCIIString} #Vector of bootstrap methods that necessitate weighting (ie a call to dBootstrapWeight!)
const validBootstrapMethodsString = ["iid","stationary", "nonoverlappingBlock", "movingBlock", "circularBlock", "taperedBlock"]::Vector{ASCIIString} #Vector of valid string representations for bootstrap methods
const validBlockLengthMethodString = ["PPW2009", "PP2002"]::Vector{ASCIIString} #Vector of valid string representations for block length detection methods





#----------------------------------------------------------
#TYPE (ABSTRACT)
#	BootstrapMethod
#	BlockLengthMethod
#PURPOSE
#	An abstract type to nest all bootstrap method types or block length selection method types
#NOTES
#----------------------------------------------------------
abstract BootstrapMethod
abstract BlockLengthMethod





#----------------------------------------------------------
#TYPES
#   IIDBoot <: BootstrapMethod
#	StationaryBoot <: BootstrapMethod
#	MovingBlockBoot <: BootstrapMethod
#	NonoverlappingBlockBoot <: BootstrapMethod
#	CircularBlockBoot <: BootstrapMethod
#	TaperedBlockBoot <: BootstrapMethod
#FIELDS
#	Fields are unique to each type and hold the parameters necessary to implement that particular bootstrap
#PURPOSE
#	This type stores the information needed to select and run a particular statistical bootstrap method
#CONSTRUCTORS
#	Outer constructors are unique to the type and are generally used to provide popular default values.
#METHODS
#	copy(b::BootstrapMethod)
#	string(b::BootstrapMethod)
#	show(b::BootstrapMethod)
#	toBootstrapMethod(x::String): Converts x to corresponding bootstrap method
#NOTES
#	Block length values can be non-positive. This is understood by the functions in this module to imply that the block-length needs to be detected.
#----------------------------------------------------------
#---- TYPE DEFINTIIONS ----
#iid Bootstrap
type BootstrapIID <: BootstrapMethod #This bootstrap method does not have any associated parameters
end
#Stationary bootstrap
type BootstrapStationary <: BootstrapMethod
	expectedBlockLength::Float64
end
BootstrapStationary() = BootstrapStationary(-1.0)
BootstrapStationary(expectedBlockLength::Int) = BootstrapStationary(convert(Float64, expectedBlockLength))
#Moving block bootstrap
type BootstrapMovingBlock <: BootstrapMethod
	blockLength::Int
end
BootstrapMovingBlock() = BootstrapMovingBlock(-1)
BootstrapMovingBlock(blockLength::Float64) = BootstrapMovingBlock(convert(Int, ceil(blockLength)))
#Nonoverlapping block bootstrap
type BootstrapNonoverlappingBlock <: BootstrapMethod
	blockLength::Int
end
BootstrapNonoverlappingBlock() = BootstrapNonoverlappingBlock(-1)
BootstrapNonoverlappingBlock(blockLength::Float64) = BootstrapNonoverlappingBlock(convert(Int, ceil(blockLength)))
#Circular block bootstrap
type BootstrapCircularBlock <: BootstrapMethod
	blockLength::Int
end
BootstrapCircularBlock() = BootstrapCircularBlock(-1)
BootstrapCircularBlock(blockLength::Float64) = BootstrapCircularBlock(convert(Int, ceil(blockLength)))
#Tapered block bootstrap
type BootstrapTaperedBlock <: BootstrapMethod
	blockLength::Int
	kernelFunction::KernelFunction
	function BootstrapTaperedBlock(blockLength::Int, kernelFunction::KernelFunction)
		if typeof(kernelFunction) != KernelPP2002Trap && typeof(kernelFunction) != KernelPP2002Smooth
			error("Invalid kernel function for tapered block bootstrap method")
		end
		new(blockLength, kernelFunction)
	end
end
BootstrapTaperedBlock() = BootstrapTaperedBlock(-1, KernelPP2002Trap())
BootstrapTaperedBlock(blockLength::Int) = BootstrapTaperedBlock(blockLength, KernelPP2002Trap())
BootstrapTaperedBlock(blockLength::Float64) = BootstrapTaperedBlock(convert(Int, ceil(blockLength)), KernelPP2002Trap())
#---- METHODS ----
#copy methods
copy(b::BootstrapIID) = BootstrapIID()
copy(b::BootstrapStationary) = BootstrapStationary(copy(b.expectedBlockLength))
copy(b::BootstrapMovingBlock) = BootstrapMovingBlock(copy(b.blockLength))
copy(b::BootstrapCircularBlock) = BootstrapCircularBlock(copy(b.blockLength))
copy(b::BootstrapNonoverlappingBlock) = BootstrapNonoverlappingBlock(copy(b.blockLength))
copy(b::BootstrapTaperedBlock) = BootstrapTaperedBlock(copy(b.blockLength), copy(b.kernelFunction))
#string
string(x::BootstrapIID) = "iid"
string(x::BootstrapStationary) = "stationary"
string(x::BootstrapNonoverlappingBlock) = "nonoverlappingBlock"
string(x::BootstrapMovingBlock) = "movingBlock"
string(x::BootstrapCircularBlock) = "circularBlock"
string(x::BootstrapTaperedBlock) = "taperedBlock"
#show method
show(io::IO, b::BootstrapIID) = println(io, "bootstrap method = " * string(b))
function show(io::IO, b::BootstrapStationary)
	println(io, "bootstrap method = " * string(b))
	println(io, "    expected block length = " * b.expectedBlockLength)
end
function show{T<:Union(BootstrapMovingBlock, BootstrapCircularBlock, BootstrapNonoverlappingBlock)}(io::IO, b::T)
	println(io, "bootstrap method = " * string(b))
	println(io, "    block length = " * b.blockLength)
end
function show(io::IO, b::BootstrapTaperedBlock)
	println(io, "bootstrap method = " * string(b))
	println(io, "    block length = " * b.blockLength)
	println(io, "    kernel function = " * string(b.kernelFunction))
end
#show wrapper on STDOUT
show(b::BootstrapMethod) = show(STDOUT, b)
#Retrieve the block length from a BootstrapMethod
function getBlockLength(bm::BootstrapMethod)
	if typeof(bm) == BootstrapIID
		return(1.0)
	elseif typeof(bm) == BootstrapStationary
		return(bm.expectedBlockLength)
	else
		return(convert(Float64, bm.blockLength))
	end
end
#toBootstrapMethod: This function is used to convert string representations of bootstrap methods into <: BootstrapMethod
function toBootstrapMethod(x::ASCIIString) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "iid"
		return(BootstrapIID())
	elseif x == "stationary"
		return(BootstrapStationary())
	elseif x == "nonoverlappingBlock"
		return(BootstrapNonoverlappingBlock())
	elseif x == "movingBlock"
		return(BootstrapMovingBlock())
	elseif x == "circularBlock"
		return(BootstrapCircularBlock())
	elseif x == "taperedBlock"
		return(BootstrapTaperedBlock())
	else
		if any(x .== validBootstrapMethodsString)
			error("You must also supply a block length for the input bootstrap method")
		else
			error("Invalid string representation bootstrap method ")
		end
	end
end
function toBootstrapMethod(x::ASCIIString, blockLength::Number) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "iid"
		return(BootstrapIID())
	elseif x == "stationary"
		return(BootstrapStationary(convert(Float64, ceil(blockLength))))
	elseif x == "nonoverlappingBlock"
		return(BootstrapNonoverlappingBlock(convert(Int, ceil(blockLength))))
	elseif x == "movingBlock"
		return(BootstrapMovingBlock(convert(Int, ceil(blockLength))))
	elseif x == "circularBlock"
		return(BootstrapCircularBlock(convert(Int, ceil(blockLength))))
	elseif x == "taperedBlock"
		return(BootstrapTaperedBlock(convert(Int, ceil(blockLength)), KernelPP2002Trap()))
	else
		error("Invalid string representation of bootstrap method")
	end
end
function toBootstrapMethod(x::ASCIIString, blockLength::Number, p1::Any)
	if x == "taperedBlock"
		if typeof(p1) <: KernelFunction
			if !(typeof(p1) == KernelPP2002Trap || typeof(p1) == KernelPP2002Smooth)
				error("Invalid choice of kernel function for use with tapered block bootstrap")
			else
				return(BootstrapTaperedBlock(convert(Int, ceil(blockLength)), kernelFunc))
			end
		elseif typeof(p1) == ASCIIString
			if !(p1 == "PP2002Trap" || p1 == "PP2002Smooth")
				error("Invalid choice of kernel function for use with tapered block bootstrap")
			else
				return(BootstrapTaperedBlock(convert(Int, ceil(blockLength)), toKernelFunction(kernelFunc)))
			end
		else
			error("Invalid type for third input given method is tapered block bootstrap")
		end
	else
		if any(x .== validBootstrapMethodsString)
			error("A third input is not necessary for the input bootstrap method string")
		else
			error("Invalid string representation of bootstrap method")
		end
	end
end



#----------------------------------------------------------
#TYPES
#   BlockLengthPPW2009 <: BlockLengthMethod
#	BlockLengthPP2002 <: BlockLengthMethod
#FIELDS
#	Fields are unique to each type and hold the parameters necessary to implement that particular block length selection method
#PURPOSE
#	This type stores the information needed to select and run a particular block length selection method
#CONSTRUCTORS
#	Outer constructors are unique to the type and are generally used to provide popular default values.
#NOTES
#----------------------------------------------------------
#----------- TYPE DEFINITIONS -----------
#Block length detection procedure of Patton, Politis, White (2009)
type BlockLengthPPW2009 <: BlockLengthMethod
	bandwidthMethod::BandwidthMethod
	bootstrapMethod::ASCIIString #Use ASCIIString here instead of BootstrapMethod because we don't care about block length (it is only needed because there are some minor differences in method depending on whether dealing with stationary boostrap versus other block bootstrap methods)
	function BlockLengthPPW2009{T<:BandwidthMethod}(bandwidthMethod::T, bootstrapMethod::ASCIIString)
		if !(bootstrapMethod == "stationary" || bootstrapMethod == "circularBlock" || bootstrapMethod == "movingBlock")
			error("BlockLengthPPW2009 detection method only viable for stationary, circular block, and moving block bootstraps")
		end
		new(bandwidthMethod, bootstrapMethod)
	end
end
BlockLengthPPW2009() = BlockLengthPPW2009(BandwidthP2003(), "stationary")
BlockLengthPPW2009(bootstrapMethod::ASCIIString) = BlockLengthPPW2009(BandwidthP2003(), bootstrapMethod)
BlockLengthPPW2009(bootstrapMethod::ASCIIString, numObs::Int) = BlockLengthPPW2009(BandwidthP2003(numObs), bootstrapMethod)
BlockLengthPPW2009(bandwidthMethod::BandwidthMethod) = BlockLengthPPW2009(bandwidthMethod, "stationary")
BlockLengthPPW2009(bandwidthMethod::ASCIIString, bootstrapMethod::ASCIIString) = BlockLengthPPW2009(KernelStat.toBandwidthMethod(bandwidthMethod), "stationary")
#Block length detection procedure of Paparoditis, Politis (2002) [NOTE: ONLY FOR TAPERED BLOCK BOOTSTRAP]
type BlockLengthPP2002 <: BlockLengthMethod
	bandwidthMethod::BandwidthMethod
	kernelFunction::KernelFunction
	function BlockLengthPP2002(bandwidthMethod::BandwidthMethod, kernelFunction::KernelFunction)
		if !(typeof(kernelFunction) == KernelPP2002Trap || typeof(kernelFunction) == KernelPP2002Smooth)
			error("Tapered block bootstrap method can only be implemented with KernelPP2002Trap or KernelPP2002Smooth kernel functions")
		end
		new(bandwidthMethod, kernelFunction)
	end
end
BlockLengthPP2002() = BlockLengthPP2002(BandwidthP2003(), KernelPP2002Trap())
BlockLengthPP2002(numObs::Int) = BlockLengthPP2002(BandwidthP2003(numObs), KernelPP2002Trap())
BlockLengthPP2002(kernelFunction::KernelFunction) = BlockLengthPP2002(BandwidthP2003(), kernelFunction)
BlockLengthPP2002(kernelFunction::KernelFunction, numObs::Int) = BlockLengthPP2002(BandwidthP2003(numObs), kernelFunction)
BlockLengthPP2002(bandwidthMethod::BandwidthMethod) = BlockLengthPP2002(bandwidthMethod, KernelPP2002Trap())
BlockLengthPP2002(bandwidthMethod::ASCIIString) = BlockLengthPP2002(KernelStat.toBandwidthMethod(bandwidthMethod), KernelPP2002Trap())
BlockLengthPP2002(bandwidthMethod::ASCIIString, kernelFunction::ASCIIString) = BlockLengthPP2002(KernelStat.toBandwidthMethod(bandwidthMethod), toKernelFunction(kernelFunction))
#------------ METHODS -------------
#copy function
copy(b::BlockLengthPPW2009) = BlockLengthPPW2009(copy(b.bandwidthMethod), copy(b.bootstrapMethod))
copy(b::BlockLengthPP2002) = BlockLengthPP2002(copy(b.bandwidthMethod), copy(b.kernelFunction))
#string function for string representation of block length method
string(b::BlockLengthPPW2009) = "PPW2009"
string(b::BlockLengthPP2002) = "PP2002"
#show method for block length methods
function show(io::IO, b::BlockLengthPPW2009)
	println(io, "block length method = " * string(b))
	println(io, "    bandwidth method = " * string(b.bandwidthMethod))
	println(io, "    bootstrap method = " * b.bootstrapMethod)
end
function show(io::IO, b::BlockLengthPP2002)
	println(io, "block length method = " * string(b))
	println(io, "    bandwidth method = " * string(b.bandwidthMethod))
	println(io, "    kernel function = " * string(b.kernelFunction))
end
#show wrapper for STDOUT
show(b::BlockLengthMethod) = show(STDOUT, b)
#toBlockLengthMethod: This function is used to convert string representations of block length methods into <: BlockLengthMethod
function toBlockLengthMethod(x::ASCIIString) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "PPW2009"
		return(BlockLengthPPW2009())
	elseif x == "PP2002"
		return(BlockLengthPP2002())
	else
		error("Invalid string representation of block length method")
	end
end
function toBlockLengthMethod(x::ASCIIString, bandwidthMethod::BandwidthMethod) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "PPW2009"
		return(BlockLengthPPW2009(bandwidthMethod))
	elseif x == "PP2002"
		return(BlockLengthPP2002(bandwidthMethod))
	else
		error("Invalid string representation of block length method")
	end
end
function toBlockLengthMethod(x::ASCIIString, bandwidthMethod::ASCIIString) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "PPW2009"
		return(BlockLengthPPW2009(KernelStat.toBandwidthMethod(bandwidthMethod)))
	elseif x == "PP2002"
		return(BlockLengthPP2002(KernelStat.toBandwidthMethod(bandwidthMethod)))
	else
		error("Invalid string representation of block length method")
	end
end
function toBlockLengthMethod(x::ASCIIString, bandwidthMethod::BandwidthMethod, p1::ASCIIString) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "PPW2009"
		return(BlockLengthPPW2009(bandwidthMethod, p1))
	elseif x == "PP2002"
		return(BlockLengthPP2002(bandwidthMethod, p1))
	else
		error("Invalid string representation of block length method")
	end
end
function toBlockLengthMethod(x::ASCIIString, bandwidthMethod::ASCIIString, p1::ASCIIString) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "PPW2009"
		return(BlockLengthPPW2009(bandwidthMethod, p1))
	elseif x == "PP2002"
		return(BlockLengthPP2002(bandwidthMethod, p1))
	else
		error("Invalid string representation of block length method")
	end
end








#----------------------------------------------------------
#TYPE
#   BootstrapParam
#FIELDS
#	Fields are:
#		numObsData::Int: Number of observations in the dataset that is to be bootstrapped
#		numObsResample::Int #Number of observations to be drawn (with replacement) per resample. For the vast majority of cases this field will equal numObsData.
#		numResample::Int #Number of resamples
#		bootstrapMethod::BootstrapMethod: Bootstrap method to use (see bootstrap method type definitions for more detail)
#		blockLengthMethod::BlockLengthMethod: Method to use to detect block length (only applies if current block length is non-positive) (see block length method types for more detail)
#		statistic::Union(Function, ASCIIString): The statistic of interest. If a function, then must be able to be applied to vector of observations and return a single value. If a string, must be one of the following recognized functions:
#			"mean": The sample mean
#			"median": The sample median
#			"variance": The sample variance
#			"var": Identical to variance
#			"std": The sample standard deviation
#			"sum": A sum of the observations (can be thought of as an unscaled mean - useful in cases where the variance of the data naturally shrinks as the number of observations grows, e.g. infill asymptotics)
#			"quantileXX": The sample quantile where "XX" is a two digit number corresponding to the quantile percentage, e.g. "quantile95" yields the 95% quantile
#		distributionParam::Union(Function, ASCIIString): The distribution parameter of the distribution of the statistic of interest. If a function, must be able to be applied to a vector of observations and return a single value.
#			"mean": The expected value of the distribution (This option is irrelevant in many cases)
#			"median": The median of the distirbution (This option is irrelevant in many cases)
#			"variance": The variance of the distribution (A popular option)
#			"var": Identical to variance
#			"std": The stadard deviation of the distribution
#			"quantileXX": A quantile of the distribution
#			"conf": A 95% confidence interval from the distribution (If you want an interval other than 95% you'll need to input it as a function)
#PURPOSE
#	This type stores the information needed to select and run a particular statistical bootstrap method
#CONSTRUCTORS
#	Many different constructors (see below). Keyword ones are particuarly important as all functions that accept keywords essentially are just employing the keyword constructor on a BootstrapParam
#METHODS
#	copy(b::BootstrapParam)
#	show(b::BootstrapParam)
#	getBlockLength(b::BootstrapParam): Returns the block length (or expected block length) associated with a BootstrapParam
#	replaceBlockLength!(b::BootstrapParam, newBlockLength<:Number): Replace the block length in a BootstrapParam with a new one
#	replaceNumObsData!(b::BootstrapParam, newObs::Int): Replace the number of observations of data field in a BootstrapParam.
#	replaceNumObsResample!(b::BootstrapParam, newObs::Int): Replace the number of observations per resample field in a BootstrapParam.
#	replaceNumResample!(b::BootstrapParam, newObs::Int): Replace the number of resamples field in a BootstrapParam.
#NOTES
#	In the vast majority of cases, numObsData = numObsResample. The constructors reflect this.
#----------------------------------------------------------
#----------- TYPE DEFINITION --------------------
type BootstrapParam
	numObsData::Int #Number of observations in the dataset that is to be bootstrapped
	numObsResample::Int #Number of observations to be drawn (with replacement) per resample
	numResample::Int #Number of resamples
	bootstrapMethod::BootstrapMethod #Bootstrap method to use
	blockLengthMethod::BlockLengthMethod #Method to use to detect block length (only applies if current block length is non-positive)
	statistic::Union(Function, ASCIIString) #Statistic of interest
	distributionParam::Union(Function, ASCIIString) #Parameter of distribution of statistic of interest
	function BootstrapParam(numObsData::Int, numObsResample::Int, numResample::Int, bootstrapMethod::BootstrapMethod, blockLengthMethod::BlockLengthMethod, statistic::Union(Function, ASCIIString), distributionParam::Union(Function, ASCIIString))
		if numObsData < 1
			error("Number of observations in dataset must be greater than zero to apply a bootstrap")
		end
		if numObsResample < 1
			error("Number of observations per resample must be greater than zero")
		end
		if numResample < 1
			error("Number of resamples to bootstrap must be greater than zero")
		end
		new(numObsData, numObsResample, numResample, bootstrapMethod, blockLengthMethod, statistic, distributionParam)
	end
end
#------------ OUTER CONSTRUCTORS ----------------
#Constructors with numObsData as first input set block length equal to defaultBlockLength (-1)
BootstrapParam(numObsData::Int) = BootstrapParam(numObsData, numObsData, defaultNumResample, BootstrapStationary(defaultBlockLength), BlockLengthPPW2009(BandwidthP2003(), "stationary"), "mean", "variance")
BootstrapParam(numObsData::Int, statistic::Union(Function, ASCIIString), distributionParam::Union(Function, ASCIIString)) = BootstrapParam(numObsData, numObsData, defaultNumResample, BootstrapStationary(defaultBlockLength), BlockLengthPPW2009(BandwidthP2003(), "stationary"), statistic, distributionParam)
BootstrapParam(numObsData::Int, numResample::Int) = BootstrapParam(numObsData, numObsData, numResample, BootstrapStationary(defaultBlockLength), BlockLengthPPW2009(BandwidthP2003(), "stationary"), "mean", "variance")
BootstrapParam(numObsData::Int, numResample::Int, statistic::Union(Function, ASCIIString), distributionParam::Union(Function, ASCIIString)) = BootstrapParam(numObsData, numObsData, numResample, BootstrapStationary(defaultBlockLength), BlockLengthPPW2009(BandwidthP2003(), "stationary"), statistic, distributionParam)
#Constructors with x (ie underlying data) as first input will attempt to estimate the appropriate block length from that data
BootstrapParam{T<:Number}(x::Vector{T}) = BootstrapParam(length(x), length(x), defaultNumResample, BootstrapStationary(dBootstrapBlockLength(x, BlockLengthPPW2009(BandwidthP2003(), "stationary"))), BlockLengthPPW2009(BandwidthP2003(), "stationary"), "mean", "variance")
BootstrapParam{T<:Number}(x::Vector{T}, statistic::Union(Function, ASCIIString), distributionParam::Union(Function, ASCIIString)) = BootstrapParam(length(x), length(x), defaultNumResample, BootstrapStationary(dBootstrapBlockLength(x, BlockLengthPPW2009(BandwidthP2003(), "stationary"))), BlockLengthPPW2009(BandwidthP2003(), "stationary"), statistic, distributionParam)
BootstrapParam{T<:Number}(x::Vector{T}, numResample::Int) = BootstrapParam(length(x), length(x), numResample, BootstrapStationary(dBootstrapBlockLength(x, BlockLengthPPW2009(BandwidthP2003(), "stationary"))), BlockLengthPPW2009(BandwidthP2003(), "stationary"), "mean", "variance")
BootstrapParam{T<:Number}(x::Vector{T}, numResample::Int, statistic::Union(Function, ASCIIString), distributionParam::Union(Function, ASCIIString)) = BootstrapParam(length(x), length(x), numResample, BootstrapStationary(dBootstrapBlockLength(x, BlockLengthPPW2009(BandwidthP2003(), "stationary"))), BlockLengthPPW2009(BandwidthP2003(), "stationary"), statistic, distributionParam)
#Constructors that employ keyword arguments follow
function BootstrapParam(numObsData::Int; numObsResample::Int=numObsData, numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean", distributionParam::ASCIIString="variance")
	if blockLengthMethod == "PPW2009"
		blm = toBlockLengthMethod(blockLengthMethod, bandwidthMethod, bootstrapMethod)
	else
		blm = toBlockLengthMethod(blockLengthMethod, bandwidthMethod)
	end
	bm = toBootstrapMethod(bootstrapMethod, blockLength)
	return(BootstrapParam(numObsData, numObsResample, numResample, bm, blm, statistic, distributionParam))
end
function BootstrapParam{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean", distributionParam::ASCIIString="variance")
	numObsData = length(x)
	if blockLengthMethod == "PPW2009"
		blm = toBlockLengthMethod(blockLengthMethod, bandwidthMethod, bootstrapMethod)
	else
		blm = toBlockLengthMethod(blockLengthMethod, bandwidthMethod)
	end
	if blockLength <= 0 #Update block length using x and blm (if necessary)
		if bootstrapMethod != "iid"
			blockLength = dBootstrapBlockLength(x, blm)
		end
	end
	bm = toBootstrapMethod(bootstrapMethod, blockLength)
	return(BootstrapParam(numObsData, numObsResample, numResample, bm, blm, statistic, distributionParam))
end
#------------ METHODS ------------------
#copy
copy(p::BootstrapParam) = BootstrapParam(copy(p.numObsData), copy(p.numObsResample), copy(p.numResample), copy(p.bootstrapMethod), copy(p.blockLengthMethod), copy(p.statistic), copy(p.distributionParam))
#show methods
function show(io::IO, b::BootstrapParam)
	println(io, "dependent bootstrap parameter values:")
	println(io, "    number of observations in dataset = " * string(b.numObsData))
	println(io, "    number of observations per resample = " * string(b.numObsResample))
	println(io, "    number of resamples = " * string(b.numResample))
	println(io, "    bootstrap method = " * string(b.bootstrapMethod))
	println(io, "    block length method = " * string(b.blockLengthMethod))
	println(io, "    current block length = " * string(getBlockLength(b)))
	println(io, "    statistic of interest = " * string(b.statistic))
	println(io, "    distribution parameter of statistic = " * string(b.distributionParam))
end
#show wrapper on STDOUT
show(b::BootstrapParam) = show(STDOUT, b)
#Retrieve the block length
getBlockLength(bp::BootstrapParam) = getBlockLength(bp.bootstrapMethod)
#Replace particular fields and incorporate an error trap
function replaceNumObsData!(param::BootstrapParam, newNumObsData::Int)
	if newNumObsData < 1
		error("Number of observations in dataset must be greater than zero to apply a bootstrap")
	end
	param.numObsData = newNumObsData
	return(param)
end
function replaceNumObsResample!(param::BootstrapParam, newNumObsResample::Int)
	if newNumObsResample < 1
		error("Number of observations per resample must be greater than zero")
	end
	param.numObsResample = newNumObsResample
	return(param)
end
function replaceNumResample!(param::BootstrapParam, newNumResample::Int)
	if newNumResample < 1
		error("Number of resamples to bootstrap must be greater than zero")
	end
	param.numResample = newNumResample
	return(param)
end
replaceBlockLength!{T<:Number}(bp::BootstrapParam, newBlockLength::T) = replaceBlockLength!(bp.bootstrapMethod, newBlockLength)
replaceBlockLength!{T<:Number}(bm::BootstrapIID, newBlockLength::T) = error("iid bootstrap does not require a block length")
replaceBlockLength!{T<:Number}(bm::BootstrapStationary, newBlockLength::T) = (bm.expectedBlockLength = convert(Float64, newBlockLength))
replaceBlockLength!{TBoot<:Union(BootstrapMovingBlock, BootstrapNonoverlappingBlock, BootstrapCircularBlock, BootstrapTaperedBlock), T<:Number}(bm::TBoot, newBlockLength::T) = (bm.blockLength = convert(Int, ceil(newBlockLength)))





#----------------------------------------------------------
#FUNCTION
#	dBootstrapBlockLength
#	dBootstrapBlockLength!
#INPUT
#	Let T<:Any and let x::Vector{T} denote the vector of original observed data. Function accepts:
#		(x::Vector{T}, bm::BlockLengthPPW2009): Perform block length selection procedure of Politis, White (2004) and Patton, Politis, and White (2009), where x is the vector of observed data. Note, computation differs for stationary versus cicular block and moving block bootstraps. Fields of BlockLengthPPW2009 reflect this.
#		(x::Vector{T}, bm::BlockLengthPP2002): Perform block length selection procedure of Paproditis, Politis (2002), where x is the vector of observed data. Note, computation is only valid for the tapered block bootstrap.
#		(x::Vector{T}, bp::BootstrapParam): Wrapper on a BootstrapParam
#		(bp::BootstrapParam, x::Vector{T}): Wrapper so order doesn't matter
#KEYWORD INPUT METHOD
#		dBootstrapBlockLength{T<:Number}(x::Vector{T}; blockLengthMethod::ASCIIString="PPW2009", bootstrapMethod::ASCIIString="stationary", bandwidthMethod::ASCIIString="P2003")
#IN-PLACE UPDATE METHODS, ie for dBootstrapBlockLength!
#		(x::Vector{T}, bp::BootstrapParam): Updates bp with the new block length in-place
#		(bp::BootstrapParam, x::Vector{T}): Wrapper so order doesn't matter
#OUTPUT
#	Estimated block length for dependent bootstrap procedure, expressed as Float64.
#PURPOSE
#	The purpose of this function is to estimate the block length to use with a given dependent bootstrap procedure.
#NOTES
#	This function uses the non-exported function dBootstrapBlockLength_MAndCorVec that performs some routines that are common to several block length methods
#	This function uses the non-exported function dBootstrapBlockLength_KernelCorVec that performs some routines that are common to several block length methods
#----------------------------------------------------------
#Block length selection method of Patton, Politis, White (2009) "Correction to Automatic Block Length Selection For the Dependent Bootstrap"
function dBootstrapBlockLength{T<:Number}(x::Vector{T}, bm::BlockLengthPPW2009)
	if length(x) < 3
		error("You must have at least 3 observations to estimate block length")
	end
	(M, xVar, covVec) = dBootstrapBlockLength_MAndCorVec(x, bm)
	kernelCovVec = dBootstrapBlockLength_KernelCovVec(covVec, KernelP2003FlatTop(), M)
	gHat = 0.0
	for k = 1:M
		gHat += 2 * k * kernelCovVec[k] #note, "2*" is because sum is from -M and M, but is symmetric about 0. Note, 0 term in sum vanishes since k=0 -> |k|=0
	end
	if bm.bootstrapMethod == "stationary"
		dHat = 2 * (xVar + 2*sum(kernelCovVec))^2 #note, in expression (1 + 2*sum(kernelCovVec)), "2*" is because sum is from -M to M, but is symmetric around 0. "1+" is the zero term of the sum which is simply equal to unity.
	elseif bm.bootstrapMethod == "circularBlock" || bm.bootstrapMethod == "movingBlock"
		dHat = (4/3) * (xVar + 2*sum(kernelCovVec))^2 #note, in expression (1 + 2*sum(kernelCovVec)), "2*" is because sum is from -M to M, but is symmetric around 0. "1+" is the zero term of the sum which is simply equal to unity.
	else
		error("Optimal block length given BlockLengthPPW2009 method is only known for stationary, circular, or moving block bootstraps")
	end
	blockLength = (2 * gHat^2 * length(x) / dHat)^(1/3)
	blockLength = min(blockLength, ceil(min(3*sqrt(length(x)), length(x) / 3))) #Enforce upper bound on block length suggested by Patton
	return(blockLength) #Equation 9 and 14 from Politis, White (2004)
end
#Block length selection method of Paparoditis, Politis (2002) "The Tapered Block Bootstrap for General Statistics From Stationary Sequences"
function dBootstrapBlockLength{T<:Number}(x::Vector{T}, bm::BlockLengthPP2002)
	if length(x) < 3
		error("You must have at least 3 observations to estimate block length")
	end
	(M, xVar, covVec) = dBootstrapBlockLength_MAndCorVec(x, bm)
	kernelCovVec = dBootstrapBlockLength_KernelCovVec(covVec, KernelP2003FlatTop(), M)
	deltaUnknown = (xVar + 2*sum(kernelCovVec))^2 #Unknown parameter in Delta expression, section 4, Paparoditis, Politis (2002). (note, "1+" is for k=0 term in summation, "2*" is because summation is from -M to M, but is symmetric about 0)
	gammaUnknown = 0.0 #Unknown parameter in Gamma expression, start of section 4
	for k = 1:M
		gammaUnknown += 2 * k^2 * kernelCovVec[k] #note, "2*" is because sum is from -M and M, but is symmetric about 0. Note, 0 term in sum vanishes since k=0 -> k^2=0
	end
	if typeof(bm.kernelFunction) == KernelPP2002Trap
		gammaHat = -5.45 * gammaUnknown #-5.45 = (1/2) * -10.9 = (1/2) * (w*w)''(0) / (w*w)(0) [OPTIMAL TRAP VALUES], see Paparoditis, Politis (2002)
		deltaHat = 1.099 * deltaUnknown #1.099 = 2 * 0.5495 = 2 * int_{-1}^{1} ((w*w)^2(x) / (w*w)^2(0)) dx [OPTIMAL TRAP VALUES], see Paparoditis, Politis (2002)
	elseif typeof(bm.kernelFunction) == KernelPP2002Smooth
		gammaHat = -5.175 * gammaUnknown #-5.175 = (1/2) * -10.35 = (1/2) * (w*w)''(0) / (w*w)(0) [OPTIMAL SMOOTH VALUES], see Paparoditis, Politis (2002)
		deltaHat = 1.1312 * deltaUnknown #1.1312 = 2 * 0.5656 = 2 * int_{-1}^{1} ((w*w)^2(x) / (w*w)^2(0)) dx [OPTIMAL SMOOTH VALUES], see Paparoditis, Politis (2002)
	else
		error("Optimal parameters for input kernel function for BlockLengthPP2002 selection procedure are unknkown. Please use KernelPP2002Trap or KernelPP2002Smooth.")
	end
	blockLength = (4 * gammaHat^2 * length(x) / deltaHat)^(1/5)
	blockLength = min(blockLength, ceil(min(3*sqrt(length(x)), length(x) / 3))) #Enforce upper bound on block length suggested by Patton
	return(blockLength) #Equation 20 from Paraproditis, Politis (2002)
end
#Wrapper methods
dBootstrapBlockLength{T<:Number}(x::Vector{T}, bp::BootstrapParam) = dBootstrapBlockLength(x, bp.blockLengthMethod)
dBootstrapBlockLength{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrapBlockLength(x, bp)
#Keyword wrapper methods
function dBootstrapBlockLength{T<:Number}(x::Vector{T}; blockLengthMethod::ASCIIString="PPW2009", bootstrapMethod::ASCIIString="stationary", bandwidthMethod::ASCIIString="P2003")
	if bootstrapMethod == "iid"
		error("It makes no sense to call this function when bootstrapMethod is set to iid")
	end
	bp = BootstrapParam(length(x), bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, bandwidthMethod=bandwidthMethod) #Use BootstrapParam constructor to deal with keyword arguments
	return(dBootstrapBlockLength(x, bp))
end
#In-place update of BootstrapParam version of function
dBootstrapBlockLength!{T<:Number}(x::Vector{T}, bp::BootstrapParam) = replaceBlockLength!(bp, dBootstrapBlockLength(x, bp.blockLengthMethod))
dBootstrapBlockLength!{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrapBlockLength!(x, bp)
#Non-exported functions used exclusively by dBootstrapBlockLength
function dBootstrapBlockLength_MAndCorVec{T<:Number}(x::Vector{T}, bm::BlockLengthMethod)
	if !(typeof(bm) == BlockLengthPPW2009 || typeof(bm) == BlockLengthPP2002)
		error("This non-exported function is not defined for the input block length method type")
	end
	(M, xVar, covVec) = bandwidth(x, bm.bandwidthMethod)
	if M > 0.5 * length(x)
		println("WARNING: Bandwidth in parameter estimation section of dBootstrapBlockLength forced to half total number of observations. Data may contain excessive dependence.")
		M = convert(Int, round(0.5 * length(x)))
	end
	if M < 2 #Even though M output of bandwidth function is always greater than 2, the above check for excessively large M can result in M being reduced below 2 (admittedly only in very unusual circumstances)
		M = 2
	end
	if length(covVec) < M #Get any additional autocovariances that we might need
		append!(covVec, autocov(x, length(covVec)+1:M))
	end
	return(M, xVar, covVec)
end
function dBootstrapBlockLength_KernelCovVec(covVec::Vector{Float64}, kF::KernelFunction, M::Int)
	if length(covVec) < M
		error("Logic fail in module. It should have been impossible for length(covVec) < M")
	end
	kernelCovVec = Array(Float64, M)
	for k = 1:M
		kernelCovVec[k] = evaluate(k / M, kF) * covVec[k]
	end
	return(kernelCovVec)
end









#----------------------------------------------------------
#FUNCTION
#	dBootstrapIndex
#	dBootstrapIndex!
#INPUT
#	Function accepts:
#		(bm::BootstrapIID, numObsData::Int, numObsResample::Int, numResample::Int): Perform the IID bootstrap
#		(bm::BootstrapStationary, numObsData::Int, numObsResample::Int, numResample::Int): Perform the stationary bootstrap
#		(bm::BootstrapMovingBlock, numObsData::Int, numObsResample::Int, numResample::Int): Perform the moving block bootstrap
#		(bm::BootstrapNonoverlappingBlock, numObsData::Int, numObsResample::Int, numResample::Int): Perform the non-overlapping block bootstrap
#		(bm::BootstrapCircularBlock, numObsData::Int, numObsResample::Int, numResample::Int): Perform the circular block bootstrap
#		(bm::BootstrapTaperedBlock, numObsData::Int, numObsResample::Int, numResample::Int): Gets indices for use with tapered-block bootstrap which is just the indices generated by the moving block bootstrap
#		(bm::BootstrapIID, numObsData::Int, numResample::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapIID, numObsData::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapStationary, numObsData::Int, numResample::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapStationary, numObsData::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapNonoverlappingBlock, numObsData::Int, numResample::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapNonoverlappingBlock, numObsData::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapMovingBlock, numObsData::Int, numResample::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapMovingBlock, numObsData::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapCircularBlock, numObsData::Int, numResample::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapCircularBlock, numObsData::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapTaperedBlock, numObsData::Int, numResample::Int): Wrapper that uses default values for some inputs
#		(bm::BootstrapTaperedBlock, numObsData::Int): Wrapper that uses default values for some inputs
#		(bp::BootstrapParam): Wrapper over a BootstrapParam
#		{T<:Number}(x::Vector{T}, bp::BootstrapParam): Wrapper that estimates the block length from data supplied, but only if current block length is non-positive
#KEYWORD METHOD (BLOCK LENGTH NOT SUPPLIED)
#		dBootstrapIndex{T<:Number}(numObsData::Int, blockLength::T; numObsResample::Int=numObsData, numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary")
#KEYWORD METHOD (BLOCK LENGTH ESTIMATED IF NECESSARY)
#		dBootstrapIndex{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003")
#METHOD FOR IN-PLACE UPDATE OF BLOCK LENGTH
#		dBootstrapIndex!{T<:Number}(x::Vector{T}, bp::BootstrapParam)
#		dBootstrapIndex!{T<:Number}(bp::BootstrapParam, x::Vector{T})
#OUTPUT
#	Matrix{Int} of indices that can be used to index into the data and construct re-sampled data.
#PURPOSE
#	The purpose of this function is to provide a set of bootstrap indices that can be used to index into a data vector to create re-sampled data
#NOTES
#	The statistic the user is bootstrapping is irrelevant to this function
#	Block length, or expected block length, should be determined (by data-driven methods or direct assignment) before calling this function, unless you use a method that will estimate block length
#----------------------------------------------------------
#iid bootstrap
function dBootstrapIndex(bm::BootstrapIID, numObsData::Int, numObsResample::Int, numResample::Int)
	if numObsData < 3
		error("You must have at least 3 observations to use dependent bootstrap routines")
	end
	return(rand(1:numObsData, numObsResample, numResample))
end
#Stationary bootstrap
function dBootstrapIndex(bm::BootstrapStationary, numObsData::Int, numObsResample::Int, numResample::Int)
	if numObsData < 3
		error("You must have at least 3 observations to use dependent bootstrap routines")
	end
	if bm.expectedBlockLength <= 0.0 || isnan(bm.expectedBlockLength)
		error("Expected block length is non-positive or NaN. Try a method that includes the underlying data (these methods typically auto-estimate the block length) or else supply the block length explicitly.")
	end
	if bm.expectedBlockLength <= 1.0
		inds = dBootstrapIndex(BootstrapIID(), numObsData, numObsResample, numResample)
	else
		inds = Array(Int, numObsResample, numResample) #Pre-allocate output
		for n in 1:numResample #Fill first row of inds with discrete uniform draws
			inds[1, n] = rand(1:numObsData)
		end
		geo1 = Geometric(1 / bm.expectedBlockLength)
		for n = 1:numResample #Fill remainder using geometric draws
			geoDraw = rand(geo1) + 1
			c = 1
			for s = 2:numObsResample
				if c == geoDraw #Start a new block
					inds[s, n] = rand(1:numObsData)
					geoDraw = rand(geo1) + 1
					c = 1
				else #Next obs in existing block
					if inds[s-1, n] == numObsData
						inds[s, n] = 1
					else
						inds[s, n] = inds[s-1, n] + 1
					end
					c += 1
				end
			end
		end
	end
	return(inds)
end
#Moving block bootstrap
function dBootstrapIndex(bm::BootstrapMovingBlock, numObsData::Int, numObsResample::Int, numResample::Int)
	if numObsData < 3
		error("You must have at least 3 observations to use dependent bootstrap routines")
	end
	if bm.blockLength <= 0 || isnan(bm.blockLength)
		error("Expected block length is non-positive or NaN. Try a method that includes the underlying data (these methods typically auto-estimate the block length) or else supply the block length explicitly.")
	end
	if bm.blockLength == 1
		inds = dBootstrapIndex(BootstrapIID(), numObsData, numObsResample, numResample)
	else
		inds = Array(Int, numObsResample, numResample) #Pre-allocate output
		blockStartUB = max(1, numObsData-bm.blockLength+1)
		for n = 1:numResample
			for s = 1:bm.blockLength:numObsResample
				inds[s, n] = rand(1:blockStartUB) #Start of block
				for sInner = s+1:min(s+bm.blockLength-1, numObsResample) #Iterate through block (use of min avoids bounds error)
					inds[sInner, n] = inds[sInner-1, n] + 1
				end
			end
		end
	end
	return(inds)
end
#Nonoverlapping block bootstrap (note, no block length detection routines available for this bootstrap method)
function dBootstrapIndex(bm::BootstrapNonoverlappingBlock, numObsData::Int, numObsResample::Int, numResample::Int)
	if numObsData < 3
		error("You must have at least 3 observations to use dependent bootstrap routines")
	end
	if bm.blockLength <= 0 || isnan(bm.blockLength)
		error("Expected block length is non-positive or NaN. Try a method that includes the underlying data (these methods typically auto-estimate the block length) or else supply the block length explicitly.")
	end
	if bm.blockLength == 1
		inds = dBootstrapIndex(BootstrapIID(), numObsData, numObsResample, numResample)
	else
		inds = Array(Int, numObsResample, numResample) #Pre-allocate output
		blockStartUB = max(1, numObsData-bm.blockLength+1)
		blockStartValues = [1:bm.blockLength:blockStartUB] #Build valid set of start indices for any block
		if blockStartValues[end] + bm.blockLength - 1 > numObsData
			if length(blockStartValues) > 1
				pop!(blockStartValues)
			end
		end
		for n = 1:numResample
			for s = 1:bm.blockLength:numObsResample
				inds[s, n] = blockStartValues[rand(1:length(blockStartValues))] #Start of block
				for sInner = s+1:min(s+bm.blockLength-1, numObsResample) #Iterate through block (use of min avoids bounds error)
					inds[sInner, n] = inds[sInner-1, n] + 1
				end
			end
		end
	end
	return(inds)
end
#Circular block bootstrap
function dBootstrapIndex(bm::BootstrapCircularBlock, numObsData::Int, numObsResample::Int, numResample::Int)
	if numObsData < 3
		error("You must have at least 3 observations to use dependent bootstrap routines")
	end
	if bm.blockLength <= 0 || isnan(bm.blockLength)
		error("Expected block length is non-positive or NaN. Try a method that includes the underlying data (these methods typically auto-estimate the block length) or else supply the block length explicitly.")
	end
	if bm.blockLength == 1
		inds = dBootstrapIndex(BootstrapIID(), numObsData, numObsResample, numResample)
	else
		inds = Array(Int, numObsResample, numResample) #Pre-allocate output
		for n = 1:numResample
			for s = 1:bm.blockLength:numObsResample
				inds[s, n] = rand(1:numObsData) #Start of block
				for sInner = s+1:min(s+bm.blockLength-1, numObsResample) #Iterate through block (use of min avoids bounds error)
					if inds[sInner-1, n] == numObsData #If this is true, then wrap back to 1
						inds[sInner, n] = 1
					else
						inds[sInner, n] = inds[sInner-1, n] + 1 #Otherwise, next index is one plus previous index
					end
				end
			end
		end
	end
	return(inds)
end
#TaperedBlock just calls moving block for this function (since Tapered block applies weighting to data obtained using movingBlock boostrap)
dBootstrapIndex(bm::BootstrapTaperedBlock, numObsData::Int, numObsResample::Int, numResample::Int) = dBootstrapIndex(BootstrapMovingBlock(bm.blockLength), numObsData, numObsResample, numResample)
#Basic wrapper methods
dBootstrapIndex(bm::BootstrapIID, numObsData::Int, numResample::Int) = dBootstrapIndex(bm, numObsData, numObsData, numResample)
dBootstrapIndex(bm::BootstrapIID, numObsData::Int) = dBootstrapIndex(bm, numObsData, numObsData, defaultNumResample)
dBootstrapIndex(bm::BootstrapStationary, numObsData::Int, numResample::Int) = dBootstrapIndex(bm, numObsData, numObsData, numResample)
dBootstrapIndex(bm::BootstrapStationary, numObsData::Int) = dBootstrapIndex(bm, numObsData, numObsData, defaultNumResample)
dBootstrapIndex(bm::BootstrapNonoverlappingBlock, numObsData::Int, numResample::Int) = dBootstrapIndex(bm, numObsData, numObsData, numResample)
dBootstrapIndex(bm::BootstrapNonoverlappingBlock, numObsData::Int) = dBootstrapIndex(bm, numObsData, numObsData, defaultNumResample)
dBootstrapIndex(bm::BootstrapMovingBlock, numObsData::Int, numResample::Int) = dBootstrapIndex(bm, numObsData, numObsData, numResample)
dBootstrapIndex(bm::BootstrapMovingBlock, numObsData::Int) = dBootstrapIndex(bm, numObsData, numObsData, defaultNumResample)
dBootstrapIndex(bm::BootstrapCircularBlock, numObsData::Int, numResample::Int) = dBootstrapIndex(bm, numObsData, numObsData, numResample)
dBootstrapIndex(bm::BootstrapCircularBlock, numObsData::Int) = dBootstrapIndex(bm, numObsData, numObsData, defaultNumResample)
dBootstrapIndex(bm::BootstrapTaperedBlock, numObsData::Int, numResample::Int) = dBootstrapIndex(bm, numObsData, numObsData, numResample)
dBootstrapIndex(bm::BootstrapTaperedBlock, numObsData::Int) = dBootstrapIndex(bm, numObsData, numObsData, defaultNumResample)
#Wrapper over a BootstrapParam
dBootstrapIndex(bp::BootstrapParam) = dBootstrapIndex(bp.bootstrapMethod, bp.numObsData, bp.numObsResample, bp.numResample)
#Wrappers that estimate block length from the data (but only if block length is non-positive)
function dBootstrapIndex{T<:Number}(x::Vector{T}, bp::BootstrapParam)
	if getBlockLength(bp) <= 0
		bpCopy = copy(bp) #This function should NOT update bp, so we need to make a copy of it
		dBootstrapBlockLength!(x, bpCopy)
		return(dBootstrapIndex(bpCopy))
	else
		return(dBootstrapIndex(bp))
	end
end
dBootstrapIndex{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrapIndex(x, bp)
#Keyword argument wrapper (block length must be supplied)
function dBootstrapIndex{T<:Number}(numObsData::Int, blockLength::T; numObsResample::Int=numObsData, numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary")
	if blockLength <= 0
		error("Input block length must be strictly positive")
	end
	bp = BootstrapParam(numObsData, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod)
	replaceBlockLength!(bp, blockLength)
	return(dBootstrapIndex(bp))
end
#Keyword argument wrapper (block length optionally estimated from input data if non-positive)
function dBootstrapIndex{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003")
	bp = BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLength=blockLength, blockLengthMethod=blockLengthMethod, bandwidthMethod=bandwidthMethod)
	if getBlockLength(bp) <= 0
		dBootstrapBlockLength(x, bp)
	end
	return(dBootstrapIndex(bp))
end
#Wrappers that update a BootstrapParam with block length in-place (but only if block length is non-positive)
function dBootstrapIndex!{T<:Number}(x::Vector{T}, bp::BootstrapParam)
	if getBlockLength(bp) <= 0
		dBootstrapBlockLength!(x, bp)
	end
	return(dBootstrapIndex(bp))
end
dBootstrapIndex!{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrapIndex!(x, bp)










#----------------------------------------------------------
#FUNCTION
#	dBootstrapWeight
#	dBootstrapWeight!
#INPUT
#	dBootstrapWeight accepts:
#		(bm::BootstrapTaperedBlock): Determine appropriate weightings to use with tapered-block bootstrap given the block length field in input bm
#	dBootstrapWeight! accepts
#		(x::Matrix, bm::BootstrapTaperedBlock): Apply weighting returned by dBootstrapWeight over bm in-place to matrix of re-sampled data x
#OUTPUT
#	dBootstrapWeight outputs a vector of weights to apply to re-sampled data
#	dBootstrapWeight! updates the input x in-place with the weights and so simply returns the boolean true
#PURPOSE
#	The purpose of dBootstrapWeight is to construct a weighting vector to apply to bootstrapped data. This occurs, for example, with the taperedBlock bootstrap
#	The purpose of dBootstrapWeight! is to apply the output of dBootstrapWeight to bootstrapped data in-place.
#NOTES
#	This function does not have any wrapper methods because it is not exported (i.e. it should only be used independently by those who know what they are doing)
#----------------------------------------------------------
function dBootstrapWeight(bm::BootstrapTaperedBlock)
	kernelInput = [ (1 / bm.blockLength) * (n - 0.5) for n = 1:bm.blockLength ]
	kernelWeight = evaluate(kernelInput, bm.kernelFunction)
	normTerm = sqrt(bm.blockLength) / norm(kernelWeight, 2)
	for n = 1:length(kernelWeight)
		kernelWeight *= normTerm
	end
	return(kernelWeight)
end
function dBootstrapWeight!(x::Matrix, bm::BootstrapTaperedBlock)
	if bm.blockLength == 1
		return(x)
	end
	w = dBootstrapWeight(bm)
	if length(w) != bm.blockLength
		error("Output of dBootstrapWeight function is wrong length")
	end
	S = div(size(x, 1), bm.blockLength) * bm.blockLength #Thus start of remainder block is at row S+1 (we don't apply weighting to remainder block)
	for n = 1:size(x, 2)
		c = 1
		for s = 1:S
			x[s, n] *= w[c]
			c += 1
			if c > length(w)
				c = 1
			end
		end
	end
	return(true)
end







#----------------------------------------------------------
#FUNCTION
#	dBootstrapData
#	dBootstrapData!
#INPUT
#	Let T<:Any. Function accepts:
#		(x::Vector{T}, bp::BootstrapParam): Obtain re-sampled data associated with original data vector x using the bootstrap procedure described in bp.
#		(bp::BootstrapParam, x::Vector{T}): Wrapper so order doesn't matter
#KEYWORD METHOD
#		{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003")
#IN-PLACE METHODS ie dBootstrapData!
#		(x::Vector{T}, bp::BootstrapParam): Identical to dBootstrapData but updates block length in bp in place if necessary
#		(bp::BootstrapParam, x::Vector{T}): Wrapper so order doesn't matter
#OUTPUT
#	Matrix{T} of re-sampled data
#PURPOSE
#	The purpose of this function is to build a matrix of bootstrapped data.
#NOTES
#	The statistic the user is bootstrapping is irrelevant to this function
#	Input type parameter is useful as it is used to determine the element type of the output matrix
#----------------------------------------------------------
function dBootstrapData{T<:Number}(x::Vector{T}, bp::BootstrapParam)
	if typeof(bp.bootstrapMethod) == BootstrapTaperedBlock
		if abs(mean(x)) > 1e-10
			error("Data must be de-meaned if bootstrap method is tapered block")
		end
	end
	bpCopy = copy(bp) #If we update block length, we don't want to affect original bp
	if getBlockLength(bpCopy) <= 0
		if typeof(bpCopy.bootstrapMethod) != BootstrapIID
			dBootstrapBlockLength!(x, bpCopy)
		end
	end
	xBoot = x[dBootstrapIndex(bpCopy)] #If bp.bootstrapMethod is BootstrapTaperedBlock, then bootstrap indices are generated using BootstrapMovingBlock
	if any(string(bpCopy.bootstrapMethod) .== methodsThatUseWeighting)
		dBootstrapWeight!(xBoot, bpCopy.bootstrapMethod) #Perform weighting on xBoot in-place
	end
	return(xBoot)
end
dBootstrapData{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrapData(x, bp)
#Keyword argument wrapper
dBootstrapData{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003") = dBootstrapData(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLength=blockLength, blockLengthMethod=blockLengthMethod, bandwidthMethod=bandwidthMethod))
#Update BootstrapParam in place wrapper
function dBootstrapData!{T<:Number}(x::Vector{T}, bp::BootstrapParam)
	if getBlockLength(bp) <= 0
		if typeof(bp.bootstrapMethod) != BootstrapIID
			dBootstrapBlockLength!(x, bp)
		end
	end
	xBoot = x[dBootstrapIndex(bp)] #If bp.bootstrapMethod is BootstrapTaperedBlock, then bootstrap indices are generated using BootstrapMovingBlock
	if any(string(bp.bootstrapMethod) .== methodsThatUseWeighting)
		dBootstrapWeight!(xBoot, bp.bootstrapMethod) #Perform weighting on xBoot in-place
	end
	return(xBoot)
end
dBootstrapData!{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrapData!(x, bp)







#----------------------------------------------------------
#FUNCTION
#	dBootstrapStatistic
#	dBootstrapStatistic!
#INPUT
#	Let T<:Any. Function accepts:
#		(x::Vector{T}, bp::BootstrapParam): Obtain re-sampled statistics using original data vector x and the bootstrap procedure described in bp.
#		(bp::BootstrapParam, x::Vector{T}): Wrapper so order doesn't matter
#KEYWORD METHOD
#		{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean")
#IN-PLACE METHODS ie dBootstrapStatistic!
#		(x::Vector{T}, bp::BootstrapParam): Identical to dBootstrapStatistic but updates block length in bp in place if necessary
#		(bp::BootstrapParam, x::Vector{T}): Wrapper so order doesn't matter
#OUTPUT
#	Vector{Float64} of bootstrapped statistics. Vector will be of length = numResample (a field in BootstrapParam)
#PURPOSE
#	The purpose of this function is to build a vector of bootstrapped values of the statistic of interest
#NOTES
#	Distribution parameter of statistic of interest is irrelevant to this function, as this will be estimated from the output of this function
#	This function employs the non-exported function dBootstrapStatistic_GetStatistic
#----------------------------------------------------------
dBootstrapStatistic{T<:Number}(x::Vector{T}, bp::BootstrapParam) = dBootstrapStatistic_GetStatistic(dBootstrapData(x, bp), bp)
dBootstrapStatistic{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrapStatistic(x, bp)
#Keyword argument wrapper
dBootstrapStatistic{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean") = dBootstrapStatistic(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLength=blockLength, blockLengthMethod=blockLengthMethod, bandwidthMethod=bandwidthMethod, statistic=statistic))
#Update BootstrapParam in place wrapper
dBootstrapStatistic!{T<:Number}(x::Vector{T}, bp::BootstrapParam) = dBootstrapStatistic_GetStatistic(dBootstrapData!(x, bp), bp)
dBootstrapStatistic!{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrapStatistic!(x, bp)
#Non-exported function used exclusively by dBootstrapStatistic to compute the actual statistics
function dBootstrapStatistic_GetStatistic{T<:Number}(d::Matrix{T}, bp::BootstrapParam)
	statVec = Array(Float64, bp.numResample)
	if typeof(bp.bootstrapMethod) == BootstrapTaperedBlock
		if typeof(bp.statistic) == ASCIIString
			if !(bp.statistic == "mean" || bp.statistic == "sum")
				error("Module is currently unable to perform the tapered block bootstrap for statistics other than the mean or sum")
			end
		else
			if !(bp.statistic == mean || bp.statistic == sum)
				error("Module is currently unable to perform the tapered block bootstrap for statistics other than the mean or sum")
			end
		end
	end
	if typeof(bp.statistic) == ASCIIString #bp.statistic is set to a recognized ASCIIString
		if bp.statistic == "mean"
			for n = 1:bp.numResample
				statVec[n] = convert(Float64, mean(d[:, n]))
			end
		elseif bp.statistic == "median"
			for n = 1:bp.numResample
				statVec[n] = convert(Float64, median(d[:, n]))
			end
		elseif bp.statistic == "variance" || bp.statistic == "var"
			for n = 1:bp.numResample
				statVec[n] = convert(Float64, var(d[:, n]))
			end
		elseif bp.statistic == "std"
			for n = 1:bp.numResample
				statVec[n] = convert(Float64, std(d[:, n]))
			end
		elseif bp.statistic == "sum"
			for n = 1:bp.numResample
				statVec[n] = convert(Float64, sum(d[:, n]))
			end
		elseif contains(bp.statistic, "quantile")
			length(bp.statistic) < 10 && error("Invalid string representation of quantile statistic to bootstrap")
			!(bp.statistic[1:9] == "quantile_") && error("Invalid string representation of quantile statistic to bootstrap")
			p = convert(Float64, parse("0." * bp.statistic[10:end]))
			!(0 < p < 1) && error("Invalid string representation of statistic to bootstrap")
			for n = 1:bp.numResample
				statVec[n] = convert(Float64, quantile(d[:, n], p))
			end
		else
			error("Invalid string representation of statistic to bootstrap")
		end
	else #bp.statistic is a function defined over Vector{Number} and returning a Float64 or type that can be converted
		statVec = Array(Float64, bp.numResample)
		for n = 1:bp.numResample
			statVec[n] = convert(Float64, bp.statistic(d[:, n]))
		end
	end
	return(statVec)
end











#----------------------------------------------------------
#FUNCTION
#	dBootstrap
#	dBootstrap!
#INPUT
#	Let T<:Any. Function accepts:
#		(x::Vector{T}, bp::BootstrapParam): Obtain bootstrapped distribution parameter of statistic of interest using original data vector x and the bootstrap procedure described in bp.
#		(bp::BootstrapParam, x::Vector{T}): Wrapper so order doesn't matter
#KEYWORD METHOD
#		{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean", distributionParam::ASCIIString="variance")
#IN-PLACE METHODS ie dBootstrap!
#		(x::Vector{T}, bp::BootstrapParam): Identical to dBootstrap but updates block length in bp in place if necessary
#		(bp::BootstrapParam, x::Vector{T}): Wrapper so order doesn't matter
#OUTPUT
#	Typically a Float64 for most distribution parameters specified in input BootstrapParam, but occassionally will be something else, e.g. a tuple (Float64, Float64) when the distribution parameter of interest is confidence bounds. Note, this implies the function is not type stable but that is not such a big deal since all actual work is farmed off to other type-stable functions.
#PURPOSE
#	The purpose of this function is to return the bootstrapped distribution parameter of the statistic of interest
#NOTES
#	These functions are not typically type stable since output might be confidence intervals (tuple) or variance (Float64). However, almost all the work is farmed out to other type stable functions so this is not really an issue.
#	This function employs the non-exported function dBootstrap_GetDistributionParam which has multiple methods for different distribution parameters
#----------------------------------------------------------
#Function to return bootstrapped parameter. Function is not type stable, hence it farms most of the work out to other functions
function dBootstrap{T<:Number}(x::Vector{T}, bp::BootstrapParam)
	if typeof(bp.distributionParam) == Function
		return(dBootstrap_GetDistributionParam(dBootstrapStatistic(x, bp), bp.distributionParam)) #Many possible output types
	elseif typeof(bp.distributionParam) == ASCIIString
		if contains(bp.distributionParam, "quantile")
			return(dBootstrap_GetDistributionParamQuantile(dBootstrapStatistic(x, bp), bp.distributionParam)) #Vector output
		elseif contains(bp.distributionParam, "conf")
			return(dBootstrap_GetDistributionParamConf(dBootstrapStatistic(x, bp), bp.distributionParam)) #Vector output
		else
			return(dBootstrap_GetDistributionParam(dBootstrapStatistic(x, bp), bp.distributionParam)) #Float64 output
		end
	else
		error("Logic fail. Invalid type in distributionParam field of BootstrapParam")
	end
end
dBootstrap{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrap(x, bp)
#Keyword argument wrapper
dBootstrap{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean", distributionParam::ASCIIString="variance") = dBootstrap(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLength=blockLength, blockLengthMethod=blockLengthMethod, bandwidthMethod=bandwidthMethod, statistic=statistic, distributionParam=distributionParam))
#Wrapper that updates BootstrapParam in-place
function dBootstrap!{T<:Number}(x::Vector{T}, bp::BootstrapParam)
	if getBlockLength(bp) <= 0
		dBootstrapBlockLength!(x, bp)
	end
	return(dBootstrap(x, bp))
end
dBootstrap!{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dBootstrap!(x, bp)
#Method for arbitrary function input
dBootstrap_GetDistributionParam(s::Vector{Float64}, f::Function) = f(s)
#Method for pre-set ASCIIString input (excluding quantiles and confidence intervals)
function dBootstrap_GetDistributionParam(s::Vector{Float64}, fStr::ASCIIString)
	if fStr == "mean"
		return(mean(s))
	elseif fStr == "median"
		return(median(s))
	elseif fStr == "variance" || fStr == "var"
		return(var(s))
	elseif fStr == "std"
		return(std(s))
	else
		error("Invalid string representation of distribution parameter function")
	end
end
#Function for pre-set ASCIIString input for quantiles
function dBootstrap_GetDistributionParamQuantile(s::Vector{Float64}, fStr::ASCIIString)
	fStrVec = split(fStr, '_')
	length(fStrVec) < 2 && error("Invalid string representation of distribution parameter quantile")
	fStrVec[1] != "quantile" && error("Invalid string representation of distribution parameter quantile")
	pVec = Array(Float64, length(fStrVec)-1)
	for k = 2:length(fStrVec)
		pVec[k-1] = convert(Float64, parse("0." * fStrVec[k]))
		!(0 < pVec[k-1] < 1) && error("Invalid string representation of distribution parameter quantile")
	end
	return(quantile(s, pVec))
end
function dBootstrap_GetDistributionParamConf(s::Vector{Float64}, fStr::ASCIIString)
	if fStr == "conf"
		return(quantile(s, [0.025, 0.975]))	#Default behaviour of 95% confidence bound
	else
		fStrVec = split(fStr, '_')
		length(fStrVec) != 3 && error("Invalid string representation of distribution parameter confidence interval")
		fStrVec[1] != "conf" && error("Invalid string representation of distribution parameter confidence interval")
		pVec = [convert(Float64, parse("0." * fStrVec[2])), convert(Float64, parse("0." * fStrVec[3]))]
		!(0 < pVec[1] < 1) && error("Distribution parameter confidence interval lower probability outside (0, 1)")
		!(0 < pVec[2] < 1) && error("Distribution parameter confidence interval upper probability outside (0, 1)")
		pVec[1] >= pVec[2] && error("Distribution parameter confidence interval lower bound >= upper bound")
		return(quantile(s, pVec))
	end
end







#----------------------------------------------------------
#FUNCTION
#	dBootstrapVar
#	dBootstrapStd
#	dBootstrapConf
#INPUT
#	Let T<:Any. All functions accept the following keyword method:
#		{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean")#		(x::Vector{T}, param::BootstrapParam): Obtain bootstrap data of x where bootstrap methodology is determined by param
#OUTPUT
#	Output is Float64 for dBootstrapVar and dBootstrapStd which is the distribution parameter (variance and standard deviation respectively) of the statistic of interest. Output is a tuple of confidence bounds for dBootstrapConf
#PURPOSE
#	These functions are merely wrappers that express the name of the distribution parameter of the statistic of interest that we want to bootstrap. Var, Std, and Conf refer to variance, standard deviation, and 95% confidence interval respectively.
#NOTES
#----------------------------------------------------------
dBootstrapVar{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean") = dBootstrap(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLength=blockLength, blockLengthMethod=blockLengthMethod, bandwidthMethod=bandwidthMethod, statistic=statistic, distributionParam="variance")
dBootstrapStd{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean") = dBootstrap(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLength=blockLength, blockLengthMethod=blockLengthMethod, bandwidthMethod=bandwidthMethod, statistic=statistic, distributionParam="std")
dBootstrapConf{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean") = dBootstrap(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLength=blockLength, blockLengthMethod=blockLengthMethod, bandwidthMethod=bandwidthMethod, statistic=statistic, distributionParam="conf")







end # module
