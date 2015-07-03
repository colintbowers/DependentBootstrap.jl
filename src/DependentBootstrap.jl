module DependentBootstrap
#-----------------------------------------------------------
#PURPOSE
#	Colin T. Bowers module for dependent bootstraps
#NOTES
#	Throughout this module, automatic block length estimation (and hence the type BlockLengthMethod) is only ever used if the current block length is <= 0 (note, default block lengths are -1)
#	Note, all keyword methods of the functions in this module are funnelled through the BootstrapParam keyword constructor.
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
		Base.copy,
		Base.deepcopy

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
		getblocklength,
		update!,
		dbootstrapblocklength,
		dbootstrapblocklength!,
		dbootstrapindex,
		dbootstrapindex!,
		dbootstrapdata,
		dbootstrapdata!,
		dbootstrapstatistic,
		dbootstrapstatistic!,
		dbootstrap,
		dbootstrap!,
		dbootstrap_mean,
		dbootstrap_mean!,
		dbootstrap_median,
		dbootstrap_median!,
		dbootstrap_var,
		dbootstrap_var!,
		dbootstrap_std,
		dbootstrap_std!,
		dbootstrap_quantile,
		dbootstrap_quantile!,
		dbootstrap_conf,
		dbootstrap_conf!




#******************************************************************************

#----------------------------------------------------------
#SET CONSTANTS FOR MODULE
#----------------------------------------------------------
const defaultNumResample = 1000::Int #Sets the default number of re-samples used in BootstrapParam constructors
const defaultBlockLength = -1::Int #A default block length that is useful for constructors that initialize a BootstrapParam without estimating the block length
const methodsThatUseWeighting = [:taperedBlock]::Vector{Symbol} #Vector of bootstrap methods that necessitate weighting (ie a call to dbootstrapweight!)



#Two abstract types to nest all bootstrap method types and block length selection method types
abstract BootstrapMethod
abstract BlockLengthMethod


#---------- BOOTSTRAP METHOD TYPES ------------------
#PURPOSE: This type stores the information needed to select and run a particular statistical bootstrap method
#NOTES
#	Block length values can be non-positive. This is understood by the functions in this module to imply that the block-length needs to be detected.
#----------------------------------------------------------
#---- TYPE DEFINTIIONS ----
type Bootstrap_Dummy <: BootstrapMethod ; end #Dummy type useful for update! functions
type BootstrapIID <: BootstrapMethod #This bootstrap method does not have any associated parameters
end
BootstrapIID(blockLength::Number) = BootstrapIID()
type BootstrapStationary <: BootstrapMethod
	expectedBlockLength::Float64
end
BootstrapStationary() = BootstrapStationary(-1.0)
BootstrapStationary(expectedBlockLength::Number) = BootstrapStationary(convert(Float64, expectedBlockLength))
type BootstrapMovingBlock <: BootstrapMethod
	blockLength::Int
end
BootstrapMovingBlock() = BootstrapMovingBlock(-1)
BootstrapMovingBlock(blockLength::Number) = BootstrapMovingBlock(convert(Int, ceil(blockLength)))
type BootstrapNonoverlappingBlock <: BootstrapMethod
	blockLength::Int
end
BootstrapNonoverlappingBlock() = BootstrapNonoverlappingBlock(-1)
BootstrapNonoverlappingBlock(blockLength::Number) = BootstrapNonoverlappingBlock(convert(Int, ceil(blockLength)))
type BootstrapCircularBlock <: BootstrapMethod
	blockLength::Int
end
BootstrapCircularBlock() = BootstrapCircularBlock(-1)
BootstrapCircularBlock(blockLength::Number) = BootstrapCircularBlock(convert(Int, ceil(blockLength)))
#Tapered block bootstrap
type BootstrapTaperedBlock <: BootstrapMethod
	blockLength::Int
	kernelFunction::KernelFunction
	function BootstrapTaperedBlock(blockLength::Int, kernelFunction::KernelFunction)
		!(typeof(kernelFunction) == KernelPP2002Trap || typeof(kernelFunction) == KernelPP2002Smooth) && error("Invalid kernel function for tapered block bootstrap method")
		new(blockLength, kernelFunction)
	end
end
BootstrapTaperedBlock() = BootstrapTaperedBlock(-1, KernelPP2002Trap())
BootstrapTaperedBlock(blockLength::Int) = BootstrapTaperedBlock(blockLength, KernelPP2002Trap())
BootstrapTaperedBlock(blockLength::Number) = BootstrapTaperedBlock(convert(Int, ceil(blockLength)), KernelPP2002Trap())
#---- METHODS ----
#deepcopy
function deepcopy(x::BootstrapMethod)
	tempArgs = [ deepcopy(getfield(x, i)) for i = 1:length(names(x)) ]
	return(eval(parse(string(typeof(x)) * "(tempArgs...)")))
end
#string
string(x::Bootstrap_Dummy) = "dummy"
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
function getblocklength(bm::BootstrapMethod)
	typeof(bm) == BootstrapIID && return(1.0)
	typeof(bm) == BootstrapStationary && return(bm.expectedBlockLength)
	return(convert(Float64, bm.blockLength))
end
#Update fields of a bootstrap method
update!(bm::BootstrapIID; blockLength::Number=-999) = true
function update!(bm::BootstrapStationary; blockLength::Number=-999)
	blockLength != -999 && (bm.expectedBlockLength = convert(Float64, blockLength))
	return(bm)
end
function update!{T<:Union(BootstrapMovingBlock, BootstrapCircularBlock, BootstrapNonoverlappingBlock, BootstrapTaperedBlock)}(bm::T; blockLength::Number=-999)
	blockLength != -999 && (bm.blockLength = convert(Int, blockLength))
	return(bm)
end
function update!(bm::BootstrapTaperedBlock; blockLength::Number=-999, kernelFunction::KernelFunction=KernelDummy())
	blockLength != -999 && (bm.blockLength = convert(Int, blockLength))
	typeof(kernelFunction) != KernelDummy && (bm.kernelFunction = kernelFunction)
	return(bm)
end




#------- BLOCK LENGTH METHOD TYPES -----------------------
#PURPOSE: This type stores the information needed to select and run a particular block length selection method
#NOTES
#----------------------------------------------------------
#----------- TYPE DEFINITIONS -----------
#Dummy block length method (useful for update! functions)
type BlockLength_Dummy <: BlockLengthMethod; end
#Block length detection procedure of Patton, Politis, White (2009)
type BlockLengthPPW2009 <: BlockLengthMethod
	bandwidthMethod::BandwidthMethod
	bootstrapMethod::Symbol
	function BlockLengthPPW2009{T<:BandwidthMethod}(bandwidthMethod::T, bootstrapMethod::Symbol)
		if !(bootstrapMethod == :stationary || bootstrapMethod == :circularBlock || bootstrapMethod == :movingBlock)
			error("BlockLengthPPW2009 detection method only viable for stationary, circular block, and moving block bootstraps")
		end
		new(bandwidthMethod, bootstrapMethod)
	end
end
BlockLengthPPW2009() = BlockLengthPPW2009(BandwidthP2003(), :stationary)
BlockLengthPPW2009(numObs::Int) = BlockLengthPPW2009(BandwidthP2003(numObs), :stationary)
BlockLengthPPW2009(numObs::Int, bootstrapMethod::Symbol) = BlockLengthPPW2009(BandwidthP2003(numObs), bootstrapMethod)
BlockLengthPPW2009(bandwidthMethod::BandwidthMethod) = BlockLengthPPW2009(bandwidthMethod, :stationary)
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
BlockLengthPP2002(numObs::Int, kernelFunction::KernelFunction) = BlockLengthPP2002(BandwidthP2003(numObs), kernelFunction)
BlockLengthPP2002(bandwidthMethod::BandwidthMethod) = BlockLengthPP2002(bandwidthMethod, KernelPP2002Trap())
#------------ METHODS -------------
#copy function
function deepcopy(x::BlockLengthMethod)
	tempArgs = [ deepcopy(getfield(x, i)) for i = 1:length(names(x)) ]
	return(eval(parse(string(typeof(x)) * "(tempArgs...)")))
end
#string function for string representation of block length method
string(b::BlockLength_Dummy) = "dummy"
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
function update!(bm::BlockLengthPPW2009; bandwidthMethod::BandwidthMethod=BandwidthDummy(), bootstrapMethod::Symbol=:none)
	typeof(bandwidthMethod) != BandwidthDummy && (bm.bandwidthMethod = bandwidthMethod)
	if bootstrapMethod != :none
		!(bootstrapMethod == :stationary || bootstrapMethod == :circularBlock || bootstrapMethod == :movingBlock) && error("BlockLengthPPW2009 detection method only viable for stationary, circular block, and moving block bootstraps")
		bm.bootstrapMethod = bootstrapMethod
	end
	return(bm)
end
function update!(bm::BlockLengthPP2002; bandwidthMethod::BandwidthMethod=BandwidthDummy(), kernelFunction::KernelFunction=KernelDummy())
	typeof(bandwidthMethod) != BandwidthDummy && (bm.bandwidthMethod = bandwidthMethod)
	if typeof(kernelFunction) != KernelDummy
		!(typeof(kernelFunction) == KernelPP2002Trap || typeof(kernelFunction) == KernelPP2002Smooth) && error("Tapered block bootstrap method can only be implemented with KernelPP2002Trap or KernelPP2002Smooth kernel functions")
		bm.kernelFunction = kernelFunction
	end
	return(bm)
end



#----------------------------------------------------------
#TYPE
#   BootstrapParam
#FIELDS
#	Fields are:
#		numObsData::Int: Number of observations in the dataset that is to be bootstrapped
#		numObsResample::Int: Number of observations to be drawn (with replacement) per resample. For the vast majority of cases this field will equal numObsData.
#		numResample::Int: Number of bootstrap resamples
#		bootstrapMethod::BootstrapMethod: Bootstrap method to use (see bootstrap method type definitions for more detail)
#		blockLengthMethod::BlockLengthMethod: Method to use to detect block length (only applies if current block length is non-positive) (see block length method types for more detail)
#		statistic::Function: The statistic of interest. Must be a function that is able to be applied to vector of observations and return a single value.
#		distributionParam::Function: The distribution parameter of the distribution of the statistic of interest. Must be a function that is able to be applied to a vector of observations.
#PURPOSE
#	This type stores the information needed to select and run a particular statistical bootstrap method. Almost every function in this module accepts this types as input.
#NOTES
#	In the vast majority of cases, numObsData = numObsResample. The constructors reflect this.
#	The block length is a field of bootstrapMethod. It is understood that any block length <= 0 implies that the user wants the function to auto-detect the block length
#	The statistic field is a function which implies type instability. However, the most popular choices for this field have been hard-coded to an if statement and an explicit call to the underlying so there is no performance overhead. Use options to take advantage of this include:
#		*mean
#		*median
#		*var
#		*std
#		*DependentBootstrap.quantile_x (where x can take value 001, 01, 05, 1, 9, 95, 99, 999. These are interpreted as probabilities with a 0. in front of them)
#		*sum
#	The distributionParam field is a function which implies type instability. Hence dbootstrap and dbootstrap! are NOT TYPE STABLE. However, specialised versions of this function that are type stable are provided and take the form dbootstrap_x where x is the function name.
#----------------------------------------------------------
#----------- TYPE DEFINITION --------------------
type BootstrapParam
	numObsData::Int #Number of observations in the dataset that is to be bootstrapped
	numObsResample::Int #Number of observations to be drawn (with replacement) per resample
	numResample::Int #Number of resamples
	bootstrapMethod::BootstrapMethod #Bootstrap method to use
	blockLengthMethod::BlockLengthMethod #Method to use to detect block length (only applies if current block length is non-positive)
	statistic::Function #Statistic of interest
	distributionParam::Function #Parameter of distribution of statistic of interest
	function BootstrapParam(numObsData::Int, numObsResample::Int, numResample::Int, bootstrapMethod::BootstrapMethod, blockLengthMethod::BlockLengthMethod, statistic::Function, distributionParam::Function)
		numObsData < 1 && error("Number of observations in dataset must be greater than zero to apply a bootstrap")
		numObsResample < 1 && error("Number of observations per resample must be greater than zero")
		numResample < 1 && error("Number of resamples to bootstrap must be greater than zero")
		new(numObsData, numObsResample, numResample, bootstrapMethod, blockLengthMethod, statistic, distributionParam)
	end
end
#------------ OUTER CONSTRUCTORS ----------------
#Keyword constructor with no data provided
function BootstrapParam(numObsData::Int; numObsResample::Int=numObsData, numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean, distributionParam::Function=var)
	if convert(Int, blockLength) != -1
		blockLength <= 0 && error("blockLength keyword argument must be strictly positive")
		typeof(blockLengthMethod) != BlockLength_Dummy && error("There is no need to specify a blockLengthMethod if you have already specified a block length")
		blockLengthMethod = BlockLengthPPW2009(numObsData)
		update!(bootstrapMethod, blockLength=blockLength)
	else
		if typeof(blockLengthMethod) == BlockLength_Dummy #We need to automatically decide a block-length procedure
			if typeof(bootstrapMethod) == BootstrapStationary; blockLengthMethod = BlockLengthPPW2009(BandwidthP2003(numObsData), :stationary)
			elseif typeof(bootstrapMethod) == BootstrapCircularBlock; blockLengthMethod = BlockLengthPPW2009(BandwidthP2003(numObsData), :circularBlock)
			elseif typeof(bootstrapMethod) == BootstrapMovingBlock; blockLengthMethod = BlockLengthPPW2009(BandwidthP2003(numObsData), :movingBlock)
			elseif typeof(bootstrapMethod) == BootstrapTaperedBlock; blockLengthMethod = BlockLengthPP2002(BandwidthP2003(numObsData), KernelPP2002Trap())
			else; error("Automatic blockLengthMethod selection not possible for chosen bootstrapMethod. Please explicitly specify either a blockLength or blockLengthMethod.")
			end
		end
	end
	return(BootstrapParam(numObsData, numObsResample, numResample, bootstrapMethod, blockLengthMethod, statistic, distributionParam))
end
#Keyword constructor with data provided (incorporates automatic block length selection if block length not provided)
function BootstrapParam{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean, distributionParam::Function=var)
	bp = BootstrapParam(length(x), numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic, distributionParam=distributionParam)
	blockLength <= 0 && dbootstrapblocklength!(x, bp) #If no block-length provided, detect it from the data
	return(bp)
end
#------------ METHODS ------------------
#deepcopy
function deepcopy(x::BootstrapParam)
	tempArgs = [ deepcopy(getfield(x, i)) for i = 1:length(names(x)) ]
	return(eval(parse(string(typeof(x)) * "(tempArgs...)")))
end
#show methods
function show(io::IO, b::BootstrapParam)
	println(io, "dependent bootstrap parameter values:")
	println(io, "    number of observations in dataset = " * string(b.numObsData))
	println(io, "    number of observations per resample = " * string(b.numObsResample))
	println(io, "    number of resamples = " * string(b.numResample))
	println(io, "    bootstrap method = " * string(b.bootstrapMethod))
	println(io, "    block length method = " * string(b.blockLengthMethod))
	println(io, "    current block length = " * string(getblocklength(b)))
	println(io, "    statistic of interest = " * string(b.statistic))
	println(io, "    distribution parameter of statistic = " * string(b.distributionParam))
end
#show wrapper on STDOUT
show(b::BootstrapParam) = show(STDOUT, b)
#Retrieve the block length
getblocklength(bp::BootstrapParam) = getblocklength(bp.bootstrapMethod)
#Update specified fields with specified keyword arguments
function update!(bp::BootstrapParam; numObsData::Int=-999, numObsResample::Int=-999, numResample::Int=-999, bootstrapMethod::BootstrapMethod=Bootstrap_Dummy(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), statistic::Function=bootstrap_dummy_func, distributionParam::Function=bootstrap_dummy_func, blockLength::Number=-999)
	if numObsData != -999
		numObsData < 1 && error("Number of observations in dataset must be greater than zero to apply a bootstrap")
		bp.numObsData = numObsData
	end
	if numObsResample != -999
		numObsResample < 1 && error("Number of observations per resample must be greater than zero")
		bp.numObsResample = numObsResample
	end
	if numResample != -999
		bp.numResample = numResample
		numResample < 1 && error("Number of resamples to bootstrap must be greater than zero")
	end
	typeof(bootstrapMethod) != Bootstrap_Dummy && (bp.bootstrapMethod = bootstrapMethod)
	typeof(blockLengthMethod) != BlockLength_Dummy && (bp.blockLengthMethod = blockLengthMethod)
	statistic != bootstrap_dummy_func && (bp.statistic = statistic)
	distributionParam != bootstrap_dummy_func && (bp.distributionParam = distributionParam)
	blockLength != -999 && update!(bp.bootstrapMethod, blockLength=blockLength)
	return(bp)
end
bootstrap_dummy_func() = true #Don't delete this. It is used as a default input to update!





#----------------------------------------------------------
#FUNCTION
#	dbootstrapblocklength
#	dbootstrapblocklength!
#INPUT
#OUTPUT
#	Estimated block length for dependent bootstrap procedure, expressed as Float64.
#PURPOSE
#	The purpose of this function is to estimate the block length to use with a given dependent bootstrap procedure.
#NOTES
#----------------------------------------------------------
#Block length selection method of Patton, Politis, White (2009) "Correction to Automatic Block Length Selection For the Dependent Bootstrap"
function dbootstrapblocklength{T<:Number}(x::AbstractVector{T}, bm::BlockLengthPPW2009)
	length(x) < 3 && error("You must have at least 3 observations to estimate block length")
	(M, xVar, covVec) = dbootstrapblocklength_MAndCorVec(x, bm)
	kernelCovVec = dbootstrapblocklength_KernelCovVec(covVec, KernelP2003FlatTop(), M)
	gHat = 0.0
	for k = 1:M
		gHat += 2 * k * kernelCovVec[k] #note, "2*" is because sum is from -M and M, but is symmetric about 0. Note, 0 term in sum vanishes since k=0 -> |k|=0
	end
	if bm.bootstrapMethod == :stationary
		dHat = 2 * (xVar + 2*sum(kernelCovVec))^2 #note, in expression (1 + 2*sum(kernelCovVec)), "2*" is because sum is from -M to M, but is symmetric around 0. "1+" is the zero term of the sum which is simply equal to unity.
	elseif bm.bootstrapMethod == :circularBlock || bm.bootstrapMethod == :movingBlock
		dHat = (4/3) * (xVar + 2*sum(kernelCovVec))^2 #note, in expression (1 + 2*sum(kernelCovVec)), "2*" is because sum is from -M to M, but is symmetric around 0. "1+" is the zero term of the sum which is simply equal to unity.
	else
		error("Optimal block length given BlockLengthPPW2009 method is only known for stationary, circular, or moving block bootstraps")
	end
	blockLength = (2 * gHat^2 * length(x) / dHat)^(1/3)
	blockLength = min(blockLength, ceil(min(3*sqrt(length(x)), length(x) / 3))) #Enforce upper bound on block length suggested by Patton
	return(blockLength) #Equation 9 and 14 from Politis, White (2004)
end
#Block length selection method of Paparoditis, Politis (2002) "The Tapered Block Bootstrap for General Statistics From Stationary Sequences"
function dbootstrapblocklength{T<:Number}(x::AbstractVector{T}, bm::BlockLengthPP2002)
	length(x) < 3 && error("You must have at least 3 observations to estimate block length")
	(M, xVar, covVec) = dbootstrapblocklength_MAndCorVec(x, bm)
	kernelCovVec = dbootstrapblocklength_KernelCovVec(covVec, KernelP2003FlatTop(), M)
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
#Non-exported helper functions
function dbootstrapblocklength_MAndCorVec{T<:Number}(x::AbstractVector{T}, bm::BlockLengthMethod)
	!(typeof(bm) == BlockLengthPPW2009 || typeof(bm) == BlockLengthPP2002) && error("This non-exported function is not defined for the input block length method type")
	(M, xVar, covVec) = bandwidth(x, bm.bandwidthMethod)
	if M > 0.5 * length(x)
		println("WARNING: Bandwidth in parameter estimation section of dbootstrapblocklength forced to half total number of observations. Data may contain excessive dependence.")
		M = convert(Int, round(0.5 * length(x)))
	end
	M < 2 && (M = 2) #Even though M output of bandwidth function is always greater than 2, the above check for excessively large M can result in M being reduced below 2 (admittedly only in very unusual circumstances)
	length(covVec) < M && append!(covVec, autocov(x, length(covVec)+1:M)) #Get any additional autocovariances that we might need
	return(M, xVar, covVec)
end
function dbootstrapblocklength_KernelCovVec(covVec::AbstractVector{Float64}, kF::KernelFunction, M::Int)
	length(covVec) < M && error("Logic fail in module. It should have been impossible for length(covVec) < M")
	kernelCovVec = [ evaluate(k / M, kF) * covVec[k] for k = 1:M ]
	return(kernelCovVec)
end
#Keyword wrapper
function dbootstrapblocklength{T<:Number}(x::AbstractVector{T}; blockLengthMethod::BlockLengthMethod=BlockLengthPPW2009(), bandwidthMethod::BandwidthMethod=BandwidthP2003())
	update!(blockLengthMethod, bandwidthMethod=bandwidthMethod)
	return(dbootstrapblocklength(x, blockLengthMethod))
end
#Wrapper methods
function dbootstrapblocklength{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam)
	getblocklength(bp) > 0 && error("Your block length is already > 0. Either reset block length to < 0 or else use a different method of this function")
	return(dbootstrapblocklength(x, bp.blockLengthMethod))
end
dbootstrapblocklength{T<:Number}(bp::BootstrapParam, x::AbstractVector{T}) = dbootstrapblocklength(x, bp)
#In-place update of BootstrapParam version of function
function dbootstrapblocklength!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam)
	getblocklength(bp) > 0 && error("Your block length is already > 0. Either reset block length to < 0 or else use a different method of this function")
	update!(bp.bootstrapMethod, blockLength=dbootstrapblocklength(x, bp.blockLengthMethod))
end
dbootstrapblocklength!{T<:Number}(bp::BootstrapParam, x::AbstractVector{T}) = dbootstrapblocklength!(x, bp)





#----------------------------------------------------------
#FUNCTION
#	dbootstrapindex
#	dbootstrapindex!
#INPUT
#OUTPUT
#	Matrix{Int} of indices that can be used to index into the data and construct re-sampled data.
#PURPOSE
#	The purpose of this function is to provide a set of bootstrap indices that can be used to index into a data vector to create re-sampled data
#NOTES
#----------------------------------------------------------
#iid bootstrap
function dbootstrapindex(bm::BootstrapIID, numObsData::Int, numObsResample::Int, numResample::Int)
	numObsData < 3 && error("You must have at least 3 observations to use dependent bootstrap routines")
	return(rand(1:numObsData, numObsResample, numResample))
end
#Stationary bootstrap
function dbootstrapindex(bm::BootstrapStationary, numObsData::Int, numObsResample::Int, numResample::Int)
	numObsData < 3 && error("You must have at least 3 observations to use dependent bootstrap routines")
	bm.expectedBlockLength <= 0.0 || isnan(bm.expectedBlockLength) && error("Expected block length is non-positive or NaN. Try a method that includes the underlying data (these methods typically auto-estimate the block length) or else supply the block length explicitly.")
	if bm.expectedBlockLength <= 1.0 #IIDBootstrap case
		inds = dbootstrapindex(BootstrapIID(), numObsData, numObsResample, numResample)
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
function dbootstrapindex(bm::BootstrapMovingBlock, numObsData::Int, numObsResample::Int, numResample::Int)
	numObsData < 3 && error("You must have at least 3 observations to use dependent bootstrap routines")
	bm.blockLength <= 0 || isnan(bm.blockLength) && error("Expected block length is non-positive or NaN. Try a method that includes the underlying data (these methods typically auto-estimate the block length) or else supply the block length explicitly.")
	if bm.blockLength == 1
		inds = dbootstrapindex(BootstrapIID(), numObsData, numObsResample, numResample)
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
function dbootstrapindex(bm::BootstrapNonoverlappingBlock, numObsData::Int, numObsResample::Int, numResample::Int)
	numObsData < 3 && error("You must have at least 3 observations to use dependent bootstrap routines")
	bm.blockLength <= 0 || isnan(bm.blockLength) && error("Expected block length is non-positive or NaN. Try a method that includes the underlying data (these methods typically auto-estimate the block length) or else supply the block length explicitly.")
	if bm.blockLength == 1
		inds = dbootstrapindex(BootstrapIID(), numObsData, numObsResample, numResample)
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
function dbootstrapindex(bm::BootstrapCircularBlock, numObsData::Int, numObsResample::Int, numResample::Int)
	numObsData < 3 && error("You must have at least 3 observations to use dependent bootstrap routines")
	bm.blockLength <= 0 || isnan(bm.blockLength) && error("Expected block length is non-positive or NaN. Try a method that includes the underlying data (these methods typically auto-estimate the block length) or else supply the block length explicitly.")
	if bm.blockLength == 1
		inds = dbootstrapindex(BootstrapIID(), numObsData, numObsResample, numResample)
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
dbootstrapindex(bm::BootstrapTaperedBlock, numObsData::Int, numObsResample::Int, numResample::Int) = dbootstrapindex(BootstrapMovingBlock(bm.blockLength), numObsData, numObsResample, numResample)
#BootstrapParam wrapper
dbootstrapindex(bp::BootstrapParam) = dbootstrapindex(bp.bootstrapMethod, bp.numObsData, bp.numObsResample, bp.numResample)
#Keyword wrapper with block length provided
function dbootstrapindex(numObsData::Int, blockLength::Number; bootstrapMethod::BootstrapMethod=BootstrapStationary(), numObsResample::Int=numObsData, numResample::Int=defaultNumResample)
	update!(bootstrapMethod, blockLength=blockLength)
	return(dbootstrapindex(bootstrapMethod, numObsData, numObsResample, numResample))
end
#Keyword wrapper with data provided (let BootstrapParam constructor do the work)
function dbootstrapindex{T<:Number}(x::AbstractVector{T}; blockLength::Number=-1, bootstrapMethod::BootstrapMethod=BootstrapStationary(), numObsResample::Int=numObsData, numResample::Int=defaultNumResample)
	bp = BootstrapParam(x, blockLength=blockLength, bootstrapMethod=bootstrapMethod, numObsResample=numObsResample, numResample=numResample)
	return(dbootstrapindex(bp))
end
#Wrappers that update a BootstrapParam with block length in-place (but only if block length is non-positive)
function dbootstrapindex!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam)
	getblocklength(bp) <= 0 && update!(bp.bootstrapMethod, blockLength=dbootstrapblocklength(x, bp))
	return(dbootstrapindex(bp))
end
dbootstrapindex!{T<:Number}(bp::BootstrapParam, x::Vector{T}) = dbootstrapindex!(x, bp)





#These two non-exported functions are used to weight bootstrapped observations. Thus far, only the tapered block bootstrap uses this function.
function dbootstrapweight(bm::BootstrapTaperedBlock)
	kernelInput = [ (1 / bm.blockLength) * (n - 0.5) for n = 1:bm.blockLength ]
	kernelWeight = evaluate(kernelInput, bm.kernelFunction)
	normTerm = sqrt(bm.blockLength) / norm(kernelWeight, 2)
	for n = 1:length(kernelWeight)
		kernelWeight *= normTerm
	end
	return(kernelWeight)
end
function dbootstrapweight!(x::Matrix, bm::BootstrapTaperedBlock)
	if bm.blockLength == 1
		return(x)
	end
	w = dbootstrapweight(bm)
	if length(w) != bm.blockLength
		error("Output of dbootstrapweight function is wrong length")
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
#	dbootstrapdata
#	dbootstrapdata!
#PURPOSE
#	The purpose of this function is to build a matrix of bootstrapped data.
#NOTES
#	The statistic the user is bootstrapping is irrelevant to this function
#----------------------------------------------------------
function dbootstrapdata!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam)
	if typeof(bp.bootstrapMethod) == BootstrapTaperedBlock #Check to make sure tapered block bootstrap is being used correctly
		abs(mean(x)) > 1e-10 && error("Data must be de-meaned if bootstrap method is tapered block")
	end
	if getblocklength(bp) <= 0 #Update block length if necessary
		typeof(bp.bootstrapMethod) != BootstrapIID && dbootstrapblocklength!(x, bp)
	end
	xBoot = x[dbootstrapindex(bp)]
	typeof(bp.bootstrapMethod) == BootstrapTaperedBlock && 	dbootstrapweight!(xBoot, bp.bootstrapMethod) #tapered block requires weighting
	return(xBoot)
end
dbootstrapdata!{T<:Number}(bp::BootstrapParam, x::AbstractVector{T}) = dbootstrapdata!(x, bp)
dbootstrapdata{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = dbootstrapdata!(x, deepcopy(bp))
dbootstrapdata{T<:Number}(bp::BootstrapParam, x::AbstractVector{T}) = dbootstrapdata(x, bp)
#Keyword argument wrapper
dbootstrapdata{T<:Number}(x::AbstractVector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1) = dbootstrapdata(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength))





#----------------------------------------------------------
#FUNCTION
#	dbootstrapstatistic
#	dbootstrapstatistic!
#PURPOSE
#	The purpose of this function is to build a vector of bootstrapped values of the statistic of interest (expressed as Vector{Float64})
#----------------------------------------------------------
dbootstrapstatistic!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = dbootstrapstatistic_getstatistic(dbootstrapdata!(x, bp), bp)
dbootstrapstatistic!{T<:Number}(bp::BootstrapParam, x::AbstractVector{T}) = dbootstrapstatistic!(x, bp)
dbootstrapstatistic{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = dbootstrapstatistic_getstatistic(dbootstrapdata(x, bp), bp)
dbootstrapstatistic{T<:Number}(bp::BootstrapParam, x::AbstractVector{T}) = dbootstrapstatistic(x, bp)
#Keyword argument wrapper
dbootstrapstatistic{T<:Number}(x::AbstractVector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean) = dbootstrapstatistic(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic))
#Non-exported function used exclusively by dbootstrapstatistic to compute the actual statistics
function dbootstrapstatistic_getstatistic{T<:Number}(d::Matrix{T}, bp::BootstrapParam)
	if typeof(bp.bootstrapMethod) == BootstrapTaperedBlock
		!(bp.statistic == mean || bp.statistic == sum) && error("Module is currently unable to perform the tapered block bootstrap for statistics other than the mean or sum")
	end
	N = size(d, 1)
	#Julia is slow when functions are passed round as variables. The following if statement checks whether the user specified function is a popular choice and if so, that function is called explicitly rather than by a variable so that there is no performance overhead.
	if bp.statistic == mean; statVec = [ convert(Float64, mean(sub(d, 1:N, m))) for m = 1:size(d, 2) ]
	elseif bp.statistic == median; statVec = [ convert(Float64, median(sub(d, 1:N, m))) for m = 1:size(d, 2) ]
	elseif bp.statistic == var; statVec = [ convert(Float64, var(sub(d, 1:N, m))) for m = 1:size(d, 2) ]
	elseif bp.statistic == std; statVec = [ convert(Float64, std(sub(d, 1:N, m))) for m = 1:size(d, 2) ]
	elseif bp.statistic == sum; statVec = [ convert(Float64, sum(sub(d, 1:N, m))) for m = 1:size(d, 2) ]
	elseif bp.statistic == quantile_001; statVec = [ convert(Float64, quantile(sub(d, 1:N, m), 0.001)) for m = 1:size(d, 2) ]
	elseif bp.statistic == quantile_01; statVec = [ convert(Float64, quantile(sub(d, 1:N, m), 0.01)) for m = 1:size(d, 2) ]
	elseif bp.statistic == quantile_05; statVec = [ convert(Float64, quantile(sub(d, 1:N, m), 0.05)) for m = 1:size(d, 2) ]
	elseif bp.statistic == quantile_1; statVec = [ convert(Float64, quantile(sub(d, 1:N, m), 0.1)) for m = 1:size(d, 2) ]
	elseif bp.statistic == quantile_9; statVec = [ convert(Float64, quantile(sub(d, 1:N, m), 0.9)) for m = 1:size(d, 2) ]
	elseif bp.statistic == quantile_95; statVec = [ convert(Float64, quantile(sub(d, 1:N, m), 0.95)) for m = 1:size(d, 2) ]
	elseif bp.statistic == quantile_99; statVec = [ convert(Float64, quantile(sub(d, 1:N, m), 0.99)) for m = 1:size(d, 2) ]
	elseif bp.statistic == quantile_999; statVec = [ convert(Float64, quantile(sub(d, 1:N, m), 0.999)) for m = 1:size(d, 2) ]
	else
		#If we reach this point, bp.statistic is not a recognized function and so we have to put up with the performance overhead
		statVec = [ convert(Float64, bp.statistic(sub(d, 1:N, m))) for m = 1:size(d, 2) ]
	end
	return(statVec)
end
#We define each of the above quantile functions locally. Note, these functions don't actually get called. I prefer to let the official quantile function handle sub-array case
quantile_001{T<:Number}(x::AbstractVector{T}) = quantile(x, 0.001)
quantile_01{T<:Number}(x::AbstractVector{T}) = quantile(x, 0.01)
quantile_05{T<:Number}(x::AbstractVector{T}) = quantile(x, 0.05)
quantile_1{T<:Number}(x::AbstractVector{T}) = quantile(x, 0.1)
quantile_9{T<:Number}(x::AbstractVector{T}) = quantile(x, 0.9)
quantile_95{T<:Number}(x::AbstractVector{T}) = quantile(x, 0.95)
quantile_99{T<:Number}(x::AbstractVector{T}) = quantile(x, 0.99)
quantile_999{T<:Number}(x::AbstractVector{T}) = quantile(x, 0.999)



#----------------------------------------------------------
#FUNCTION
#	dbootstrap
#	dbootstrap!
#PURPOSE
#	The purpose of this function is to return the bootstrapped distribution parameter of the statistic of interest
#NOTES
#	WARNING: This function is not type stable. This is because the user can input any function for transforming bootstrapped statistics into a distribution parameter.
#	Versions of this function for specific, popular parameters of interest that will be type-stable are below
#----------------------------------------------------------
#Function to return bootstrapped parameter. Function is not type stable, hence it farms most of the work out to other functions
dbootstrap!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = dbootstrap_getdistributionparam(dbootstrapstatistic!(x, bp))
dbootstrap!{T<:Number}(bp::BootstrapParam, x::AbstractVector{T}) = dbootstrap!(x, bp)
dbootstrap{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = dbootstrap_getdistributionparam(dbootstrapstatistic(x, bp))
dbootstrap{T<:Number}(bp::BootstrapParam, x::AbstractVector{T}) = dbootstrap(x, bp)
#Keyword argument wrapper
dbootstrap{T<:Number}(x::AbstractVector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean, distributionParam::Function=var) = dbootstrapstatistic(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic, distributionParam=distributionParam))
#Non-exported function used exclusively by dbootstrap to compute the distribution parameter (I've made this a separate function as I might make it more complicated in the future)
dbootstrap_getdistributionparam{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = bp.distributionParam(x)



#The following functions are dedicated type-stable versions of dbootstrap for specific distribution parameters.
dbootstrap_mean{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = mean(dbootstrapstatistic(x, bp))
dbootstrap_mean!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = mean(dbootstrapstatistic!(x, bp))
dbootstrap_mean{T<:Number}(x::AbstractVector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean) = dbootstrap_mean(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic))
dbootstrap_median{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = median(dbootstrapstatistic(x, bp))
dbootstrap_median!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = median(dbootstrapstatistic!(x, bp))
dbootstrap_median{T<:Number}(x::AbstractVector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean) = dbootstrap_median(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic))
dbootstrap_var{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = var(dbootstrapstatistic(x, bp))
dbootstrap_var!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = var(dbootstrapstatistic!(x, bp))
dbootstrap_var{T<:Number}(x::AbstractVector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean) = dbootstrap_var(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic))
dbootstrap_std{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = std(dbootstrapstatistic(x, bp))
dbootstrap_std!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam) = std(dbootstrapstatistic!(x, bp))
dbootstrap_std{T<:Number}(x::AbstractVector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean) = dbootstrap_std(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic))
dbootstrap_quantile{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam, p::Float64=0.95) = quantile(dbootstrapstatistic(x, bp), p)
dbootstrap_quantile!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam, p::Float64=0.95) = quantile(dbootstrapstatistic!(x, bp), p)
dbootstrap_quantile{T<:Number}(x::AbstractVector{T}, p::Float64=0.95; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean) = dbootstrap_quantile(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic), p)
dbootstrap_quantile{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam, p::Vector{Float64}=[0.025, 0.975]) = quantile(dbootstrapstatistic(x, bp), p)
dbootstrap_quantile!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam, p::Vector{Float64}=[0.025, 0.975]) = quantile(dbootstrapstatistic!(x, bp), p)
dbootstrap_quantile{T<:Number}(x::AbstractVector{T}, p::Vector{Float64}=[0.025, 0.975]; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean) = dbootstrap_quantile(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic), p)
dbootstrap_conf{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam, p::Vector{Float64}=[0.025, 0.975]) = length(p) != 2 ? error("Confidence interval must have exactly two probability inputs") : quantile(dbootstrapstatistic(x, bp), p)
dbootstrap_conf!{T<:Number}(x::AbstractVector{T}, bp::BootstrapParam, p::Vector{Float64}=[0.025, 0.975]) = length(p) != 2 ? error("Confidence interval must have exactly two probability inputs") : quantile(dbootstrapstatistic!(x, bp), p)
dbootstrap_conf{T<:Number}(x::AbstractVector{T}, p::Vector{Float64}=[0.025, 0.975]; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::BootstrapMethod=BootstrapStationary(), blockLengthMethod::BlockLengthMethod=BlockLength_Dummy(), blockLength::Number=-1, statistic::Function=mean) = dbootstrap_conf(x, BootstrapParam(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, blockLength=blockLength, statistic=statistic), p)




end # module
