
using 	DependentBootstrap,
		KernelStat

#Load any specific variables/functions that are needed (use import ModuleName1.FunctionName1, ModuleName2.FunctionName2, etc)

#Specify the variables/functions to export (use export FunctionName1, FunctionName2, etc)
# export 	testModuleBasic


#******************************************************************************

#----------------------------------------------------------
#SET CONSTANTS FOR MODULE
#----------------------------------------------------------
const bootstrapList = [BootstrapIID(), BootstrapStationary(), BootstrapMovingBlock(), BootstrapCircularBlock(), BootstrapNonoverlappingBlock(), BootstrapTaperedBlock()]::Vector{BootstrapMethod}
const blockLengthList = [BlockLengthPPW2009(), BlockLengthPP2002()]::Vector{BlockLengthMethod}
const kernelFuncTaperedList = [KernelPP2002Trap(), KernelPP2002Smooth()]::Vector{KernelFunction}
const bootstrapStringForPPWList = [:stationary, :circularBlock, :movingBlock]::Vector{Symbol}
const bandwidthList = [BandwidthWhiteNoise(), BandwidthBartlett(), BandwidthP2003()]::Vector{BandwidthMethod}
const statisticList = [mean, sum, var, std, median]::Vector{Function}
const distributionParamList = [mean, median, var, std, quantile]::Vector{Function}
#Functions for randomly generating input
randBL() = rand() < 0.8 ? -1 : rand(0:20)
randEBL() = rand() < 0.8 ? -1.0 : rand() * 20
randBootstrap(b::BootstrapIID) = BootstrapIID()
randBootstrap(b::BootstrapStationary) = BootstrapStationary(randEBL())
randBootstrap(b::BootstrapMovingBlock) = BootstrapMovingBlock(randBL())
randBootstrap(b::BootstrapCircularBlock) = BootstrapCircularBlock(randBL())
randBootstrap(b::BootstrapNonoverlappingBlock) = BootstrapNonoverlappingBlock(randBL())
randBootstrap(b::BootstrapTaperedBlock) = BootstrapTaperedBlock(randBL(), randKernel(kernelFuncTaperedList[rand(1:length(kernelFuncTaperedList))]))
randBandwidth(bW::BandwidthWhiteNoise) = BandwidthWhiteNoise()
randBandwidth(bW::BandwidthBartlett) = BandwidthBartlett()
randBandwidth(bW::BandwidthP2003) = BandwidthP2003(1.0, 2 + 0.15 * abs(randn()), convert(Int, round(5 + 0.3 * abs(randn()))))
randKernel(kT::KernelPP2002Trap) =  KernelPP2002Trap(0.5 * rand())
randKernel(kT::KernelPP2002Smooth) =  KernelPP2002Smooth(abs(randn()) + 1.0)
randBlockLength(b::BlockLengthPPW2009) = BlockLengthPPW2009(randBandwidth(bandwidthList[rand(1:length(bandwidthList))]), bootstrapStringForPPWList[rand(1:length(bootstrapStringForPPWList))])
randBlockLength(b::BlockLengthPP2002) = BlockLengthPP2002(randBandwidth(bandwidthList[rand(1:length(bandwidthList))]), randKernel(kernelFuncTaperedList[rand(1:length(kernelFuncTaperedList))]))
randStatistic() = statisticList[rand(1:length(statisticList))]
randDistributionParam() = distributionParamList[rand(1:length(distributionParamList))]
function randBootstrapParam(A::Int)
	if A == 1
		numObsData = rand(3:200)
		numObsResample = numObsData
		numResample = rand(1:20)
	elseif A == 2
		numObsData = 2000
		numObsResample = numObsData
		numResample = 600
	else
		error("Invalid input")
	end
	bootstrapMethod = randBootstrap(bootstrapList[rand(1:length(bootstrapList))])
	blockLengthMethod = randBlockLength(blockLengthList[rand(1:length(blockLengthList))])
	if typeof(bootstrapMethod) == BootstrapTaperedBlock
		rand() > 0.5 ? (statistic = mean) : (statistic = sum)
	else
		statistic = randStatistic()
	end
	distributionParam = randDistributionParam()
	return(BootstrapParam(numObsData, numObsResample, numResample, bootstrapMethod, blockLengthMethod, statistic, distributionParam))
end

#Function for testing module
function testModuleRandom(KK::Int)
	cNonIID = 0
	cPosBL = 0
	for kk = 1:KK
		#Test bootstrap param input
		b = randBootstrapParam(1)
		x = randn(b.numObsData)
		if typeof(b.bootstrapMethod) != BootstrapIID
			if getblocklength(b) <= 0
				cNonIID += 1
				blockLength = dbootstrapblocklength(x, b)
			end
		end
		bCopy = deepcopy(b)
		if typeof(bCopy.bootstrapMethod) != BootstrapIID
			if getblocklength(b) <= 0
				dbootstrapblocklength!(x, bCopy)
			end
		end
		if getblocklength(bCopy) > 0
			cPosBL += 1
			inds1 = dbootstrapindex(bCopy)
		end
		d1 = dbootstrapdata(x, b)
		d2 = dbootstrapdata!(x, b)
		s1 = dbootstrapstatistic(x, b)
		s2 = dbootstrapstatistic!(x, b)
		y1 = dbootstrap(x, b)
		y2 = dbootstrap!(x, b)
		#Test keyword arguments
		numObsData = b.numObsData
		numObsResample = b.numObsResample
		numResample = b.numResample
		bootstrapMethod = b.bootstrapMethod
		blockLengthMethod = b.blockLengthMethod
		bandwidthMethod = b.blockLengthMethod.bandwidthMethod
		statistic = b.statistic
		distributionParam = b.distributionParam
		if typeof(bootstrapMethod) != BootstrapIID
			bL = dbootstrapblocklength(x, blockLengthMethod=blockLengthMethod, bandwidthMethod=bandwidthMethod)
		end
		indsK = dbootstrapindex(x, numObsResample=numObsResample, numResample=numResample, blockLength=rand(1:5), bootstrapMethod=bootstrapMethod)
		d = dbootstrapdata(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod)
		d = dbootstrapdata(x, numObsResample=numObsResample, numResample=numResample, blockLength=rand(1:5), bootstrapMethod=bootstrapMethod)
		s = dbootstrapstatistic(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, statistic=statistic)
		s = dbootstrapstatistic(x, numObsResample=numObsResample, numResample=numResample, blockLength=rand(1:5), bootstrapMethod=bootstrapMethod, statistic=statistic)
		s = dbootstrap(x, numObsResample=numObsResample, numResample=numResample, bootstrapMethod=bootstrapMethod, blockLengthMethod=blockLengthMethod, statistic=statistic, distributionParam=distributionParam)
		s = dbootstrap(x, numObsResample=numObsResample, numResample=numResample, blockLength=rand(1:5), bootstrapMethod=bootstrapMethod, statistic=statistic, distributionParam=distributionParam)
	end
	println("Random test passed.")
	println("NonIID bootstrap method on " * string(cNonIID) * " of " * string(KK) * " iterations")
	println("Positive block length on " * string(cPosBL) * " of " * string(KK) * " iterations")
end


#Function for testing the block lengths
function testBlockLength()
	datasetDir = "/home/" * ENV["USER"] * "/.julia/v0.3/DependentBootstrap/test/"
	d1 = readcsv(datasetDir * "d1.csv")
	d2 = readcsv(datasetDir * "d2.csv")
	d3 = readcsv(datasetDir * "d3.csv")
	d4 = readcsv(datasetDir * "d4.csv")
	d5 = readcsv(datasetDir * "d5.csv")
	d6 = readcsv(datasetDir * "d6.csv")
	d7 = readcsv(datasetDir * "d7.csv")
	bMStationary = BlockLengthPPW2009(length(d1), :stationary)
	bMCircular = BlockLengthPPW2009(length(d1), :circularBlock)
	bMPPTrap = BlockLengthPP2002(length(d1), KernelPP2002Trap())
	bMPPSmooth = BlockLengthPP2002(length(d1), KernelPP2002Smooth())
	blPatton = [1.7616 2.0165; 6.0104 6.8802; 9.9051 11.3385; 11.2333 12.8589; 7.3506 8.4144; 11.8778 13.5967; 22.6527 25.9308]
	blOriginal = [1.756310026822603 2.010473102043058; 6.034031705365532 6.907242033150295; 9.957766996351472 11.398797704751043; 11.224677287339041 12.849047958881897; 7.9092103443486925 9.05378572852609; 11.896892646671194 13.61854245477252; 22.585333734911174 25.853753199173056]
	blMe = Array(Float64, 7, 2)
	blMePP = Array(Float64, 7, 2)
	for n = 1:7
		n == 1 && (x = vec(d1))
		n == 2 && (x = vec(d2))
		n == 3 && (x = vec(d3))
		n == 4 && (x = vec(d4))
		n == 5 && (x = vec(d5))
		n == 6 && (x = vec(d6))
		n == 7 && (x = vec(d7))
		blMe[n, 1] = dbootstrapblocklength(x, bMStationary)
		blMe[n, 2] = dbootstrapblocklength(x, bMCircular)
		blMePP[n, 1] = dbootstrapblocklength(x, bMPPTrap)
		blMePP[n, 2] = dbootstrapblocklength(x, bMPPSmooth)
	end
	for n = 1:7
		println("Dataset " * string(n) * ":")
		println("    Module block length (stationary)          = " * string(blMe[n, 1]))
		println("    Patton (Matlab) block length (stationary) = " * string(blPatton[n, 1]))
		println("    Module block length (circular)            = " * string(blMe[n, 2]))
		println("    Patton (Matlab) block length (circular)   = " * string(blPatton[n, 2]))
		println("NOTE: Small deviations between this module and Pattons are expected due to minor differences in method for computing autocovariances")
	end
	if sum(sum(abs(blMe - blOriginal))) > 1e-10
		error("Block length estimates are significantly different from those originally created by this module. A recent update may have significantly altered the outcomes. Please review.")
	end
	for n = 1:7
		println("Dataset " * string(n) * ":")
		println("    Module block length (stationary)  = " * string(blMe[n, 1]))
		println("    Module block length (circular)    = " * string(blMe[n, 2]))
		println("    Module block length (taperedTrap) = " * string(blMePP[n, 1]))
		println("    Module block length (taperedSmth) = " * string(blMePP[n, 1]))
	end
end


function testBootstrapIndexVisual()
	numObs = 20
	blockLength = 5
	numResample = 5
	bpIID = BootstrapParam(numObs, bootstrapMethod=BootstrapIID())
	bpSB = BootstrapParam(numObs, bootstrapMethod=BootstrapStationary())
	bpCB = BootstrapParam(numObs, bootstrapMethod=BootstrapCircularBlock())
	bpMB = BootstrapParam(numObs, bootstrapMethod=BootstrapMovingBlock())
	bpN = BootstrapParam(numObs, bootstrapMethod=BootstrapNonoverlappingBlock(), blockLength=rand(2:8))
	#First test bootstrap indices on small datasets where output can be inspected
	for k = 1:5
		if k == 1
			bp = bpIID; println("Bootstrap method = iid. Number of observations = " * string(numObs))
		end
		if k == 2
			bp = bpSB; println("Bootstrap method = stationary (expected block length of " * string(blockLength) * "). Number of observations = " * string(numObs))
		end
		if k == 3
			bp = bpCB; println("Bootstrap method = circularBlock (block length of " * string(blockLength) * "). Number of observations = " * string(numObs))
		end
		if k == 4
			bp = bpMB; println("Bootstrap method = movingBlock (block length of " * string(blockLength) * "). Number of observations = " * string(numObs))
		end
		if k == 5
			bp = bpN; println("Bootstrap method = nonoverlappingBlock (block length of " * string(blockLength) * "). Number of observations = " * string(numObs))
		end
		if k > 1
			update!(bp, blockLength=blockLength)
		end
		inds = dbootstrapindex(bp.bootstrapMethod, numObs, numObs, numResample)
		for m = 1:numResample
			println("    Resampling indices number " * string(m) * ":" * string(vec(inds[:, m])))
		end
	end
end


#Function for testing the bootstrap
function testBootstrapVisual()
	numResample = 10
	datasetDir = "/home/" * ENV["USER"] * "/.julia/v0.3/DependentBootstrap/test/"
	d1 = readcsv(datasetDir * "d1.csv")
	d2 = readcsv(datasetDir * "d2.csv")
	d3 = readcsv(datasetDir * "d3.csv")
	d4 = readcsv(datasetDir * "d4.csv")
	d5 = readcsv(datasetDir * "d5.csv")
	d6 = readcsv(datasetDir * "d6.csv")
	d7 = readcsv(datasetDir * "d7.csv")
	d1M5 = readcsv(datasetDir * "d1_Mean5.csv")
	d2M5 = readcsv(datasetDir * "d2_Mean5.csv")
	d3M5 = readcsv(datasetDir * "d3_Mean5.csv")
	d4M5 = readcsv(datasetDir * "d4_Mean5.csv")
	d5M5 = readcsv(datasetDir * "d5_Mean5.csv")
	d6M5 = readcsv(datasetDir * "d6_Mean5.csv")
	d7M5 = readcsv(datasetDir * "d7_Mean5.csv")
	bpIID = BootstrapParam(length(d1), bootstrapMethod=BootstrapIID())
	bpSB = BootstrapParam(length(d1), bootstrapMethod=BootstrapStationary())
	bpCB = BootstrapParam(length(d1), bootstrapMethod=BootstrapCircularBlock())
	bpMB = BootstrapParam(length(d1), bootstrapMethod=BootstrapMovingBlock())
	bpT = BootstrapParam(length(d1), bootstrapMethod=BootstrapTaperedBlock())
	for n = 1:7
		println("Dataset " * string(n))
		n == 1 && (x = vec(d1))
		n == 2 && (x = vec(d2))
		n == 3 && (x = vec(d3))
		n == 4 && (x = vec(d4))
		n == 5 && (x = vec(d5))
		n == 6 && (x = vec(d6))
		n == 7 && (x = vec(d7))
		n == 1 && (x5 = vec(d1M5))
		n == 2 && (x5 = vec(d2M5))
		n == 3 && (x5 = vec(d3M5))
		n == 4 && (x5 = vec(d4M5))
		n == 5 && (x5 = vec(d5M5))
		n == 6 && (x5 = vec(d6M5))
		n == 7 && (x5 = vec(d7M5))
		for k = 1:5
			k == 1 && (bp = bpIID)
			k == 2 && (bp = bpSB)
			k == 3 && (bp = bpCB)
			k == 4 && (bp = bpMB)
			k == 5 && (bp = bpT)
			k == 1 && println("Bootstrap method = iid")
			k == 2 && println("Bootstrap method = stationary")
			k == 3 && println("Bootstrap method = circularBlock")
			k == 4 && println("Bootstrap method = movingBlock")
			k == 5 && println("Bootstrap method = taperedBlock")
			k == 5 && (x = x - mean(x))
			k == 5 && (x5 = x5 - mean(x5))
			update!(bp, numResample=numResample)
			meanStat = dbootstrapstatistic!(bp, x)
			println("    Mean = 0, mean stats:")
			if k == 5
				println("    NOTE: Data de-meaned before applying tapered block bootstrap")
			end
			for j = 1:length(meanStat)
				@printf(STDOUT, "    %1.3f", meanStat[j])
			end
			meanStat = dbootstrapstatistic!(bp, x5)
			println("")
			println("    Mean = 5, mean stats:")
			if k == 5
				println("    NOTE: Data de-meaned before applying tapered block bootstrap")
			end
			for j = 1:length(meanStat)
				@printf(STDOUT, "    %1.3f", meanStat[j])
			end
			println("")
		end
	end
end




testModuleRandom(100)
testBlockLength()
testBootstrapIndexVisual()
testBootstrapVisual()
println("All tests passed.")





