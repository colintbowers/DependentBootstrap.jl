
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
const bootstrapStringForPPWList = ["stationary", "circularBlock", "movingBlock"]::Vector{ASCIIString}
const bandwidthList = [BandwidthWhiteNoise(), BandwidthBartlett(), BandwidthP2003()]::Vector{BandwidthMethod}
const statisticList = ["mean", mean, "sum", sum, "variance", var, "std", std, "median", median, "quantile"]::Vector{Any}
const distributionParamList = ["mean", mean, "median", median, "variance", var, "std", std, "quantile", "conf"]::Vector{Any}
#Functions for randomly generating input
function randBL()
	if rand() < 0.8 #Most of the time we force the function to use routines that estimate block length
		return(-1)
	else
		return(rand(0:20))
	end
end
function randEBL()
	if rand() < 0.8 #Most of the time we force the function to use routines that estimate block length
		return(-1.0)
	else
		return(rand() * 20)
	end
end
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
function randStatistic()
	x = statisticList[rand(1:length(statisticList))]
	if typeof(x) == ASCIIString
		if x == "quantile"
			n = rand(1:99)
			if n < 10
				nStr = "0" * string(n)
			else
				nStr = string(n)
			end
			x = "quantile" * nStr
		end
	end
	return(x)
end
function randDistributionParam()
	x = distributionParamList[rand(1:length(distributionParamList))]
	if typeof(x) == ASCIIString
		if x == "quantile"
			n = rand(1:99)
			if n < 10
				nStr = "0" * string(n)
			else
				nStr = string(n)
			end
			x = "quantile" * nStr
		end
	end
	return(x)
end
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
		if rand() > 0.5
			statistic = mean
		else
			statistic = "mean"
		end
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
		if typeof(b.bootstrapMethod) == BootstrapTaperedBlock
			x = x - mean(x)
		end
		if typeof(b.bootstrapMethod) != BootstrapIID
			cNonIID += 1
			blockLength = dBootstrapBlockLength(x, b)
		end
		bCopy = copy(b)
		if typeof(bCopy.bootstrapMethod) != BootstrapIID
			dBootstrapBlockLength!(x, bCopy)
		end
		if getBlockLength(bCopy) > 0
			cPosBL += 1
			inds1 = dBootstrapIndex(bCopy)
			inds2 = dBootstrapIndex(x, bCopy)
		end
		d1 = dBootstrapData(x, b)
		d2 = dBootstrapData!(x, b)
		s1 = dBootstrapStatistic(x, b)
		s2 = dBootstrapStatistic!(x, b)
		y1 = dBootstrap(x, b)
		y2 = dBootstrap!(x, b)
		#Test keyword arguments
		numObsDataStr = b.numObsData
		numObsResampleStr = b.numObsResample
		numResampleStr = b.numResample
		bootstrapMethodStr = string(b.bootstrapMethod)
		blockLengthMethodStr = string(b.blockLengthMethod)
		bandwidthMethodStr = string(b.blockLengthMethod.bandwidthMethod)
		blockLengthIndexFuncStr = rand(-8:3)
		statisticStr = string(b.statistic)
		distributionParamStr = string(b.distributionParam)
		if bootstrapMethodStr != "iid"
			bL = dBootstrapBlockLength(x, blockLengthMethod=blockLengthMethodStr, bootstrapMethod=bootstrapMethodStr, bandwidthMethod=bandwidthMethodStr)
		end
		if blockLengthIndexFuncStr > 0
			indsK = dBootstrapIndex(x, numObsResample=numObsResampleStr, numResample=numResampleStr, blockLength=blockLengthIndexFuncStr, bootstrapMethod=bootstrapMethodStr, blockLengthMethod=blockLengthMethodStr, bandwidthMethod = bandwidthMethodStr)
		end
		d = dBootstrapData(x, numObsResample=numObsResampleStr, numResample=numResampleStr, blockLength=blockLengthIndexFuncStr, bootstrapMethod=bootstrapMethodStr, blockLengthMethod=blockLengthMethodStr, bandwidthMethod = bandwidthMethodStr)
		s = dBootstrapStatistic(x, numObsResample=numObsResampleStr, numResample=numResampleStr, blockLength=blockLengthIndexFuncStr, bootstrapMethod=bootstrapMethodStr, blockLengthMethod=blockLengthMethodStr, bandwidthMethod = bandwidthMethodStr, statistic=statisticStr)
		s = dBootstrap(x, numObsResample=numObsResampleStr, numResample=numResampleStr, blockLength=blockLengthIndexFuncStr, bootstrapMethod=bootstrapMethodStr, blockLengthMethod=blockLengthMethodStr, bandwidthMethod = bandwidthMethodStr, statistic=statisticStr, distributionParam=distributionParamStr)
	end
	println("Random test passed.")
	println("NonIID bootstrap method on " * string(cNonIID) * " of " * string(KK) * " iterations")
	println("Positive block length on " * string(cPosBL) * " of " * string(KK) * " iterations")
end


#Function for testing the block lengths
function testBlockLength()
	datasetDir = "/home/colin/Cloud/Codebase/Julia/Modules/ctbDependentBootstrap/test/"
	d1 = readcsv(datasetDir * "d1.csv")
	d2 = readcsv(datasetDir * "d2.csv")
	d3 = readcsv(datasetDir * "d3.csv")
	d4 = readcsv(datasetDir * "d4.csv")
	d5 = readcsv(datasetDir * "d5.csv")
	d6 = readcsv(datasetDir * "d6.csv")
	d7 = readcsv(datasetDir * "d7.csv")
	bMStationary = BlockLengthPPW2009("stationary", length(d1))
	bMCircular = BlockLengthPPW2009("circularBlock", length(d1))
	bMPPTrap = BlockLengthPP2002(KernelPP2002Trap(), length(d1))
	bMPPSmooth = BlockLengthPP2002(KernelPP2002Smooth(), length(d1))
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
		blMe[n, 1] = dBootstrapBlockLength(x, bMStationary)
		blMe[n, 2] = dBootstrapBlockLength(x, bMCircular)
		blMePP[n, 1] = dBootstrapBlockLength(x, bMPPTrap)
		blMePP[n, 2] = dBootstrapBlockLength(x, bMPPSmooth)
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
	bpIID = BootstrapParam(numObs, bootstrapMethod="iid")
	bpSB = BootstrapParam(numObs, bootstrapMethod="stationary")
	bpCB = BootstrapParam(numObs, bootstrapMethod="circularBlock")
	bpMB = BootstrapParam(numObs, bootstrapMethod="movingBlock")
	bpN = BootstrapParam(numObs, bootstrapMethod = "nonoverlappingBlock")
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
			replaceBlockLength!(bp, blockLength)
		end
		inds = dBootstrapIndex(bp.bootstrapMethod, numObs, numObs, numResample)
		for m = 1:numResample
			println("    Resampling indices number " * string(m) * ":" * string(vec(inds[:, m])))
		end
	end
end



#Function for testing the bootstrap
function testBootstrapVisual()
	numResample = 10
	datasetDir = "/home/colin/Cloud/Codebase/Julia/Modules/ctbDependentBootstrap/test/"
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
	bpIID = BootstrapParam(length(d1), bootstrapMethod="iid")
	bpSB = BootstrapParam(length(d1), bootstrapMethod="stationary")
	bpCB = BootstrapParam(length(d1), bootstrapMethod="circularBlock")
	bpMB = BootstrapParam(length(d1), bootstrapMethod="movingBlock")
	bpT = BootstrapParam(length(d1), bootstrapMethod="taperedBlock")
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
			replaceNumResample!(bp, numResample)
			meanStat = dBootstrapStatistic!(bp, x)
			println("    Mean = 0, mean stats:")
			for j = 1:length(meanStat)
				@printf(STDOUT, "    %1.3f", meanStat[j])
			end
			meanStat = dBootstrapStatistic!(bp, x5)
			println("")
			println("    Mean = 5, mean stats:")
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





