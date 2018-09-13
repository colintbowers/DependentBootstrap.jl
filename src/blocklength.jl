
"""
    optblocklength(data, bi::BootInput)::Float64
    optblocklength(data ; kwargs...)::Float64

Provides an estimate of the optimal block-length to use with a dependent bootstrap.

For multivariate datasets, optimal block length is estimated for each column of data, and then
bi.fblocklengthcombine, which is a function that maps Vector{Float64} to Float64, is called
to reduce the multiple estimates to a single estimates. The default value for fblocklengthcombine
is median.

Block length methods currently implemented include: \n
     - Patton, Politis, White (2009) "Correction to Automatic Block Length Selection For the Dependent Bootstrap" \n

For all methods discussed above, bandwidth is estimated following Politis (2003) "Adaptive Bandwidth Choice", using the
flat-top kernel suggested in that paper.
"""
function optblocklength(x::AbstractVector{<:Number}, blm::BLPPW2009{P2003}, bootmethod::Tbm)::Float64 where {Tbm<:BootMethod}
	length(x) < 3 && error("You must have at least 3 observations to estimate block length")
    (M, xVar, covVec) = blocklength_ma_and_cor(x) #Bandwidth method currently forced to politis (2003)
	kernelCovVec = blocklength_kernel_cov(covVec, M)
	gHat = 0.0
	for k = 1:M
		gHat += 2 * k * kernelCovVec[k] #note, "2*" is because sum is from -M and M, but is symmetric about 0. Note, 0 term in sum vanishes since k=0 -> |k|=0
	end
    dHat = optblocklength_blppw2009_dhat(bootmethod, xVar, kernelCovVec)
    #Equation 9 and 14 from Politis, White (2004)
	blocklength = (2 * gHat^2 * length(x) / dHat)^(1/3)
	blocklength = min(blocklength, ceil(min(3*sqrt(length(x)), length(x) / 3))) #Enforce upper bound on block length suggested by Patton
	blocklength = max(blocklength, 1.0)
	return blocklength
end
(optblocklength_blppw2009_dhat(a::BootStationary, xvar::Float64, kernelcovvec::Vector{Float64})::Float64) = 2 * (xvar + 2*sum(kernelcovvec))^2 #note, in expression (1 + 2*sum(kernelCovVec)), "2*" is because sum is from -M to M, but is symmetric around 0. "1+" is the zero term of the sum which is simply equal to unity.
(optblocklength_blppw2009_dhat(a::T, xvar::Float64, kernelcovvec::Vector{Float64})::Float64) where {T<:Union{BootCircular,BootMoving,BootIID,BootNoOverlap}} = (4/3) * (xvar + 2*sum(kernelcovvec))^2 #note, in expression (1 + 2*sum(kernelCovVec)), "2*" is because sum is from -M to M, but is symmetric around 0. "1+" is the zero term of the sum which is simply equal to unity.



#Block length selection method of Paparoditis, Politis (2002) "The Tapered Block Bootstrap for General Statistics From Stationary Sequences"
# function optblocklength(x::AbstractVector{<:Number}, blm::BLPP2002{P2003,Tkf}, bootmethod::Tbm)::Float64 where {Tbm<:BootMethod,Tkf<:KernelFunctionMethod}
# 	length(x) < 3 && error("You must have at least 3 observations to estimate block length")
#     Tbm != BootTapered && println("WARNING: Optimal parameter values in pp2002 block length procedure are unknown for bootstrap method $(bootmethod).")
# 	(M, xVar, covVec) = blocklength_ma_and_cor(x)
# 	kernelCovVec = blocklength_kernel_cov(covVec, M)
# 	deltaUnknown = (xVar + 2*sum(kernelCovVec))^2 #Unknown parameter in Delta expression, section 4, Paparoditis, Politis (2002). (note, "1+" is for k=0 term in summation, "2*" is because summation is from -M to M, but is symmetric about 0)
# 	gammaUnknown = 0.0 #Unknown parameter in Gamma expression, start of section 4
# 	for k = 1:M
# 		gammaUnknown += 2 * k^2 * kernelCovVec[k] #note, "2*" is because sum is from -M and M, but is symmetric about 0. Note, 0 term in sum vanishes since k=0 -> k^2=0
# 	end
#     (gammaHat, deltaHat) = optblocklength_blpp2002_param(blm.kernelfunction, gammaUnknown, deltaUnknown)
#     #Equation 20 from Paraproditis, Politis (2002)
# 	blocklength = (4 * gammaHat^2 * length(x) / deltaHat)^(1/5)
# 	blocklength = min(blocklength, ceil(min(3*sqrt(length(x)), length(x) / 3))) #Enforce upper bound on block length suggested by Patton
# 	blocklength = max(blocklength, 1.0)
# 	return blocklength
# end
# (optblocklength_blpp2002_param(kf::KernelTrap, gu::Float64, du::Float64)::Float64) = (-5.45*gu, 1.099*du) #-5.45 = (1/2) * -10.9 = (1/2) * (w*w)''(0) / (w*w)(0) [OPTIMAL TRAP VALUES], see Paparoditis, Politis (2002), 1.099 = 2 * 0.5495 = 2 * int_{-1}^{1} ((w*w)^2(x) / (w*w)^2(0)) dx [OPTIMAL TRAP VALUES], see Paparoditis, Politis (2002)
# (optblocklength_blpp2002_param(kf::KernelSmooth, gu::Float64, du::Float64)::Float64) = (-5.175*gu, 1.1312*deltaUnknown) #-5.175 = (1/2) * -10.35 = (1/2) * (w*w)''(0) / (w*w)(0) [OPTIMAL SMOOTH VALUES], see Paparoditis, Politis (2002), 1.1312 = 2 * 0.5656 = 2 * int_{-1}^{1} ((w*w)^2(x) / (w*w)^2(0)) dx [OPTIMAL SMOOTH VALUES], see Paparoditis, Politis (2002)
# (optblocklength_blpp2002_param(kf::T, gu::Float64, du::Float64)::Float64) where {T} = error("Optimal parameters not known for input kernel function: $(kf)")
(optblocklength(x::AbstractVector{<:Number}, bi::BootInput{Tbm,BLDummy})::Float64) where {Tbm<:BootMethod} = error("Logic fail. It should not have been possible to call this method. Please file an issue with full stacktrace.")
(optblocklength(x::AbstractVector{<:Number}, bi::BootInput)::Float64) = optblocklength(x, bi.blocklengthmethod, bi.bootmethod)
#Multivariate dataset wrapper
(optblocklength(data, bl::Tbl, bm::Tbm, f::Tf)::Float64) where {Tbl<:BlockLengthMethod,Tbm<:BootMethod,Tf<:Function} = f([ optblocklength(local_get_var(data, i), bl, bm) for i = 1:num_var(data) ])
(optblocklength(data, bi::BootInput)::Float64) = bi.fblocklengthcombine([ optblocklength(local_get_var(data, i), bi) for i = 1:num_var(data) ])
#Keyword input wrapper
(optblocklength(data ; kwargs...)::Float64) = optblocklength(data, BootInput(data, kwargs...))

#These two functions are used by several of the block-length selection procedures
function blocklength_ma_and_cor(x::AbstractVector{T})::Tuple{Int, Float64, Vector{Float64}} where {T<:Number}
	(M, xVar, covVec) = bandwidth_politis_2003(x)
	if M > 0.5 * length(x)
		println("WARNING: Bandwidth in parameter estimation section of blocklength forced to half total number of observations. Data may contain excessive dependence.")
		M = Int(round(0.5 * length(x)))
	end
	M < 2 && (M = 2) #Even though M output of bandwidth function is always greater than 2, the above check for excessively large M can result in M being reduced below 2 (admittedly only in very unusual circumstances)
	length(covVec) < M && append!(covVec, autocov(x, length(covVec)+1:M)) #Get any additional autocovariances that we might need
	return(M, xVar, covVec)
end
function blocklength_kernel_cov(covVec::AbstractVector{Float64}, M::Int)::Vector{Float64}
	length(covVec) < M && error("Error in blocklength_kernel_cov likely caused by logic fail in blocklength_ma_and_cor. Please lodge a github issue, preferably with reproducible example.")
	kernelCovVec = Float64[ kernel_politis_2003_flat_top(k / M) * covVec[k] for k = 1:M ]
	return(kernelCovVec)
end

#Note, the following non-exported function is called in the ForecastEval package, so if you alter the function name you will
#need to adjust that package too.
"""
    bandwidth_politis_2003(x::AbstractVector{T})::Tuple{Int, Float64, Vector{Float64}} where {T<:Number}

Implements the methodology from Politis (2003) "Adaptive Bandwidth Choice" to obtain a data-driven bandwidth estimate.

Return tuple is, in order, the bandwidth estimate, the variance of x, and the autocorrelations used to get the bandwidth estimate.

Note, most users won't be interested in the second and third output, but sometimes this routine will be called by other
functions that need these terms, so they are returned to avoid duplicate computation.
"""
function bandwidth_politis_2003(x::AbstractVector{T})::Tuple{Int, Float64, Vector{Float64}} where {T<:Number}
    length(x) < 2 && error("Input data must have at least two observations")
    adjustmentTerm = 1.0 #This term is used by me for debugging. It serves no statistical purpose.
    politis_c = 2.0 #This value is again recommended in Politis (2003)
    K = max(5, Int(ceil(sqrt(log(10, length(x)))))) #This value is recommended in Politis (2003)
	mHat = 1
	corVec = Float64[]
	append!(corVec, autocor(x, 1:min(20, length(x)-1))) #I add autocorrelations in blocks as this is more efficient than doing it one at a time
	corBound = politis_c * sqrt(log(10, length(x)) / length(x)) #Note, use of log base 10 is deliberate and recommended in Politis (2003)
	KCounter = 0
	mHatFound = false
    for k = 1:length(x)-2
        #Add one to counter if bound is satisfied, otherwise reset counter to 0
        abs(corVec[k]) < corBound ? (KCounter += 1) : (KCounter = 0)
		if KCounter >= K #We found K autocorrelations in a row that satisfy the bound, so break.
			mHat = k - K + 1
			mHatFound = true
			break
		end
        #If we run out of autocorrelations to check, add another block of them to corVec
        k == length(corVec) && append!(corVec, autocor(x, length(corVec)+1:min(length(corVec)+20, length(x)-1)))
	end
    mHatFound == false && (mHat = length(x) - 1) #Bound mHat in the case where we fail to hit the break condition
	M = Int(ceil(adjustmentTerm * 2 * mHat)) #"2*" is a standard rule, see e.g. Politis (2003).
    M > length(x) - 1 && (M = length(x) - 1) #Apply upper bound to M
    M < 2 && (M = 2) #Apply lower bound to M
	xVar = ((length(x)-1) /  length(x)) * var(x) #Used to scale autocorrelations to autocovariances for return argument
	return (M, xVar, xVar * corVec)
end

"""
    kernel_politis_2003_flat_top(x::Float64)::Float64

Implements the flat-top kernel function discussed in Politis (2003) "Adaptive Bandwidth Choice"
"""
function kernel_politis_2003_flat_top(x::Float64)::Float64
    x_abs = abs(x)
    x_abs <= 0.5 && return 1.0
    x_abs <= 1.0 && return (2.0 * (1 - x_abs))
    return 0.0
end
