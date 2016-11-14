
#Local function for checking a BootInput is valid for a call to dbootvecinds
function check_bi_for_dbootvecinds(bi::BootInput)
    bi.numobs < 1 && error("Number of observations must be strictly positive")
    bi.numresample < 1 && error("Number of resamples must be strictly positive")
    bi.numobsperresample < 1 && error("Number of observations per resample must be strictly positive")
    bi.blocklength <= 0.0 && error("Block length must be strictly positive")
    isnan(bi.blocklength) && error("Block length is set to NaN")
    isinf(bi.blocklength) && error("Block length is infinite")
    return(true)
end


"""
    dbootinds(bi::BootInput)::Vector{Vector{Int}}

Each inner vector of the returned Vector{Vector{Int}} provides indices that, when used to index the original dataset, will provide
a single resampled dataset. A keyword variant of this function is also provided.
"""
function dbootinds(bi::BootInput, bm::BootIID)::Vector{Vector{Int}}
    !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
    return(Vector{Int}[ rand(1:bi.numobs, bi.numobsperresample) for n = 1:bi.numresample ])
end
function dbootinds(bi::BootInput, bm::BootStationary)::Vector{Vector{Int}}
    !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
    bi.blocklength <= 1.0 && return(dbootinds(bi, BootIID()))
    inds = Vector{Int}[ Array(Int, bi.numobsperresample) for n = 1:bi.numresample ]
    geo1 = Geometric(1 / bi.blocklength)
    for n = 1:bi.numresample #Fill remainder using geometric draws
        (c, geoDraw) = (1, 1)
        for s = 1:bi.numobsperresample
            if c == geoDraw #Start a new block
                inds[n][s] = rand(1:bi.numobs)
                geoDraw = rand(geo1) + 1
                c = 1
            else #Next obs in existing block
                inds[n][s-1] == bi.numobs ? (inds[n][s] = 1) : (inds[n][s] = inds[n][s-1] + 1)
                c += 1
            end
        end
    end
	return(inds)
end
function dbootinds(bi::BootInput, bm::BootMoving)::Vector{Vector{Int}}
    !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
    bL = Int(ceil(bi.blocklength))
    bL == 1 && return(dbootinds(bi, BootIID()))
    inds = Vector{Int}[ Array(Int, bi.numobsperresample) for n = 1:bi.numresample ]
    blockStartUB = max(1, bi.numobs-bL+1)
    for n = 1:bi.numresample
        for s = 1:bL:bi.numobsperresample
            inds[n][s] = rand(1:blockStartUB) #Start of block
            for sInner = s+1:min(s+bL-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
                inds[n][sInner] = inds[n][sInner-1] + 1
            end
        end
    end
	return(inds)
end
function dbootinds(bi::BootInput, bm::BootNoOverlap)::Vector{Vector{Int}}
    !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
    bL = Int(ceil(bi.blocklength))
    bL == 1 && return(dbootinds(bi, BootIID()))
    inds = Vector{Int}[ Array(Int, bi.numobsperresample) for n = 1:bi.numresample ]
	blockStartUB = max(1, bi.numobs-bL+1)
	blockStartValues = collect(1:bL:blockStartUB) #Build valid set of start indices for any block
    (blockStartValues[end] + bL - 1 > bi.numobs) && (length(blockStartValues) > 1) && pop!(blockStartValues)
	for n = 1:bi.numresample
		for s = 1:bL:bi.numobsperresample
			inds[n][s] = blockStartValues[rand(1:length(blockStartValues))] #Start of block
			for sInner = s+1:min(s+bL-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
				inds[n][sInner] = inds[n][sInner-1] + 1
			end
		end
	end
	return(inds)
end
function dbootinds(bi::BootInput, bm::BootCircular)::Vector{Vector{Int}}
    !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
    bL = Int(ceil(bi.blocklength))
    bL == 1 && return(dbootinds(bi, BootIID()))
    inds = Vector{Int}[ Array(Int, bi.numobsperresample) for n = 1:bi.numresample ]
	for n = 1:bi.numresample
		for s = 1:bL:bi.numobsperresample
			inds[n][s] = rand(1:bi.numobs) #Start of block
			for sInner = s+1:min(s+bL-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
                inds[n][sInner-1] == bi.numobs ? (inds[n][sInner] = 1) : (inds[n][sInner] = inds[n][sInner-1] + 1)
			end
		end
	end
	return(inds)
end
dbootinds(bi::BootInput, bm::BootTapered)::Vector{Vector{Int}} = dbootinds(bi, BootMoving())
dbootinds(bi::BootInput)::Vector{Vector{Int}} = dbootinds(bi, bi.bootmethod)
dbootinds(data, bi::BootInput)::Vector{Vector{Int}} = dbootinds(bi) #This method is included for consistency
function dbootinds(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                    blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=numobs)::Vector{Vector{Int}}
    return(dbootinds(BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
end
