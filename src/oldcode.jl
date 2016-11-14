#CODE FOR MATRIX{INT} OUTPUT FROM DBOOTINDS
# """
#     dbootinds(bi::BootInput)::Matrix{Int}
#
# Each column of the returned Matrix{Int} provides indices that, when used to index the original dataset, will provide
# a single resampled dataset. A keyword variant of this function is also provided.
# """
# dbootinds(bi::BootInput, bm::BootIID)::Matrix{Int} = check_bi_for_dbootvecinds(bi) ? rand(1:bi.numobs, bi.numobsperresample, bi.numresample) : error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
# function dbootinds(bi::BootInput, bm::BootStationary)::Matrix{Int}
#     !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
#     bi.blocklength <= 1.0 && return(dbootinds(bi, BootIID()))
#     inds = Array(Int, bi.numobsperresample, bi.numresample)
#     for n in 1:bi.numresample #Fill first row of inds with discrete uniform draws
#         inds[1, n] = rand(1:bi.numobs)
#     end
#     geo1 = Geometric(1 / bi.blocklength)
#     for n = 1:bi.numresample #Fill remainder using geometric draws
#         geoDraw = rand(geo1) + 1
#         c = 1
#         for s = 2:bi.numobsperresample
#             if c == geoDraw #Start a new block
#                 inds[s, n] = rand(1:bi.numobs)
#                 geoDraw = rand(geo1) + 1
#                 c = 1
#             else #Next obs in existing block
#                 inds[s-1, n] == bi.numobs ? (inds[s, n] = 1) : (inds[s, n] = inds[s-1, n] + 1)
#                 c += 1
#             end
#         end
#     end
# 	return(inds)
# end
# function dbootinds(bi::BootInput, bm::BootMoving)::Matrix{Int}
#     !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
#     bL = Int(ceil(bi.blocklength))
#     bL == 1 && return(dbootinds(bi, BootIID()))
#     inds = Array(Int, bi.numobsperresample, bi.numresample) #Pre-allocate output
#     blockStartUB = max(1, bi.numobs-bL+1)
#     for n = 1:bi.numresample
#         for s = 1:bL:bi.numobsperresample
#             inds[s, n] = rand(1:blockStartUB) #Start of block
#             for sInner = s+1:min(s+bL-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
#                 inds[sInner, n] = inds[sInner-1, n] + 1
#             end
#         end
#     end
# 	return(inds)
# end
# function dbootinds(bi::BootInput, bm::BootNoOverlap)
#     !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
#     bL = Int(ceil(bi.blocklength))
#     bL == 1 && return(dbootinds(bi, BootIID()))
# 	inds = Array(Int, bi.numobsperresample, bi.numresample) #Pre-allocate output
# 	blockStartUB = max(1, bi.numobs-bL+1)
# 	blockStartValues = collect(1:bL:blockStartUB) #Build valid set of start indices for any block
#     (blockStartValues[end] + bL - 1 > bi.numobs) && (length(blockStartValues) > 1) && pop!(blockStartValues)
# 	for n = 1:bi.numresample
# 		for s = 1:bL:bi.numobsperresample
# 			inds[s, n] = blockStartValues[rand(1:length(blockStartValues))] #Start of block
# 			for sInner = s+1:min(s+bL-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
# 				inds[sInner, n] = inds[sInner-1, n] + 1
# 			end
# 		end
# 	end
# 	return(inds)
# end
# function dbootinds(bi::BootInput, bm::BootCircular)
#     !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
#     bL = Int(ceil(bi.blocklength))
#     bL == 1 && return(dbootinds(bi, BootIID()))
# 	inds = Array(Int, bi.numobsperresample, bi.numresample) #Pre-allocate output
# 	for n = 1:bi.numresample
# 		for s = 1:bL:bi.numobsperresample
# 			inds[s, n] = rand(1:bi.numobs) #Start of block
# 			for sInner = s+1:min(s+bL-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
#                 inds[sInner-1, n] == bi.numobs ? (inds[sInner, n] = 1) : (inds[sInner, n] = inds[sInner-1, n] + 1)
# 			end
# 		end
# 	end
# 	return(inds)
# end
# dbootinds(bi::BootInput, bm::BootTapered) = dbootinds(bi, BootMoving())
# dbootinds(bi::BootInput) = dbootinds(bi, bi.bootmethod)
# function dbootinds{T<:Number}(data::Array{T, N} ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
#                     blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=numobs)
#     return(dbootinds(BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample)))
# end
