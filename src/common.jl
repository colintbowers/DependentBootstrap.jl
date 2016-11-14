
#Core type methods
Base.show(io::IO, bm::BootMethod) = string(shortsymbol(bm))
Base.show(io::IO, blm::BLMethod) = string(shortsymbol(blm))
function show(io::IO, b::BootInput)
	println(io, "dependent bootstrap input:")
	println(io, "    number of observations in dataset = " * string(b.numobs))
    println(io, "    current block length = " * string(b.blocklength))
	println(io, "    number of resamples = " * string(b.numresample))
	println(io, "    bootstrap method = " * string(b.bootmethod))
	println(io, "    block length method = " * string(b.blmethod))
	println(io, "    level 1 statistic = " * string(b.flevel1))
	println(io, "    level 2 statistic = " * string(b.flevel2))
end

"""
    setblocklength!(b::BootInput, blocklength::Number)

Adjust the block length specified in a BootInput
"""
setblocklength!(b::BootInput, blocklength::Number) = (b.blocklength = Float64(blocklength))

"""
    setnumresample!(b::BootInput, numresample::Number)

Adjust the number of resamples specified in a BootInput
"""
setnumresample!(b::BootInput, numresample::Number) = (b.numresample = Int(numresample))

"""
    setflevel1!(b::BootInput, flevel1::Function)

Adjust the flevel1 function in a BootInput
"""
setflevel1!(b::BootInput, flevel1::Function) = (b.flevel1 = flevel1)

"""
    setflevel2!(b::BootInput, flevel1::Function)

Adjust the flevel2 function in a BootInput
"""
setflevel2!(b::BootInput, flevel2::Function) = (b.flevel2 = flevel2)
