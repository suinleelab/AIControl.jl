export denseblocks

mutable struct DenseBlockIterator
    readers::Array{Any}
    blockSize::Int64
    blockWidth::Int64
    block::Array{Float64,2}
    offset::Int64
    done::Bool
    constantColumn::Bool
    loop::Bool
end

function denseblocks(readers, blockSize::Int64; constantColumn=false, loop=false)
    blockWidth = constantColumn ? length(readers) + 1 : length(readers)
    block = ones(Float64, blockSize, blockWidth)
    if constantColumn
        block[:,end] .= 1.0
    end
    DenseBlockIterator(readers, blockSize, blockWidth, block, 0, false, constantColumn, loop)
end

#Depricated for Julia 1.0
#Base.start(it::DenseBlockIterator) = 0
#Base.done(it::DenseBlockIterator, nil) = it.done && !it.loop # never done with a constant column

#function Base.next(it::DenseBlockIterator, nil)
function next(it::DenseBlockIterator, nil)
    
    if it.constantColumn
        it.block[:,1:end-1] .= 0.0
    else
        it.block[:,:] .= 0.0
    end

    # Fill in the block
    if !it.done

        foundRead = false
        for i in 1:length(it.readers)
            reader = it.readers[i]

            while !eof(reader) && position(reader) <= it.offset + it.blockSize
                it.block[position(reader) - it.offset, i] += value(reader)
                advance!(reader)
                foundRead = true
            end
        end

        # See if we are really done or just found a blank block
        if !foundRead
            it.done = true
            for i in 1:length(it.readers)
                it.done = it.done && eof(it.readers[i])
            end
        end
    end
    
    # update the offset
    it.offset += it.blockSize

    it.block, 0
end

function Base.iterate(it::DenseBlockIterator, state=0)
    if it.done && !it.loop
        return nothing
    else
        return next(it, state)
    end
end
