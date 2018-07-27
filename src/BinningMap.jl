import Base: eof, close
export BinningMap, close, value, position, eof, advance!

type BinningMap
    reader::BamReader
    binSize::Int64
    position::Int64
    value::Float64
    skipDup::Bool
end

function BinningMap(reader::BamReader, binSize; skipDup=true)
    fm = BinningMap(reader, binSize, 0, 0.0, skipDup)
    
    advance!(fm)
    fm
end
close(fm::BinningMap) = close(fm.reader)
value(fm::BinningMap) = fm.value
position(fm::BinningMap) = fm.position
eof(fm::BinningMap) = fm.position <= 0

function advance!(fm::BinningMap)
    fm.position = floor((fm.reader.position-1)/fm.binSize) + 1
    binEnd = fm.position*fm.binSize
    
    # Fill in the bin
    fm.value = 0.0
    lastPos = -1
    while fm.reader.position != -1 && fm.reader.position <= binEnd
        if !fm.skipDup || fm.reader.position != lastPos
            fm.value += 1
        end
        lastPos = fm.reader.position
        advance!(fm.reader)
    end
end