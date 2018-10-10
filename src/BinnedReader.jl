import Base: eof, close
export BinnedReader, close, position, value, eof, advance!, write_binned

mutable struct BinnedReader
    fileStream
    pair::Array{UInt32}
end

function BinnedReader(fileName::String)
    f = open(fileName)
    br = BinnedReader(f, zeros(UInt32, 2))
    advance!(br)
    br
end
close(br::BinnedReader) = close(br.fileStream)
value(br::BinnedReader) = br.pair[2]
position(br::BinnedReader) = br.pair[1]
eof(br::BinnedReader) = br.pair[1] == 0

function advance!(br::BinnedReader)
    if !eof(br.fileStream)
        read!(br.fileStream, br.pair)
    else
        br.pair[1] = 0 # mark that we are at eof
    end
end

function write_binned(bamFile::String, binSize::Int64, readOrientation; skipDup=true)
    bm = BinningMap(BamReader(bamFile, readOrientation, ReferenceContigs_hg38), binSize, skipDup=skipDup)
    out = open(bamFile*"."*(readOrientation == :reverse ? "r" : (readOrientation == :forward ? "f" : "a"))*"bin$binSize", "w")
    while !eof(bm)
        write(out, UInt32(bm.position))
        write(out, UInt32(bm.value))
        advance!(bm)
    end
    close(out)
end


function write_binned(bamFile::String, target::String, binSize::Int64, readOrientation; skipDup=true)
    bm = BinningMap(BamReader(bamFile, readOrientation, ReferenceContigs_hg38), binSize, skipDup=skipDup)
    out = open(target, "w")
    while !eof(bm)
        write(out, UInt32(bm.position))
        write(out, UInt32(bm.value))
        advance!(bm)
    end
    close(out)
end
