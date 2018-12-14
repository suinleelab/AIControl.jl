using GZip

import Base: eof, close, position
export BamReader, close, value, eof, advance!, eachposition

mutable struct BamReader
    bamStream
    readOrientation #useReverseReads::Bool
    done::Bool
    position::Int64
    contigs::ReferenceContigs
end

function BamReader(bamFileName::String, readOrientation, contigs)
    f = GZip.open(bamFileName)

    # make sure this is a BAM file
    code = read!(f, Array{UInt8}(undef, 4))
    @assert code == b"BAM\1"

    # get through the header data
    l_text = read(f, Int32)
    skip(f, l_text)

    # make sure the contigs match our reference
    n_ref = read(f, Int32)
    if !(n_ref == contigs.count)
        println("Your bam files is not aligned to the UCSC hg38 genome.")
        println("See the step 3.1 at https://github.com/hiranumn/AIControl.jl to realign your genome")
        println("to the specific version of hg38 using bowtie2.")
        exit()
    end
        
    for j in 1:n_ref
        l_name = read(f, Int32)
        refName = String(read(f, Array{UInt8}(undef, l_name))[1:end-1]) # ignore the null terminator
        l_ref = read(f, Int32)
        if !(l_ref == contigs.sizes[j]) || !(refName == contigs.names[j])
            println("Your bam files is not aligned to the UCSC hg38 genome.")
            println("See the step 3.1 at https://github.com/hiranumn/AIControl.jl to realign your genome")
            println("to the specific version of hg38 using bowtie2.")
            exit()
        end
    end

    r = BamReader(f, readOrientation, false, 1, contigs)
    advance!(r)
    r
end
close(reader::BamReader) = GZip.close(reader.bamStream)
value(reader::BamReader) = 1
position(reader::BamReader) = reader.position
eof(reader::BamReader) = reader.position == -1

function advance!(r::BamReader)
    f = r.bamStream
    while !r.done
        if peek(f) == -1 # eof does not work on the BAM files either in C++ or here (BAM vs. gzip issue?)
            r.done = true
            r.position = -1
            return
        end

        buf = Array{Int32}(undef, 5) # [block_size, refID, pos, bin_mq_nl, flag_nc]
        gzread(f, pointer(buf), 20)
        block_size = buf[1]
        refID = buf[2] + 1 # the reference contig this read maps to

        # get the read position
        if refID != 0
            r.position = buf[3] + r.contigs.offsets[refID] + 1; # convert to 1 based indexing
        end

        forward = (buf[5] & 1048576) == 0 # see if we are reverse complemented
        skip(f, block_size-16) # skip the rest of the entry

        # break if we found a read in the right direction
        if refID != 0 && (r.readOrientation == :any || (forward && r.readOrientation == :forward) || (!forward && r.readOrientation == :reverse))
            return
        end
    end
end

# here we want to update the reader
#eachposition(r::BamReader) = BamReaderIterator(r)
#struct BamReaderIterator
#	reader::BamReader
#end
#Base.start(it::BamReaderIterator) = it.reader.position
#Base.done(it::BamReaderIterator, position) = position == -1
#function Base.next(it::BamReaderIterator, position)
#	pos = it.reader.position
#	advance!(it.reader)
#	pos,it.reader.position
#end
