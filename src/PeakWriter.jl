import Base: close
export PeakWriter, close, writePeak, sortPeaks, PeakWriter_unfused, WritePeak_unfused

mutable struct PeakWriter
    #a output stream object to write to
    outstream
    contigs::ReferenceContigs
    cur_ref::Int64
    peakid::Int64
    lastpos::Int64
end

function PeakWriter(output_stream, contigs)
    pw = PeakWriter(output_stream, contigs, 1, 1, 0)
    pw
end

close(pw::PeakWriter) = close(pw.outstream)

function writePeak(pw::PeakWriter, binSize::Int64, binPosStart::Int64, binPosEnd::Int64, pval, fold; prescision=2)
    
    # Some assertion to prevent misuse.
    @assert pw.lastpos<binPosStart
    
    # Get the starting position and ending position
    startPos = binSize*(binPosStart-1)+1
    endPos = binSize*binPosEnd
    
    # If your starting position goes over the size of current chromosome.
    while startPos > pw.contigs.offsets[pw.cur_ref]+pw.contigs.sizes[pw.cur_ref]
        pw.cur_ref += 1
    end
    
    # Current chrom name
    chr = pw.contigs.names[pw.cur_ref]
    # Get the current position
    startPos = startPos-pw.contigs.offsets[pw.cur_ref]
    endPos = endPos-pw.contigs.offsets[pw.cur_ref]
    # Name peak
    peakname = "peak_$(pw.peakid)"
    # Calculate score (cutoff fold at 1000)
    score = minimum([1000, fold])
    
    # Write it to file
    output = "$(chr)\t$(startPos)\t$(endPos)\t$(peakname)\t$(round(score,digits=prescision))\t.\t$(round(fold,digits=prescision))\t$(round(pval,digits=prescision))\t-1\t-1"
    println(pw.outstream, output)
    
    # Update some data
    pw.peakid += 1
    pw.lastpos = binPosEnd
    
    output
end


#########################################################################
# Preprocesses peak array into list of peaks for peak writer to take in #
#########################################################################
function sortPeaks(pvals, folds, th::Float64)
    flag = false
    pval_list = []
    fold_list = []
    startPos = 0
    peaks = []
    
    # parse through pvals
    for i in 1:length(pvals)
        if pvals[i] > th
            
            if flag == false
                # record starting position
                startPos = i
            end
            
            # set flag
            flag = true
            # record pvals
            push!(pval_list, pvals[i])
            push!(fold_list, folds[i])
        elseif pvals[i] < th && flag
            # set off flag
            flag = false
            # determine max pvals
            maxpval = maximum(pval_list)
            maxfold = maximum(fold_list)
            # reset lists
            pval_list = []
            fold_list = []
            # get ending pos
            endPos = i-1
            # record peaks
            # println([startPos, endPos, maxpval, maxfold])
            push!(peaks, (Int(startPos), Int(endPos), maxpval, maxfold))
        end
    end
    
    # finish off a peak
    if flag
        maxpval = maximum(pval_list)
        maxfold = maximum(fold_list)
        push!(peaks, (Int(startPos), length(pvals), maxpval, maxfold))
    end
        
    peaks
end

mutable struct PeakWriter_unfused
    #a output stream object to write to
    Outstream
    contigs::ReferenceContigs
    cur_ref::Int64
    id::Int64
end

function PeakWriter_unfused(output_stream, contigs)
    sw = PeakWriter_unfused(output_stream, contigs, 1, 1)
    sw
end

function WritePeak_unfused(sw::PeakWriter_unfused, binSize::Int64, binPosStart::Int64, binPosEnd::Int64, pval, fold; prescision=3)
    startPos = binSize*(binPosStart-1)+1
    endPos = binSize*binPosEnd
    
    while startPos > sw.contigs.offsets[sw.cur_ref]+sw.contigs.sizes[sw.cur_ref]
        sw.cur_ref += 1
    end
    
    chr = sw.contigs.names[sw.cur_ref]
    startPos = startPos-sw.contigs.offsets[sw.cur_ref]
    endPos = endPos-sw.contigs.offsets[sw.cur_ref]
    peakname = "peak_$(sw.id)"
    sw.id += 1
    score = fold
    
    output = "$(chr)\t$(startPos)\t$(endPos)\t$(peakname)\t$(round(score,digits=prescision))\t.\t$(round(fold,digits=prescision))\t$(round(pval,digits=prescision))\t-1\t-1"

    println(sw.Outstream, output) 
end
