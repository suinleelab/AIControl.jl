export load_narrowpeak, filterPreds, combinePeaks, binarizePeaks, window_bed_file

function load_narrowpeak(stream, contigs, index; binSize=1000, loadp=true, mlogt=false, verbose=0)
    # Get number of bins
    numBins = ceil(Int64, sum(contigs.sizes) / binSize)
    # Create offsets difctionary
    chrOffsets = Dict{String,Int64}()
    for i in 1:contigs.count
        chrOffsets[contigs.names[i]] = contigs.offsets[i]
    end

    # mark all bins that are touched with 1
    binValues = falses(numBins)
    # also record the highest p-value
    pValues = zeros(numBins)
    # also keep track how many peaks were evaluated
    count = 0 
    for line in eachline(stream)
        parts = split(rstrip(line), '\t')
        if haskey(chrOffsets, parts[1])
            count += 1
            if count < verbose
                println(parts)
            end
            startPos = ceil(Int64, (chrOffsets[parts[1]]+Int(parse(Float64, parts[2])))/binSize)
            endPos = ceil(Int64, (chrOffsets[parts[1]]+Int(parse(Float64, parts[3])))/binSize)
            for i in startPos:endPos
                assert(i!=0)
                # record bin that it was touched
                binValues[i] = true
                # loadp
                if loadp
                    # minus log 10 p-value transformation if necessary.
                    if mlogt
                        pval = -1*log(10, parse(Float64, parts[index]))
                    else
                        pval = parse(Float64, parts[index])
                    end
                    # record p-value if more significant.
                    if pValues[i] < pval
                        pValues[i] = pval
                    end
                end
            end
        end
    end
    close(stream)
    
    #assert(sum(binValues) == sum([i>0 for i in pValues]))
    binValues, pValues
end

function window_bed_file(stream, contigs; binSize=1000)
    numBins = ceil(Int64, sum(contigs.sizes) / binSize)
    chrOffsets = Dict{String,Int64}()
    for i in 1:contigs.count
        chrOffsets[contigs.names[i]] = contigs.offsets[i]
    end

    # mark all bins that are touched with 1
    binValues = falses(numBins)
    for line in eachline(stream)
        parts = split(line, '\t')
        if haskey(chrOffsets, parts[1])
            startPos = ceil(Int64, (chrOffsets[parts[1]]+parse(Int64, parts[2]))/binSize)
            endPos = ceil(Int64, (chrOffsets[parts[1]]+parse(Int64, parts[3]))/binSize)
            for i in startPos:endPos
                binValues[i] = true
            end
        end
    end

    binValues
end

# Get top x predictions with its truth values.
# Returns truth and preds. 
function filterPreds(top, truth, pred)
    temp = Any[]
    for i in 1:length(truth)
        push!(temp, (truth[i], pred[i]))
    end
    temp = sort(temp, by= x -> x[2], rev=true)[1:top]
    [i[1] for i in temp], [i[2] for i in temp]
end

# Combine signal vector a and b by taking minimum element wise.
# Returns a vector.
function combinePeaks(a, b)
    ret = zeros(length(a))
    for i in 1:length(a)
        ret[i] = minimum((a[i], b[i])) 
    end
    ret
end

# Extract top X peaks and output binarized vector of same size.
# Returns a vector.
function binarizePeaks(top, pred)
    temp = Any[]
    for i in 1:length(pred)
        push!(temp, (i, pred[i])) 
    end
    sort!(temp, by=x->x[2], rev=true)
    
    ret = falses(length(pred))
    for i in 1:top
        ret[temp[i][1]] = true
    end
    ret
end