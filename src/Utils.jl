export covariance, cov2cor!, filter2d, smooth, addConstColumn, weights_masker

##################################################
# Calculates covariance matrix for binary matrix 
##################################################
function covariance(bigData::BitArray{2}; chunkSize=100000, quiet=false)
    P,N = size(bigData)
    XtX = zeros(Float64, P, P)
    varSums = zeros(Float64, P, 1)
 
    # force the chunk size to line up with 64 bit word boundaries,
    # this is most important for loading from a file, but we also use it here.
    # we try and keep the size close to what was requested
    chunkSize = max(round(Int64, chunkSize/64),1)*64
    # build XtX incrementally and also the totals of every variable.
    chunk = Array(Float32, P, chunkSize)
    numChunks = round(Int64, ceil(N/chunkSize))
    for i in 1:numChunks-1
        chunk[:,:] = bigData[:,(i-1)*chunkSize+1:i*chunkSize]
        XtX .+= A_mul_Bt(chunk,chunk) # using a float array is important to get LAPACK speed
        varSums .+= sum(chunk,2)
        #if !quiet println(STDERR, "processed $(i*chunkSize*1000) bp...") end
    end
 
    # get the last unevenly sized chunk
    chunk = Array(Float32, P, N - (numChunks-1)*chunkSize)
    chunk[:,:] = bigData[:,(numChunks-1)*chunkSize+1:end]
    XtX .+= A_mul_Bt(chunk,chunk)
    varSums .+= sum(chunk,2)

    # convert XtX to a covariance matrix
    XtX .-= varSums*varSums'/N
    XtX ./= (N-1)
end

#####################################
# Converts covariance to correlation 
#####################################
function cov2cor!(M)
    for i in 1:size(M)[1]
        val = sqrt(M[i,i])
        if val > 0.0
            M[i,:] ./= val
            M[:,i] ./= val
        end
    end
    M
end

######################################
# Takes in matrix and true/fase masks 
######################################
function filter2d(matrix, ymask, xmask)
    temp = filter(x->x[2], [(matrix[i,1:end], ymask[i]) for i in 1:size(matrix)[1]])
    temp = [i[1] for i in temp]
    temp = hcat(temp...)
    
    temp = filter(x->x[2], [(temp[i,1:end], xmask[i]) for i in 1:size(temp)[1]])
    temp = [i[1] for i in temp]
    temp = hcat(temp...)
end

######################################
# Smoothing function for peak calling 
######################################
function smooth(a, width)
    if width > length(a)
        width = length(a)
    end
    
    ret = copy(a)
    counts = ones(length(a))
    
    #Aggregating the values
    for i in 1:width
        for j in 1:length(a)-i
            ret[j] += a[j+i]
            ret[end-j+1] += a[end-j-i+1]
            counts[j] += 1
            counts[end-j+1] += 1
        end 
    end
    ret./counts
end

######################################
# Adds constant column to data matrix
######################################
function addConstColumn(M::Array{Int64, 2})
    _const = ones(Int64, (size(M)[1], 1))
    hcat(M, _const)
end

##########################################################################
# masks weight based on some conditions 
# mask = weights_masker("ENCFF000YRS", 0, JLD.load(ctrllistdata)["ctrls"])
##########################################################################
function weights_masker(expID::String, mode::Int, clist; metadata::String="../data/metadata.csv", constant::Bool=true)
    m = CSV.read(metadata)
    
    # mode 0 removes matched control
    if mode == 0
        cexp = convert(String, m[m[:ID].==expID, :CTRL1][1])
        ignorelist = m[m[:EXP].==cexp, :ID]
        
    # mode 1 removes conrol form the same cellline
    elseif mode == 1
        ct = convert(String, m[m[:ID].==expID, :CELLTYPE][1])
        ignorelist = m[(m[:IFCTRL].==true).&(m[:CELLTYPE].==ct), :ID]
        
    # mode 2 removes conrol form the same lab
    elseif mode == 2
        lab = convert(String, m[m[:ID].==expID, :LAB][1])
        ignorelist = m[(m[:IFCTRL].==true).&(m[:LAB].==lab), :ID]
    
    # combination of mode 1 and 2
    elseif mode == 3
        ct = convert(String, m[m[:ID].==expID, :CELLTYPE][1])
        ignorelist1 = m[(m[:IFCTRL].==true).&(m[:CELLTYPE].==ct), :ID]
        lab = convert(String, m[m[:ID].==expID, :LAB][1])
        ignorelist2 = m[(m[:IFCTRL].==true).&(m[:LAB].==lab), :ID]
        ignorelist = append!(ignorelist1,ignorelist2)
    end
    
    ignorelist = collect(Set(ignorelist))

    mask = []
    for c in clist
        if c in ignorelist
            push!(mask, false)
        else
            push!(mask, true)
        end
    end

    # some assertion
    assert(length(ignorelist) >= length(mask)-sum(mask))

    if constant
        push!(mask, true)
    end

    convert(BitArray{1}, mask)
end
