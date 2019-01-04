export computeXtX, computeBeta, computeFits, estimateD, callPeaks, generatePeakFile, generateUnfusedPeakFile

#########################################################################################################
# Computes XtX for linear regression
#
# Example usage:
# out = computeXtX(MatrixReader("/scratch/hiranumn/forward.data100", 10000))
# out1 = computeXtX(MatrixReader("/scratch/hiranumn/reverse.data100", 10000))
# JLD.save("../data/xtxs.jld", "XtX1-f", out[1], "XtX2-f", out[2], "XtX1-r", out1[1], "XtX2-r", out1[2])
#########################################################################################################

function computeXtX(mr::MatrixReader; num_chroms=0, verbose=0, binsize=100)
    
    ##############################
    # Use all chroms for default #
    ##############################
    if num_chroms > length(ReferenceContigs_hg38.sizes) || num_chroms < 1
        num_chroms = length(ReferenceContigs_hg38.sizes)
    end
    training_limit = Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize))
    
    ################
    # Compute XtXs #
    ################
    XtX1 = zeros(Float64, mr.expsize+1, mr.expsize+1)
    XtX2 = zeros(Float64, mr.expsize+1, mr.expsize+1)
    
    advance!(mr)
    count = 0
    while !eof(mr) && count*mr.blocksize < training_limit
        
        # compute
        ctrl = convert(Array{Float64,2}, addConstColumn(mr.data)')
        if count%2==0
            BLAS.syrk!('U', 'N', 1.0, ctrl, 1.0, XtX1)
        else
            BLAS.syrk!('U', 'N', 1.0, ctrl, 1.0, XtX2)
        end
        
        # report progress
        if verbose>0 && count % (verbose) == 0
            println(count, ":", training_limit)
        end
        
        # update
        advance!(mr)
        count += 1
    end
    
    #############################
    # Converting to full matrix #
    #############################
    XtX1 = XtX1+XtX1'
    XtX2 = XtX2+XtX2'
    for i in 1:size(XtX1)[1]
       XtX1[i, i] ./= 2 
       XtX2[i, i] ./= 2 
    end
    
    XtX1, XtX2
end

#########################################################################################################
# Computes Beta for specific target
#
# Example usage:
# verbosity = 100
#
# _mr = MatrixReader("/scratch/hiranumn/forward.data100", 10000)
# _br = BinnedReader("/scratch/hiranumn/target_data_nodup/ENCFF000YRS.bam.fbin100")
#
# w_f = computeBeta(_mr, _br, "f", verbose=verbosity, xtxfile="../data/xtxs.jld")
#
# _mr = MatrixReader("/scratch/hiranumn/reverse.data100", 10000)
# _br = BinnedReader("/scratch/hiranumn/target_data_nodup/ENCFF000YRS.bam.rbin100")
#
# w_r = computeBeta(_mr, _br, "r", verbose=verbosity, xtxfile="../data/xtxs.jld")
#
# JLD.save("ENCFF000YRS.jld", "w1-f", w_f[1], "w2-f", w_f[2], "w1-r", w_r[1], "w2-r", w_r[2])
#########################################################################################################

function computeBeta(mr::MatrixReader, br::BinnedReader, direction::String; binsize=100, num_chroms=0, verbose=0, mask=[], xtxfile="../data/xtx.jld")

    ##############################
    # Use all chroms for default #
    ##############################
    if num_chroms > length(ReferenceContigs_hg38.sizes) || num_chroms < 1
        num_chroms = length(ReferenceContigs_hg38.sizes)
    end
    training_limit = Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize))
    
    ###########################################
    # Prepare matrices for weight calculation #
    ###########################################
    datapath = joinpath(@__DIR__, "..", "data")
    XtX1 = load(joinpath(datapath, xtxfile))["XtX1-$(direction)"]
    XtX2 = load(joinpath(datapath, xtxfile))["XtX2-$(direction)"]
    if length(mask)>0
        # mask input if necessary
        @assert mr.expsize+1==length(mask)
        XtX1 = filter2d(XtX1, mask, mask)
        XtX2 = filter2d(XtX2, mask, mask)
    end
    Xty1 = zeros((size(XtX1)[1],1))
    Xty2 = zeros((size(XtX2)[1],1))
    
    ###################
    # Compute Weights #
    ###################
        
    count = 0
    for target in denseblocks([br], mr.blocksize, constantColumn=false, loop=true)
        
        # update
        count += 1
        if count*mr.blocksize > training_limit break end
        
        # load control
        advance!(mr)
        ctrl = convert(Array{Float64,2}, addConstColumn(mr.data)')
        if length(mask)>0 ctrl = ctrl[mask, :] end
        
        # compute
        if count % 2 == 0
            Xty1 .+= ctrl*target
        else
            Xty2 .+= ctrl*target
        end
		
        # report progress
	if verbose>0 && count % (verbose) == 0
        	progress = Int(floor((count*mr.blocksize/training_limit)*1000))/10
        	printString = "$(progress)% completed ($(count*mr.blocksize)/$(training_limit))"
        	if direction=="f"
            		printString = printString*" on forward signals."
        	else
            		printString = printString*" on reverse signals."
        	end
        	println(printString)
    	end
    end
    
    m = size(XtX1)[1]
    beta1 = inv(XtX1 + 0.00001*Matrix(1.0I, m, m))*Xty1
    beta2 = inv(XtX2 + 0.00001*Matrix(1.0I, m, m))*Xty2
    
    beta1, beta2
end

###############################################################################
# Computes fits for specific target
#
# Example usage:
# verbosity = 100
#
# _mr = MatrixReader("/scratch/hiranumn/forward.data100", 10000)
# ff = computeFits(_mr, "ENCFF000YRS.jld", "f", verbose=verbosity)
# 
# _mr = MatrixReader("/scratch/hiranumn/reverse.data100", 10000)
# fr = computeFits(_mr, "ENCFF000YRS.jld", "r", verbose=verbosity)
# 
# JLD.save("ENCFF000YRS_fit.jld", "fit-f", ff, "fit-r", fr)
###############################################################################

function computeFits(mr::MatrixReader, weightfile::String, direction::String; binsize=100, num_chroms=0, verbose=0, mask=[])
    
    ##############################
    # Use all chroms for default #
    ##############################
    if num_chroms > length(ReferenceContigs_hg38.sizes) || num_chroms < 1
        num_chroms = length(ReferenceContigs_hg38.sizes)
    end
    training_limit = Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize))
    
    
    weight1 = load(weightfile)["w1-$(direction)"]
    weight2 = load(weightfile)["w2-$(direction)"]
    
    ##########################
    # Compute regression fit #
    ##########################
    regression_fit = zeros(Float16, training_limit)
     
    advance!(mr) 

    count = 0
    binpos = 0
    while !eof(mr) && count*mr.blocksize < training_limit
        
        # get data
        ctrl = convert(Array{Float64,2}, addConstColumn(mr.data)')
        if length(mask)>0 ctrl = ctrl[mask, :] end
        
        # compute
        if count % 2 == 0 
            pred = ctrl'*weight2
        else
            pred = ctrl'*weight1
        end
        
        # record
        for j in 1:length(pred)
            binpos +=1
            try 
                regression_fit[binpos] = pred[j]
            catch
            end
        end   
        
        # report progress
	if verbose>0 && (count+1) % (verbose) == 0
        	progress = Int(floor(((count+1)*mr.blocksize/training_limit)*1000))/10
        	printString = "$(progress)% completed ($((count+1)*mr.blocksize)/$(training_limit))"
        	if direction=="f"
            		printString = printString*" on forward signals."
        	else
            		printString = printString*" on reverse signals."
        	end
        	println(printString)
    	end
        
        advance!(mr)   
        # update
        count += 1
    end
    
    regression_fit
end

############################################################################
# Estimates distance between forward and reverse reads for specific target #
############################################################################

function estimateD(forwardtarget, reversetarget; binsize=100)
    
    #########################
    # Load in target values #
    #########################
    b = BinnedReader(forwardtarget)
    targetf = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes)/binsize)))
    while !eof(b) && position(b) < length(targetf)
        targetf[position(b)] = value(b)
        advance!(b)
    end
    close(b)
    b = BinnedReader(reversetarget)
    targetr = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes)/binsize)))
    while !eof(b) && position(b) < length(targetr)
        targetr[position(b)] = value(b)
        advance!(b)
    end
    close(b)

    ##################################################
    # Figure out distance that gives minimum overlap #
    ##################################################
    d = []
    for i in 0:4
        f = vcat([0 for j in  1:i], targetf[1:end-i])
        @assert f!=targetr
        @assert length(f)==length(targetr)
        push!(d, sum(abs.(f-targetr)))
    end
    argmin(d)-1
end

# Calls peak for both forward and reverse strands
function callPeaks(br::BinnedReader, fitfile::String, direction::String; num_chroms=0, verbose=0, binsize=100, base=1, smoothing=true)
    
    if num_chroms > length(ReferenceContigs_hg38.sizes) || num_chroms < 1
        num_chroms = length(ReferenceContigs_hg38.sizes)
    end
    training_limit = Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize))
    
    # fill in target vector
    target = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)))
    while !eof(br) && position(br) < length(target)
        target[position(br)] = value(br)
        advance!(br)
    end
    close(br)
    regfit = load(fitfile)["fit-$(direction)"]
    
    if verbose>0 println("Loaded peak signals.") end
    
    # Do smoothing if necessary
    if smoothing
        smooth1000 = smooth(regfit, 10)
        smooth5000 = smooth(regfit, 50)
        smooth10000 = smooth(regfit, 100)
    end
    m = mean(target)
    
    #if verbose>0 println("smoothed.") end
    
    # Recording vector
    pvals = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)))
    folds = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)))
    lambdas = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)))
    
    # Compute p-values and folds
    for i in 1:length(regfit)
        if smoothing
            lambda = maximum([regfit[i], smooth1000[i], smooth5000[i], smooth10000[i], base])
        else
            lambda = maximum([regfit[i], base])
        end
        pval = -1*log(10, 1-cdf(Poisson(lambda), target[i]))
        fold = target[i]/lambda
        if pval == Inf
            pval = Float32(typemax(Int64))
        end
        
        #This version of the code will assign pval to individual bins.
        pvals[i] = pval
        folds[i] = fold
        lambdas[i] = lambda
		
	# report progress
	if verbose>0 && (i/10000) % (verbose) == 0
        	progress = Int(floor((i/training_limit)*1000))/10
        	printString = "$(progress)% completed ($(i)/$(training_limit))"
        	if direction=="f"
            		printString = printString*" on forward signals."
        	else
            		printString = printString*" on reverse signals."
        	end
        	println(printString)
    	end
		
    end
    
    pvals, folds, target, lambdas
end

# combines reverse and forward peaks
function generatePeakFile(pfile::String, name::String; th=1.5, binsize=100)
    
    # create a vector final folds and pvals 
    pvals = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes)/binsize)))
    folds = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes)/binsize)))
    
    # get signals for each direction
    saved_data = load(pfile)
    forward_p = saved_data["p-f"]
    reverse_p = saved_data["p-r"]
    folds_f = saved_data["fold-f"]
    folds_r = saved_data["fold-r"]
    offset = saved_data["offset"]
    
    # Figuring out forward and reverse offset
    forward_offset = Int(ceil(offset/2))
    reverse_offset = Int(floor(offset/2))
    
    #Write into pvals and folds vectors from forward signals
    for i in 1:length(forward_p)
        try
            pvals[i+forward_offset] = forward_p[i]
            folds[i+forward_offset] = folds_f[i]
        catch
        end
    end
    
    #For reverse signals
    for i in 1:length(reverse_p)
        try
            if reverse_p[i] < pvals[i-reverse_offset] 
                pvals[i-reverse_offset] = reverse_p[i]
            end
            if folds_r[i] < folds[i-reverse_offset] 
                folds[i-reverse_offset] = folds_r[i]
            end
        catch
        end
    end

    
    #Get peaks to write on files
    processed_peaks = sortPeaks(pvals, folds, th)
    pw = PeakWriter(open("$(name).narrowPeak", "w"), ReferenceContigs_hg38)
    written = []
    for p in processed_peaks
        push!(written, writePeak(pw, 100, Int(p[1]), Int(p[2]), p[3], p[4]))
    end
    close(pw)
    written
end


function generateUnfusedPeakFile(pfile::String, name::String; th=1.5, binsize=100)
    
    # Data loading
    saved_data = load(pfile)
    forward_p = saved_data["p-f"]
    reverse_p = saved_data["p-r"]
    folds_f = saved_data["fold-f"]
    folds_r = saved_data["fold-r"]
    offset = saved_data["offset"]
        
    pvals = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes)/binsize)))
    folds = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes)/binsize)))
    
    forward_offset = Int(ceil(offset/2))
    reverse_offset = Int(floor(offset/2))
    
    for i in 1:length(forward_p)
        try
            pvals[i+forward_offset] = forward_p[i]
            folds[i+forward_offset] = folds_f[i]
        catch
        end
    end
    
    for i in 1:length(reverse_p)
        try
            if reverse_p[i] < pvals[i-reverse_offset] 
                pvals[i-reverse_offset] = reverse_p[i]
            end
            if folds_r[i] < folds[i-reverse_offset] 
                folds[i-reverse_offset] = folds_r[i]
            end
        catch
        end
    end
    
    fout = open("$(name).narrowPeak","w")
    pw = PeakWriter_unfused(fout, ReferenceContigs_hg38)
    for i in 1:length(pvals)
        #This version of the code will assign pval to individual bins.
        pval = pvals[i]
        if pval >= th
            WritePeak_unfused(pw, 100, i, i, pval, folds[i]) 
        end
    end
    close(fout)
    pvals
end
