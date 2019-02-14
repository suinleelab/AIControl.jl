#!/usr/bin/env julia

## help statement

function printUsage()
    println("==================USAGE===================")
    println("julia aicontrolScript.jl [bamfile] [option1] [option2] ...")
    println("\t\t --ctrlfolder=[path]: path to a control folder [default:.]")
    println("\t\t --name=[string]: prefix for output files [default:bamfile_prefix]")
    println("\t\t --p=[float]: -log10 pvalue threshold [default: 0.03 (or -log10(0.03)=1.5)]")
    println("\t\t --disableParallel: a flag to disable parallel processing [default:false]")
    println("\t\t --dup: a flag to use duplicate reads [default:false]")
    println("\t\t --reduced: a flag to use subsampled control datasets [default:false]")
    println("\t\t --fused: a flag to fuse consecutive peaks [default:false]")
    println("")
    println("julia aicontrolScript.jl example.sorted.bam --ctrlfolder=. --name=test")
end

if "--help" in ARGS || "--h" in ARGS || length(ARGS)==0
    printUsage()
    exit()
end

## check for file existance
bamfilepath = ARGS[1]
if !isfile(bamfilepath)
    println(stderr, "Input bam file does not exist.")
    printUsage()
    exit()
end

isDup = false
dupstring = ".nodup"

isFull = true
fullstring = ""

isFused = false

isParallelDisabled = false 

name = ""
#xtxfolder = ""
ctrlfolder = ""
contigpath = ""

mlog10p = 1.5

try
    ## parsing arguments
    if "--dup" in ARGS
        global isDup = true
        global dupstring = ".dup"
    end
    
    if "--fused" in ARGS
        global isFused = true
    end
    
    if "--reduced" in ARGS
        global isFull = false
        global fullstring = ".reduced"
    end
    
    if "--disableParallel" in ARGS
        global isParallelDisabled = true
    end

    global name = split(split(bamfilepath, "/")[end], ".")[1]
    global fileprefix = split(split(bamfilepath, "/")[end], ".")[1]
    temp = filter(x->occursin("--name", x), ARGS)
    if length(temp)>0
        global name = split(temp[1], "=")[2]
        global fileprefix = split(split(temp[1], "=")[2], "/")[end]
    end
    
    global ctrlfolder = "."
    temp = filter(x->occursin("--ctrlfolder", x), ARGS)
    if length(temp)>0
        global ctrlfolder = split(temp[1], "=")[2]
    end
    
    temp = filter(x->occursin("--p", x), ARGS)
    if length(temp)>0
        pthreshold = float(split(temp[1], "=")[2])
        global mlog10p = -1*log10(pthreshold)
    end
    
catch
    printUsage()
    exit()
end

println("============PARAMETERS====================")
println("path to control data: ", ctrlfolder)
println("prefix: ", name)
println("p-value (-log10)    : ", mlog10p)
println("Parallel: ", !isParallelDisabled)
println("isDup : ", isDup)
println("isFull: ", isFull)
println("isFused: ", isFused)
println("=========================================")

#check for file existance

if !isfile("$(ctrlfolder)/forward.data100$(fullstring)$(dupstring)")
    println(stderr, "$(ctrlfolder)/forward.data100$(fullstring)$(dupstring) missing.")
    println(stderr, "Please specify its location by --ctrlfolder=[path to the folder]")
    println(stderr, "Please read the step3 at https://github.com/hiranumn/AIControl.jl")
    printUsage()
    exit()
end

if !isfile("$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring)")
    println(stderr, "$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring) missing.")
    println(stderr, "Please specify its location by --ctrlfolder=[path to the folder]")
    println(stderr, "Please read the step3 at https://github.com/hiranumn/AIControl.jl")
    printUsage()
    exit()
end 

using Distributed
using JLD2
using FileIO

# Adding processes unless specifically told to disable parallel computing
if !isParallelDisabled 
    addprocs(2)
end
@everywhere using AIControl

# Making a directory
# Check if there is a directory already
if isdir("$(name)")
    println(stderr, "A directory with name \"$(name)\" already exist.")
    println(stderr, "AIControl will try to use intermediate files from the last run with the same name.")
    println(stderr, "If this is not a desired behavior, please give a different name with the --name option.")
else
    mkdir("$(name)")
end

# Checking progress
progress = 0
if isfile("$(name)/$(fileprefix).jld2")
    tempdata = load("$(name)/$(fileprefix).jld2")
    if "offset" in keys(tempdata)
        progress = 4
    elseif "fold-r" in keys(tempdata)
        progress = 3
    elseif "fit-r" in keys(tempdata)
        progress = 2
    elseif "w2-r" in keys(tempdata)
        progress = 1
    end
end

println("Progress: ", progress)

if !(isfile("$(name)/$(fileprefix).fbin100") && isfile("$(name)/$(fileprefix).fbin100"))
    println("Binning files ...")
    write_binned(bamfilepath, "$(name)/$(fileprefix).fbin100", 100, :forward)
    write_binned(bamfilepath, "$(name)/$(fileprefix).rbin100", 100, :reverse)
end

if progress < 1
    # Computing weights
    @everywhere function wrapper2(args)
        verbosity = 320
        _mr = MatrixReader(args[1], 10000)
        _br = BinnedReader(args[2])
        w = computeBeta(_mr, _br, args[3], verbose=verbosity, xtxfile=args[4])
    end
    println("Computing weights ...")

    outcome = pmap(wrapper2, [["$(ctrlfolder)/forward.data100$(fullstring)$(dupstring)","$(name)/$(fileprefix).fbin100","f", "xtxs$(fullstring)$(dupstring).jld2"],["$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring)","$(name)/$(fileprefix).rbin100","r", "xtxs$(fullstring)$(dupstring).jld2"]])

    tempdata = Dict()
    tempdata["w1-f"] = outcome[1][1]
    tempdata["w2-f"] = outcome[1][2]
    tempdata["w1-r"] = outcome[2][1]
    tempdata["w2-r"] = outcome[2][2] 
    save("$(name)/$(fileprefix).jld2", tempdata)
end

if progress < 2
    # Computing fits
    @everywhere function wrapper3(args)
        verbosity = 320
        _mr = MatrixReader(args[1], 10000)
        f = computeFits(_mr, args[3], args[2], verbose=verbosity)
    end
    println("Computing fits ...")

    outcome = pmap(wrapper3, [["$(ctrlfolder)/forward.data100$(fullstring)$(dupstring)","f", "$(name)/$(fileprefix).jld2"],["$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring)","r", "$(name)/$(fileprefix).jld2"]])

    tempdata = load("$(name)/$(fileprefix).jld2")
    tempdata["fit-f"] = outcome[1]
    tempdata["fit-r"] = outcome[2]
    save("$(name)/$(fileprefix).jld2", tempdata)
end

if progress < 3
    # Calling peaks
    @everywhere function wrapper4(args)
        verbosity = 320
        _br = BinnedReader(args[1])
        p, fold, t, l = callPeaks(_br, args[3], args[2], verbose=verbosity)
        p, fold
    end
    println("Calling peaks ...")

    outcome = pmap(wrapper4, [["$(name)/$(fileprefix).fbin100","f", "$(name)/$(fileprefix).jld2"],["$(name)/$(fileprefix).rbin100","r", "$(name)/$(fileprefix).jld2"]])

    tempdata = load("$(name)/$(fileprefix).jld2")
    tempdata["p-f"] = outcome[1][1]
    tempdata["fold-f"] = outcome[1][2]
    tempdata["p-r"] = outcome[2][1]
    tempdata["fold-r"] = outcome[2][2]
    save("$(name)/$(fileprefix).jld2", tempdata)
end

if progress < 4
    # Learning offset
    println("Estimating peak distance ...")
    offset = estimateD("$(name)/$(fileprefix).fbin100", "$(name)/$(fileprefix).rbin100")
    tempdata = load("$(name)/$(fileprefix).jld2")
    tempdata["offset"] = offset
    save("$(name)/$(fileprefix).jld2", tempdata)
end

###############
# Write peaks #
###############
println("Writing peaks out ...")
if !isFused
    test = generateUnfusedPeakFile("$(name)/$(fileprefix).jld2", String("$(name)/$(fileprefix)"), th=mlog10p)
else
    test = generatePeakFile("$(name)/$(fileprefix).jld2", String("$(name)/$(fileprefix)"), th=mlog10p)
end

println("Done. Peaks written to $(name)/$(fileprefix).narrowPeak")







