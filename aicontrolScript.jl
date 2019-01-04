## help statement

function printUsage()
    println("==================USAGE===================")
    println("julia aicontrolScript.jl [bamfile] [option1] [option2] ...")
    println("\t\t --dup: using duplicate reads [default:false]")
    println("\t\t --reduced: using subsampled control datasets [default:false]")
    println("\t\t --fused: fusing consecutive peaks [default:false]")
    #println("\t\t --xtxfolder=[path]: path to a folder with xtx.jld2 [default:./data]")
    println("\t\t --ctrlfolder=[path]: path to a control folder [default:./data]")
    println("\t\t --name=[string]: prefix for output files [default:bamfile_prefix]")
    println("\t\t --p=[float]: pvalue threshold [default:0.15]")
    println("")
    println("Example: julia aicontrolScript.jl test.bam --ctrlfolder=/scratch --name=test")
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

    global name = split(split(bamfilepath, "/")[end], ".")[1]
    temp = filter(x->occursin("--name", x), ARGS)
    if length(temp)>0
        global name = split(temp[1], "=")[2]
    end
    
    global ctrlfolder = "./data"
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
println("isDup : ", isDup)
println("isFull: ", isFull)
println("isFused: ", isFused)
println("prefix: ", name)
println("p-value (-log10)    : ", mlog10p)
println("path to control data: ", ctrlfolder)
#println("path to other data  : ", xtxfolder)
println("=========================================")

#check for file existance

if !isfile("$(ctrlfolder)/forward.data100$(fullstring)$(dupstring)")
    println(stderr, "$(ctrlfolder)/forward.data100$(fullstring)$(dupstring) missing.")
    println(stderr, "Please specify its location by --ctrlfolder=[path to the folder]")
    println(stderr, "Please read the step4 at https://github.com/hiranumn/AIControl.jl")
    printUsage()
    exit()
end

if !isfile("$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring)")
    println(stderr, "$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring) missing.")
    println(stderr, "Please specify its location by --ctrlfolder=[path to the folder]")
    println(stderr, "Please read the step4 at https://github.com/hiranumn/AIControl.jl")
    printUsage()
    exit()
end 

using Distributed
using JLD2
using FileIO
addprocs(2)
@everywhere using AIControl

# Checking progress
progress = 0
if isfile("$(name).jld2")
    tempdata = load("$(name).jld2")
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

if !(isfile("$(name).fbin100") && isfile("$(name).fbin100"))
    # Binning code
    @everywhere function wrapper1(args)
        write_binned(args[1], args[2], 100, args[3])
    end
    println("Binning files ...")
    pmap(wrapper1, [[bamfilepath, "$(name).fbin100", :forward], [bamfilepath, "$(name).rbin100", :reverse]])
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

    outcome = pmap(wrapper2, [["$(ctrlfolder)/forward.data100$(fullstring)$(dupstring)","$(name).fbin100","f", "xtxs$(fullstring)$(dupstring).jld2"],["$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring)","$(name).rbin100","r", "xtxs$(fullstring)$(dupstring).jld2"]])

    tempdata = Dict()
    tempdata["w1-f"] = outcome[1][1]
    tempdata["w2-f"] = outcome[1][2]
    tempdata["w1-r"] = outcome[2][1]
    tempdata["w2-r"] = outcome[2][2] 
    save("$(name).jld2", tempdata)
end

if progress < 2
    # Computing fits
    @everywhere function wrapper3(args)
        verbosity = 320
        _mr = MatrixReader(args[1], 10000)
        f = computeFits(_mr, args[3], args[2], verbose=verbosity)
    end
    println("Computing fits ...")

    outcome = pmap(wrapper3, [["$(ctrlfolder)/forward.data100$(fullstring)$(dupstring)","f", "$(name).jld2"],["$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring)","r", "$(name).jld2"]])

    tempdata = load("$(name).jld2")
    tempdata["fit-f"] = outcome[1]
    tempdata["fit-r"] = outcome[2]
    save("$(name).jld2", tempdata)
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

    outcome = pmap(wrapper4, [["$(name).fbin100","f", "$(name).jld2"],["$(name).rbin100","r", "$(name).jld2"]])

    tempdata = load("$(name).jld2")
    tempdata["p-f"] = outcome[1][1]
    tempdata["fold-f"] = outcome[1][2]
    tempdata["p-r"] = outcome[2][1]
    tempdata["fold-r"] = outcome[2][2]
    save("$(name).jld2", tempdata)
end

if progress < 4
    # Learning offset
    println("Estimating peak distance ...")
    offset = estimateD("$(name).fbin100", "$(name).rbin100")
    tempdata = load("$(name).jld2")
    tempdata["offset"] = offset
    save("$(name).jld2", tempdata)
end

###############
# Write peaks #
###############
println("Writing peaks out ...")
if !isFused
    test = generateUnfusedPeakFile("$(name).jld2", String("$(name)"), th=mlog10p)
else
    test = generatePeakFile("$(name).jld2", String("$(name)"), th=mlog10p)
end

println("Done. Peaks written to $(name).narrowPeak")







