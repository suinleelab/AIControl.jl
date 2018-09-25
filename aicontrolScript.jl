## help statement

function printUsage()
    println("==================USAGE===================")
    println("julia aicontrolScript.jl [bamfile] [option1] [option2] ...")
    println("\t\t --dup: using duplicate reads [default:false]")
    println("\t\t --reduced: using subsampled control datasets [default:false]")
    println("\t\t --xtxfolder=[path]: path to a folder with xtx.jld [default:./data]")
    println("\t\t --ctrlfolder=[path]: path to a control folder [default:./data]")
    println("\t\t --name=[string]: prefix for output files [default:bamfile_prefix]")
    println("\t\t --p=[float]: pvalue threshold [default:0.15]")
    println("")
    println("Example: julia aicontrolScript.jl test.bam --dup --reduced --ctrlfolder=/scratch --name=test")
end

if "--help" in ARGS || "--h" in ARGS || length(ARGS)==0
    printUsage()
    quit()
end

## check for file existance
bamfilepath = ARGS[1]
if !isfile(bamfilepath)
    println("Input bam file does not exist.")
    println()
    printUsage()
    quit()
end

isDup = false
dupstring = ".nodup"

isFull = true
fullstring = ""

name = ""
xtxfolder = ""
ctrlfolder = ""
contigpath = ""

mlog10p = 1.5

try
    ## parsing arguments
    if "--dup" in ARGS
        isDup = true
        dupstring = ".dup"
    end

    if "--reduced" in ARGS
        isFull = false
        fullstring = ".reduced"
    end

    name = split(split(bamfilepath, "/")[end], ".")[1]
    temp = filter(x->contains(x, "--name"), ARGS)
    if length(temp)>0
        name = split(temp[1], "=")[2]
    end
    
    ctrlfolder = "./data"
    temp = filter(x->contains(x, "--ctrlfolder"), ARGS)
    if length(temp)>0
        ctrlfolder = split(temp[1], "=")[2]
    end
    
    xtxfolder = "./data"
    temp = filter(x->contains(x, "--xtxfolder"), ARGS)
    if length(temp)>0
        xtxfolder = split(temp[1], "=")[2]
    end
    
    # This option is currently not necessary...
    contigpath = "./aic.contig"
    temp = filter(x->contains(x, "--contig"), ARGS)
    if length(temp)>0
        contigpath = split(temp[1], "=")[2]
    end
    
    try
        temp = filter(x->contains(x, "--p"), ARGS)
        if length(temp)>0
            pthreshold = float(split(temp[1], "=")[2])
            mlog10p = -1*log10(pthreshold)
        end
    end
    
catch
    printUsage()
    quit()
end

println("============PARAMETERS====================")
println("isDup : ", isDup)
println("isFull: ", isFull)
println("prefix: ", name)
println("p-value (-log10)    : ", mlog10p)
println("path to control data: ", ctrlfolder)
println("path to other data  : ", xtxfolder)
println("=========================================")

#check for file existance
if !isfile("$(xtxfolder)/xtxs$(fullstring)$(dupstring).jld")
    println("xtx.jld file missing.")
    quit()
end

if !isfile("$(ctrlfolder)/forward.data100$(fullstring)$(dupstring)")
    println("forward.data100$(fullstring)$(dupstring) missing.")
    quit()
end

if !isfile("$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring)")
    println("forward.data100$(fullstring)$(dupstring) missing.")
    quit()
end 

addprocs(2)
@everywhere using AIControl

# Checking progress
progress = 0
if isfile("$(name).jld")
    tempdata = JLD.load("$(name).jld")
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
        verbosity = 100
        _mr = MatrixReader(args[1], 10000)
        _br = BinnedReader(args[2])
        w = computeBeta(_mr, _br, args[3], verbose=verbosity, xtxfile=args[4])
    end
    println("Computing weights ...")
    
    println("$(xtxfolder)/xtxs$(fullstring)$(dupstring).jld")

    outcome = pmap(wrapper2, [["$(ctrlfolder)/forward.data100$(fullstring)$(dupstring)","$(name).fbin100","f", "$(xtxfolder)/xtxs$(fullstring)$(dupstring).jld"],["$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring)","$(name).rbin100","r", "$(xtxfolder)/xtxs$(fullstring)$(dupstring).jld"]])

    JLD.save("$(name).jld", "w1-f", outcome[1][1], "w2-f",  outcome[1][2], "w1-r",  outcome[2][1], "w2-r",  outcome[2][2])
end

if progress < 2
    # Computing fits
    @everywhere function wrapper3(args)
        verbosity = 100
        _mr = MatrixReader(args[1], 10000)
        f = computeFits(_mr, args[3], args[2], verbose=verbosity)
    end
    println("Computing fits ...")

    outcome = pmap(wrapper3, [["$(ctrlfolder)/forward.data100$(fullstring)$(dupstring)","f", "$(name).jld"],["$(ctrlfolder)/reverse.data100$(fullstring)$(dupstring)","r", "$(name).jld"]])

    tempdata = JLD.load("$(name).jld")
    tempdata["fit-f"] = outcome[1]
    tempdata["fit-r"] = outcome[2]
    JLD.save("$(name).jld", tempdata)
end

if progress < 3
    # Calling peaks
    @everywhere function wrapper4(args)
        _br = BinnedReader(args[1])
        p, fold, t, l = callPeaks(_br, args[3], args[2], verbose=0)
        p, fold
    end
    println("Calling peaks ...")

    outcome = pmap(wrapper4, [["$(name).fbin100","f", "$(name).jld"],["$(name).rbin100","r", "$(name).jld"]])

    tempdata = JLD.load("$(name).jld")
    tempdata["p-f"] = outcome[1][1]
    tempdata["fold-f"] = outcome[1][2]
    tempdata["p-r"] = outcome[2][1]
    tempdata["fold-r"] = outcome[2][2]
    JLD.save("$(name).jld", tempdata)
end

if progress < 4
    # Learning offset
    println("Estimating peak distance ...")
    offset = estimateD("$(name).fbin100", "$(name).rbin100")
    tempdata = JLD.load("$(name).jld")
    tempdata["offset"] = offset
    JLD.save("$(name).jld", tempdata)
end

###############
# Write peaks #
###############
println("Writing peaks out ...")
test = generateUnfusedPeakFile("$(name).jld", String("$(name)"), th=mlog10p)

println("Done. Peaks written to $(name).peaks")







