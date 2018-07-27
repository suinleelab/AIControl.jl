function checkControlAvailability(basedir, binsize=100, excludes=[])
    controlFiles = collect(Set([i[1:end-8] for i in readdir(basedir)]))
    availableControls = String[]
    mask = Bool[]
    for c in controlFiles
        if isfile("$(basedir)$(c).fbin$(binsize)") && isfile("$(basedir)$(c).rbin$(binsize)")
            if !(c[1:end-4] in excludes)
                push!(availableControls, c)
                push!(mask, true)
            else
                push!(mask, false)
            end
        end
    end
    availableControls, mask
end

function loadControls(controlFiles, basedir, direction, suffix)
    readers = BinnedReader[]
    control_names = Any[]
    for i in 1:length(controlFiles)
        filename = "$(basedir)$(controlFiles[i])$(direction)$(suffix)"
        push!(readers, BinnedReader(filename))
        push!(control_names, filename)
    end
    readers, control_names
end