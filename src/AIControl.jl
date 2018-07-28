module AIControl
    using Distributions
    using JLD
    using CSV
    using DataFrames
    using GZip

    include("./ReferenceContigs.jl")
    include("./BamReader.jl")
    include("./BinningMap.jl")
    include("./BinnedReader.jl")
    include("./DenseBlockIterator.jl")
    include("./MatrixReader.jl")
    include("./Utils.jl")
    include("./PeakWriter.jl")
    include("./PeakCaller.jl")
end