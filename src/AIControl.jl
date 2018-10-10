module AIControl
    using Distributions
    using JLD2
    using FileIO
    using CSV
    using DataFrames
    using GZip
    using Statistics
    using LinearAlgebra

    include("./ReferenceContigs.jl")
    include("./BamReader.jl")
    include("./BinningMap.jl")
    include("./BinnedReader.jl")
    include("./DenseBlockIterator.jl")
    include("./MatrixReader.jl")
    include("./Utils.jl")
    include("./PeakWriter.jl")
    include("./PeakCaller.jl")
    include("./EvaluationUtil.jl")
end
