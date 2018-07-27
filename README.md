# AIControl

## Required libraries for AIControl
AIControl module is coded in **Julia 0.6**.

Following modules are required.
Use `Pkg.add()` to install libraries.
- Pkg.add("DataFrames")
- Pkg.add("JLD")
- Pkg.add("Distributions")
- Pkg.add("CSV")

## Required data files
AIControl requires some precomputed data files.
- `forward.data100.tar.bz2` (2.3GB): compressed binned (100bps) signals of forward reads for all 440 control experiments.
- `reverse.data100.tar.bz2` (2.3GB): compressed binned (100bps) signals of reverse reads for all 440 control experiments.
- `xtxs.jld` (6.0MB): A pre-computed transpose(X)\*X matrix of size 441 by 441. 

When you uncompress the tar.bz2 files, each requires around 6.5GB of your HDD.  
We have a subsampled version (1.4GB compressed and 4.6GB uncompressed) for those who need smaller data.

## Paper
We have a paper in BioRxiv evaluating and comparing the performance of AIControl to other peak callers.  
**AIControl:  Replacing matched control experiments with machine learning improves ChIP-seq peak identification** ([BioRxiv](https://www.biorxiv.org/content/early/2018/03/08/278762?rss=1))


