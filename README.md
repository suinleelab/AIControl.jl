# AIControl

AIControl removes need of performing Input/control exeperiments for ChIP-seq assays by imputing them from public data.

![alt text](images/concept.png)

*Figure 1: (a) Comparison of AIControl to other peak calling algorithms. (left) AIControl
learns appropriate combinations of publicly available control ChIP-seq datasets to impute background
noise distributions at a fine scale. (right) Other peak calling algorithms use only one
control dataset, so they must use a broader region (typically within 5,000-10,000 bps) to estimate
background distributions. (bottom) The learned fine scale Poisson (background) distributions are
then used to identify binding activities across the genome. (b) An overview of the AIControl
approach. A single control dataset may not capture all sources of background noise. AIControl
more rigorously removes background ChIP-seq noise by using a large number of publicly available
control ChIP-seq datasets*

## Required libraries for AIControl

AIControl module is coded in **Julia 0.6.4**.

Before you start, make sure your have the following required modules.
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

## How to use

**1. Map your FASTQ file from ChIP-seq to the `hg38` assembly from UCSC.**  
   We have tested our pipeline with `bowtie2`. You can download the genome assembly data from [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)  
   *Example command:*\n
   `bowtie2 -x hg38 -q -p 10 -U ERR231591.fastq -S ERR231591.sam`  
   
**2. Convert the resulting sam file into a bam format.**  
*Example command:*  
`samtools view -Sb ERR231591.sam > ERR231591.bam`  
   
**3. Sort your bam file in lexicographical order.**  
   If you went through step 1 with the UCSC hg38 assembly, sorting with `samtools sort` will do its job.  
   *Example command:*  
   `samtools sort -o ERR231591.bam.sorted ERR231591.bam`  
   
   Sometimes your bam file is mapped to hg38, but to a slightly differet version or different ordering of chromosomes (a.k.a. non-lexicographic). For example, if you download a bam file directly from ENCODE portal, this is unfortunately the case. A recommended way of avoiding this problem is to extract a fastq file from your bam file, go back to step 1, and remap it with bowtie2 using the UCSC hg38 assembly. `bedtools` provide a way to generate a fastq file from your bam file.
   *Example command:*  
   `bedtools bamtofastq [OPTIONS] -i <BAM> -fq <FASTQ>`  
   
**4. Download data files and locate them in the right places.**  

