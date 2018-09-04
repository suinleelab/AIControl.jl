# AIControl

AIControl makes ChIP-seq assays **easier**, **cheaper**, and **more accurate** by imputing background data from mass control data available in public.

Here is an overview of AIControl framework form our paper. 
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

## Installing utility softwares
These external softwares are required to run AIControl pipeline from a `.fastq` file to a `.narrowpeak` file.
- bowtie2
- samtools
- bedtools

The recommended way of installing these softwares is to use `conda`.  
Download annaconda Python distribution from [here](https://anaconda.org/anaconda/python).  
Install `annaconda` and run the following commands.
- ` conda install -c bioconda bowtie2 `
- ` conda install -c bioconda samtools `
- ` conda install -c bioconda bedtools `

## Julia modules required for AIControl

AIControl module is coded in **Julia 1.0 and 0.7**.  
You can download Julia from [here](https://julialang.org/).  

Before you start, make sure your have the following required modules.
Use `Pkg.add()` to install libraries.
- `Pkg.add("DataFrames")`
- `Pkg.add("JLD")`
- `Pkg.add("Distributions")`
- `Pkg.add("CSV")`

Then, you can simplely install AIControl module by typing 
- `Pkg.clone("https://github.com/hiranumn/AIControl.git")`

## Data files required for AIControl
AIControl uses a mass amount of public control data for ChIP-seq (roughly 450 chip-seq runs). We have done our best to compress them so that you only need to download about **4.6GB** (can be smaller with the `--reduced` option). These files require approximately **13GB** of free disk space to unfold. You can unfold them to anywhere you want as long as you specify the location with the `--ctrlfolder` option. **The default location is `./data`**

**[Here](https://drive.google.com/open?id=1Xh6Fjah1LoRMmbaJA7_FzxYcbqmpNUPZ) is a link to a Google Drive folder that contains all compressed data.** The link contains 8 files and you need to download different files depending on which option you are planning to run with.

### Option 1: with 440 ENCODE controls with non-duplicate reads
**Using these files will recreate the main results of the paper.**  
- `forward.data100.nodup.tar.bz2` (2.3GB):   
- `reverse.data100.nodup.tar.bz2` (2.3GB):  
When extracted, both files together occupy approximately 13GB of space. 

### Option 2: with 273 ENCODE controls with non-duplicate reads
We also have a subsampled version, which is validated to have comperable performance. 
- `forward.data100.reduced.nodup.tar.bz2` (1.5GB):   
- `reverse.data100.reduced.nodup.tar.bz2` (1.5GB):  
When extracted, both files together occupy approximately 8GB of space. 

We also have versions where duplicates are not removed (indicated with `.dup`, and used with the `--dup` option).

## Paper
We have an accompanying paper in BioRxiv evaluating and comparing the performance of AIControl to other peak callers in various metrics and settings. **AIControl: Replacing matched control experiments with machine learning improves ChIP-seq peak identification** ([BioRxiv](https://www.biorxiv.org/content/early/2018/03/08/278762?rss=1))

## How to use AIControl (step by step)

**Step 1: Map your FASTQ file from ChIP-seq to the `hg38` assembly from the UCSC database.**  
   We have validated our pipeline with `bowtie2`. You can download the genome assembly data from [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz). In case you need the exact reference database that we used for bowtie2, they are available through our [Google Drive](https://drive.google.com/open?id=1Xh6Fjah1LoRMmbaJA7_FzxYcbqmpNUPZ) as a zip file named `bowtie2ref.zip`.  
*Example command:*  
`bowtie2 -x hg38 -q -p 10 -U example.fastq -S example.sam`  
   
**Step 2: Convert the resulting sam file into a bam format.**  
*Example command:*  
`samtools view -Sb example.sam > example.bam`  
   
**Step 3: Sort the bam file in lexicographical order.**  
If you go through step 1 with the UCSC hg38 assembly, sorting with `samtools sort` will do its job.  
*Example command:*  
`samtools sort -o example.bam.sorted example.bam`  
   
Sometimes your bam file is mapped to hg38, but to a slightly differet version or different ordering of chromosomes (a.k.a. non-lexicographic). For instance, if you download a bam file directly from the ENCODE website, it is mapped to a slightly different version of hg38. A recommended way of resolving this issue is to extract a fastq file from your bam file, go back to step 1, and remap it with bowtie2 using the UCSC hg38 assembly. `bedtools` provides a way to generate a fastq file from your bam file. `samtools sort` does not quite do its job here sadly.
*Example command:*  
`bedtools bamtofastq  -i example.bam -fq example.fastq`  
   
**4. Download data files and locate them in the right places.**  
As stated, AIControl requires you to download precomputed data files. Please download and extract them to the `./data` folder, or otherwise specify the location with `--ctrlfolder` option.   

**5. Run AIControl as julia script.**  
Here is a sample command  
`julia aicontrolScript.jl test.bam --ctrlfolder=/scratch/hiranumn/data --xtxfolder=./data --name=test`

Do `julia aicontrolScript.jl --help` or `-h` for help.

We currently support the following flags. 

- `--dup`: using duplicate reads \[default:false\]
- `--reduced`: using subsampled control datasets \[default:false\]
- `--xtxfolder=[path]`: path to a folder with xtx.jld \[default:./data\]
- `--ctrlfolder=[path]`: path to a control folder \[default:./data\]
- `--name=[string]`: prefix for output files \[default:bamfile_prefix\]
- `--p=[float]`: pvalue threshold \[default:0.15\]

## Simple trouble shooting
- You are using Julia 0.6 and above.
- You downloaded necessary files for `--reduced` or `--dup` if you are running with those flags.
- You sorted the input bam files according to the UCSC hg38 assembly as specified in Step 3.  

If you have questions, please e-mail to hiranumn at cs dot washington do edu.

