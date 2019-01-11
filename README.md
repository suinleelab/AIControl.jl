# AIControl.jl

[![Build Status](https://travis-ci.org/hiranumn/AIControl.jl.svg?branch=master)](https://travis-ci.org/hiranumn/AIControl.jl)

AIControl makes ChIP-seq assays **easier**, **cheaper**, and **more accurate** by imputing background data from a massive amount of publicly available control data.

Here is an overview of the AIControl framework from our paper. 
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

## Major Updates
- (12/14/2018) Cleared all deprecations. AIControl now works with Julia 1.0. Please delete the precompiled cache from the previous versions of AIControl. You may do so by deleting the `.julia` folder. 
- (12/15/2018) Updated some error messages to better direct users (12/13/2018).
- (1/7/2019) Made AIControl Pkg3 compatible for Julia 1.0.3

## Installation 
AIControl can be used on any **Linux** or **macOS** machine. While we tested and validated that AIControl works on **Windows** machines, we believe that it is easier for you to set up the AIControl pipeline on the Unix based systems.

AIControl expects a sorted `.bam` file as an input and outputs a `.narrowpeak` file. Typically, for a brand new ChIP-seq experiment, you start with a `.fastq` file, and you will need some external softwares to convert the `.fastq` file to a sorted `.bam` file. Thus, the whole AIControl pipeline needs the following sets of programs and packages installed on your local machine. We will explain how to install them in sections below.
- `Julia (Julia 1.0 and above)`
- `bowtie2`: for aligning a `.fastq` file to the hg38 genome
- `samtools`: for sorting an alinged bam file
- `bedtools`: for converting a bam file back to a fastq file (OPTIONAL for Step 3.1)

### 1a. Installing Julia 1.0 for a Linux machine
The terminal commands below will install julia 1.0.3 on a linux machine. Please change the url accordingly. You can also download julia [here](https://julialang.org/downloads/). **[CAUTION:] We highly recommend avoiding the conda version of julia** as it currently known to have a problem locating libLLVM.so in many environments.
```
cd
wget https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.3-linux-x86_64.tar.gz
tar xvzf julia-1.0.3-linux-x86_64.tar.gz
echo  'export PATH=$PATH:~/julia-1.0.3/bin' >> ~/.bashrc
source ~/.bashrc
```

### 1b. Installing Julia 1.0 for a mac OS machine
Please first download the `.dmg` file for mac OS from the [julia website](https://julialang.org/downloads/), double-click to open it, and drag the icon to the Applications folder. Then, the following terminal command will put julia in your `PATH` and make it executable from command line. **[CAUTION:] We highly recommend avoiding the conda version of julia** as it currently known to have a problem locating libLLVM.so in many environments.

```
echo 'export PATH="/Applications/Julia-1.0.app/Contents/Resources/julia/bin/:${PATH}"' >> ~/.bash_profile
source ~/.bash_profile
```
See [this](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started#On_macOS_X) for more trouble shooting. 

### 2. Installing Julia Packages
The terminal command below will install required julia packages and AIControl.
```
julia -e 'using Pkg; Pkg.add(["FileIO", "JLD2"]); Pkg.add(PackageSpec(url = "https://github.com/hiranumn/AIControl.jl")); using AIControl'
```

### 3. Installing external softwares with miniconda
Please download and install miniconda from [here](https://conda.io/miniconda.html). The terminal command below will install required external softwares using conda package management system.
```
conda install -c bioconda bowtie2 samtools bedtools
```

## Control data files required for AIControl
AIControl uses a massive amount of public control data for ChIP-seq (roughly 450 chip-seq runs). We have done our best to compress them so that you only need to download about **4.6GB**. These files require approximately **13GB** of free disk space to unfold. The following terminal commands will download and decompress the compressed control data.
```
wget https://dada.cs.washington.edu/aicontrol/forward.data100.nodup.tar.bz2
tar xvjf forward.data100.nodup.tar.bz2
wget https://dada.cs.washington.edu/aicontrol/reverse.data100.nodup.tar.bz2
tar xvjf reverse.data100.nodup.tar.bz2
```
You can also obtain the control files from [our data repository](https://dada.cs.washington.edu/aicontrol/) or [Google Drive](https://drive.google.com/open?id=1Xh6Fjah1LoRMmbaJA7_FzxYcbqmpNUPZ).

## Paper
We have an accompanying paper in BioRxiv evaluating and comparing the performance of AIControl to other peak callers in various metrics and settings. **AIControl: Replacing matched control experiments with machine learning improves ChIP-seq peak identification** ([BioRxiv](https://www.biorxiv.org/content/early/2018/03/08/278762?rss=1)). You can find the supplementary data files and peaks files generated by the competing peak callers on [Google Drive](https://drive.google.com/open?id=1Xh6Fjah1LoRMmbaJA7_FzxYcbqmpNUPZ).

## Running AIControl (step by step)

### Step 0: Download a toy example.
The terminal command below will download a `.fastq` file that you may use as a toy example.  
They are also available at [our data repository](https://dada.cs.washington.edu/aicontrol/).
```
wget https://dada.cs.washington.edu/aicontrol/example.fastq
```

### Step 1: Map your FASTQ file from ChIP-seq to the `hg38` assembly from the UCSC database.
The following terminal commands will a) download and untar the reference database file for `bowtie2` and b) run `bowtie2` to map a `.fastq` file to the UCSC hg38 genome, which is available at [the UCSC repository](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz). 
```
wget https://dada.cs.washington.edu/aicontrol/bowtie2ref.tar.bz2
tar xvjf bowtie2ref.tar.bz2
bowtie2 -x bowtie2ref/hg38 -q -p 10 -U example.fastq -S example.sam
````  
Unlike other peak callers, the core idea of AIControl is to leverage all available control datasets. This requires all data (your target and public control datasets) to be mapped to the exact same reference genome. Our control datasets are currently mapped to the hg38 assembly from [the UCSC repository](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz). **So please make sure that your data is also mapped to the same assembly**. Otherwise, our pipeline will report an error.
   
### Step 2: Convert the resulting sam file into a bam format.  
```
samtools view -Sb example.sam > example.bam
```  
   
### Step 3: Sort the bam file in lexicographical order.
```
samtools sort -o example.bam.sorted example.bam
```  

### [Optional] Step 3.1: If AIControl reports an error for a mismatch of genome assembly.
You are likely here, because the AIControl script raised an error that look like the screenshot below. Otherwise, please move on to Step 4.

<img src="images/error3_1.png" alt="alt text" width="500"/>

The error is most likely caused by a mismatch of genome assembly that your dataset and control datasets are mapped to. Our control datasets are mapped to the hg38 from [the UCSC repository](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz). On the other hand, your bam file is probably mapped to a slightly differet version of the hg38 assembly or different ordering of chromosomes (a.k.a. non-lexicographic). For instance, if you download a `.bam` file directly from the ENCODE website, it is mapped to a slightly different chromosome ordering of hg38. A recommended way of resolving this issue is to extract a `.fastq` file from your `.bam` file, go back to step 1, and remap it with `bowtie2` using the UCSC hg38 assembly. `bedtools` provides a way to generate a `.fastq` file from your `.bam` file.  
```
bedtools bamtofastq  -i example.bam -fq example.fastq
```  
We will regularly update the control data when a new major version of the genome becomes available; however, covering for all versions with small changes to the existing version is not realistic.
   
### Step 4: Download the AIControl julia script.
The following terminal command will download the AIControl julia script and make it executable. You can also find it within this github repository.
```
wget https://github.com/hiranumn/AIControl.jl/raw/master/aicontrolScript.jl
```
Please also place the downloaded control data files to the same folder, or otherwise specify their location with `--ctrlfolder` option.    
### Step 5: Run AIControl. 
The terminal command below will run AIControl. 
```
julia aicontrolScript.jl example.bam.sorted --ctrlfolder=. --name=test
```
Do `julia aicontrolScript.jl --help` or `-h` for help.

We support the following flags. 

- `--dup`: using duplicate reads \[default:false\]
- `--reduced`: using subsampled control datasets \[default:false\]
- `--ctrlfolder=[path]`: path to a control folder \[default:./data\]
- `--name=[string]`: prefix for output files \[default:bamfile_prefix\]
- `--p=[float]`: pvalue threshold \[default:0.15\]

If you would like to use the `--dup` or `--reduced` options, please download appropriate versions of compressed control data indicated with `.dup` or `.reduced`.

## Simple trouble shooting
Make sure that:
- You are using Julia 1.0.
- You downloaded necessary control files for `--reduced` or `--dup` if you are running with those flags.
- You sorted the input bam files according to the UCSC hg38 assembly as specified in Step 1 (and 3.1).

## We have tested our implementation on ...
- macOS Sierra (2.5GHz Intel Core i5 & 8GB RAM)
- Ubuntu 18.04 
- Windows 8.0

If you have any question, please e-mail to hiranumn at cs dot washington dot edu.
