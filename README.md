# README

cisVar is a pipeline to estimate pre-ChIP frequencies from pooled post-ChIP
frequencies via a regression on the genotypes of the individuals in the pool.
To use this code you must have a mapped BAM file for each pool (a single pool is
fine) as well as genotype info for all individuals in the pool in VCF format.
Additional information can be provided also, see below for information.

For an overview of this method, or to cite this work, please see the following:

Tehranchi, Ashley K., Marsha Myrthil, Trevor Martin, Brian L. Hie, David Golan,
and Hunter B. Fraser. “[Pooled ChIP-Seq Links Variation in Transcription Factor
Binding to Complex Disease Risk](https://doi.org/10.1016/j.cell.2016.03.041).”
Cell 165, no. 3 (April 21, 2016): 730–41.

This code was written by Ashley Tehranchi with minor modifications by Mike
Dacre. It was produced at Stanford University, but is released here under the
MIT license.

Current version: 2.0.0b2

## Overview

`cisVar.py` is the pipeline written in python3 and uses `regression_qtls.R`
(`regression_qtls.R` is called inside of `cisVar.py` so you will not run it
directly).

In addition, the code at `scripts/combine.py` can be used to create combined and
merged (per-SNP) pandas dataframes from the outputs of `cisVar.py` for multiple
samples.

A complete example Snakemake pipeline is provided in `pipeline`, you should be
able to run it by modifying only the `cisvar_config.json` file, see the
Snakemake section below.

Right now the default minor allele frequency filter is `0.1>MAF<0.99`.
To change these and other regression constants, edit the `regressions_qtls.R`
script.

## Installation

Install via PyPI:

```shell
pip install cisVar
```

Or install from github:

```shell
pip install https://github.com/TheFraserLab/cisVar/tarball/master
```

Alternatively, you can just clone the repo and use it directly from there.

```shell
git clone https://github.com/TheFraserLab/cisVar.git
```

It should work with python 2 or 3, but python 3 is highly recommended.


### Prereqs

Requires `samtools v1.9` and `bedtools v2.26.0` as well as the following python
modules (installed automatically by pip):

`pandas`, `numpy`, `psutil`, `snakemake`, `wget`

Additionally requires that `R` with `Rscript` are installed, with the `ggplot2`
module installed.

## Usage

To run this code on your own data you need:

- VCF(s) with genotypes and ref/alt alleles for every SNP you wish to test
- A mapped BAM file, ideally this will make use of
    [Hornet](https://github.com/TheFraserLab/Hornet) to remove reference bias.

In addition, the following data can be helpful:

- A list of individuals to filter genotypes with (one newline separated file of
    individuals per pool)
- A file with ref and alt alleles for your SNPs of interest, to modify those in
    the genotype VCF files (BED/VCF/txt)
- A file to limit the SNPs to consider to a subset of those in the VCF
    (BED/VCF/txt)

Example pipeline:

```shell
cisVar.py prep -F SampleName -i individuals.txt.gz --chrom-format chr /path/to/geno/*vcf.gz

cisVar.py mpileup -F SampleName -f fastaFile -B sortedBam

cisVar.py post -F SampleName -r readDepth -a allelesFile

cisVar.py geno -F SampleName -r readDepth -i individualsFile -g genotypesFile

cisVar.py qtls -F SampleName -r readDepth -n numberIndividuals

cisVar.py tidy -F SampleName -r readDepth

scripts/combine.py {sample}.readDepth.regression.pd
```

A readDepth of 20 is generally optimal, the sample name can be whatever you
want, but should be unique per sample. `{sample}` is a placeholder to allow any
number of samples to be combined in the last step.


### Snakemake

The above pipeline can be automated with
[Snakemake](https://snakemake.readthedocs.io/en/stable/).

To use, install cisVar, navigate to the root of your project, and run `cisVar.py
get_snake` to copy the Snakefile and config file over . Then edit the
`cisvar_config.json` file to match your needs.

You will also need to edit the `Snakefile` to set the `script_prep` string to
match what is needed by your system.

The following are the config options for that file:

| Option       | Description                                                                                                                                                                                                   |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| name         | A general name for this run, file prefix will be <name>.<sample>.<read_depth>                                                                                                                                 |
| sample_name  | The name of the sample, default is population. Used only in the combination of multiple samples.                                                                                                              |
| samples      | A list of samples, can just be a list, or a dictionary of {sample:group}, the 'group' in this case allows the use of the same genotype files for multiple samples, can also be a path to a separate json file |
| read_depth   | An integer depth to require for each SNP to be considered                                                                                                                                                     |
| max_cores    | Used only when parsing VCFs, if you have multiple VCF files (e.g. per chromosome), they will be parsed in parallel up to this many cores (or max avaialable on machine)                                       |
| sort_vcfs    | Either 1 or 0, if 1 assumes that VCF files contain a `chr#` string in the file name, and sorts the order of files to be chr1->22,X,Y,MT. Don't use if your VCFs don't have `chr#` in the name                 |
| chrom_format | 'chr', 'num', 'ignore': Force format of chromosome name to be `chr#` or `#`. This ensures that all input files have the same format. Use ignore to do nothing.                                                |
| bams         | A path to the mapped BAM files, must contain the `{sample}` string (unless you only have one bam), e.g. `/path/to/{sample}.sorted.bam`, `{sample}` must be in samples                                         |
| cisVar       | Path to the cisVar repository                                                                                                                                                                                 |
| vcfs         | Can be a single path (for one vcf), a list of vcfs, or a glob string (e.g. `/path/to/vcfs/*.vcf.comm.gz`)                                                                                                     |
| genome_fa    | Path to a FastA file of the genome you mapped to, single file only.                                                                                                                                           |
| inds         | Optional: used to filter VCFs so that the genotype files contain only the individuals in the sample, e.g. `/path/to/inds/{sample}.ind.txt.gz`. Newline separated file of individuals.                         |
| locs         | Optional: a BED/VCF/text file of SNP locations to consider, used to limit the total to be a subset of the genotype file.                                                                                      |
| alleles      | Optional: a BED/VCF/text file of alternate ref/alt alleles. Must be a subset of the genotype VCFs. If there is an entry in this file, it's ref/alt alleles will be used instead of those in the genotype file |

Note the last three files are optional, also if `samples` is a dict, then the
value will be used in place of the sample. For example, if you have two samples
for the same population that are `yri1` and `yri2`, but they both use the same
genotype file `yri.geno.vcf`, you can make samples `{'yri1': 'yri', 'yri2':
'yri'}` and then `yri` will be used to pick the ind, loc, and allele files


#### Cluster

To run on a cluster, run `cisVar.py get_snake` with `-x` and edit the
`cluster.json` file to match your cluster environment, then run e.g.:

```shell
snakemake -j 100 --cluster-config cluster.json --cluster sbatch -n {threads} \
-t {params.time} --mem={resources.mem_mb} -p {cluster.queue} \
-o {cluster.out} -e {cluster.err}
```

or

```shell
snakemake -j 100 --cluster-config cluster.json --cluster qsub \
-l nodes=1:ppn={threads} -l walltime={params.time} -l mem={resources.mem_mb} \
-o {cluster.out} -e {cluster.err}
```

To set the maximum allowed memory per job, add the argument
`--resources mem_mb=32000`.

To also combine files, add `combine` to the end of the command.


## Script help

Below are help options available on the command line for cisVar, all these steps
are run by the above snakemake pipeline.

### cisVar.py

```
usage: cisVar.py [-h] {prep,mpileup,post,geno,qtls,tidy,get_snake} ...

cisVar: Find cis QTLs based on an experimental selection method

Ashley Tehranchi <tehranchi576@gmail.com>

Stanford University

Version: 2.0.0b1
Created: 2015-12-12
Updated: 2018-05-16

Example usage:
cisVar prep -F test_new -i individuals.txt.gz --chrom-format chr
cisVar mpileup -F <SampleName> -f <fastaFile> -B <sortedBam>
cisVar post -F <SampleName> -r <readDepth>
cisVar geno -F <SampleName> -r <readDepth> -i <individualsFile>
cisVar qtls -F <SampleName> -r <readDepth> -n <numberIndividuals>
cisVar tidy -F <SampleName> -r <readDepth> -p out.dataframe.pandas -t out.dataframe.txt

Note:
The qtls regression step will use approximately 32GB of memory on an averaged-
sized dataset.

The geno step will use approximately 20GB of memory on the same dataset.

positional arguments:
  {mpileup,post,geno,qtls}
    prep                Prepare genotype files
    mpileup (m)         Run mpileup
    post (p)            Run POST frequency calculation
    geno (g)            Munge genotypes to prepare for regression
    qtls (q, regression, r)
                        Run the regression
    tidy (t)            Tidy up regression, call open/clsoed

optional arguments:
  -h, --help            show this help message and exit
```

#### prep

This step converts VCFs into genotype and individual files that can be used by
the pipeline.

```
usage: cisVar.py prep [-h] [-F PREFIX_NAME] [-r TRIAL_DEPTHS] [-i ALL_INDS]
                      [-l LIMIT_FILE] [-a ALLELE_FILE]
                      [--chrom-format {chr,num,ignore}] [--include-indels]
                      [-c CORES]
                      vcf_files [vcf_files ...]

Prepare genotype files

optional arguments:
  -h, --help            show this help message and exit

Run Options:
  -F PREFIX_NAME, --SampleName PREFIX_NAME
                        sample/population name (default: cis_var)
  -r TRIAL_DEPTHS, --readDepth TRIAL_DEPTHS
                        minimum read depth per variant (default: 20)

Prep Options:
  -i ALL_INDS, --all-inds ALL_INDS
                        File of individuals in all groups, one per line
  -l LIMIT_FILE, --limit-file LIMIT_FILE
                        BED/VCF/txt file of SNPs to consider
  -a ALLELE_FILE, --allele-file ALLELE_FILE
                        BED/VCF/txt file of alleles to override VCF allels
                        (subset of vcf)
  --chrom-format {chr,num,ignore}
                        chr: make format "chr#", num: make format "#", ignore:
                        do nothing (default: ignore)
  --include-indels      Do not skip indels
  -c CORES, --cores CORES
                        Number of cores to use (default: all)
  vcf_files             VCF files with genotypes
```

#### mpileup

This is just a simple wrapper for samtools mpileup

```
usage: cisVar.py mpileup [-h] [-F PREFIX_NAME] [-r TRIAL_DEPTHS] -f
                         ALLCHRFASTA -B SORTEDBAM [-p MPILEUPBEDFILE]

Run mpileup

optional arguments:
  -h, --help            show this help message and exit

Run Options:
  -F PREFIX_NAME, --SampleName PREFIX_NAME
                        sample/population name (default: cis_var)
  -r TRIAL_DEPTHS, --readDepth TRIAL_DEPTHS
                        minimum read depth per variant (default: 20)

mpileup Options:
  -f ALLCHRFASTA, --fasta ALLCHRFASTA
                        fasta file with all chromosomes (Required)
  -B SORTEDBAM, --BAMfile SORTEDBAM
                        sorted BAM file (Required)
  -p MPILEUPBEDFILE, --mpileupBEDfile MPILEUPBEDFILE
                        BED to use instead of the BED generated in the prep
                        phase (Do not use if possible, use prep with limit
                        instead)
```

#### post

This step actually calculates the POST-frequencies for the data.

```
usage: cisVar.py post [-h] [-F PREFIX_NAME] [-r TRIAL_DEPTHS] [-a GENOSFILE]

Run POST frequency calculation

optional arguments:
  -h, --help            show this help message and exit

Run Options:
  -F PREFIX_NAME, --SampleName PREFIX_NAME
                        sample/population name (default: cis_var)
  -r TRIAL_DEPTHS, --readDepth TRIAL_DEPTHS
                        minimum read depth per variant (default: 20)

POST Options (Deprecated):
  -a GENOSFILE, --allelesFile GENOSFILE
                        The genotypes file, (Optional, default is file created
                        in prep)
```

#### geno

This step converts the genotype file made in the prep step into a matrix that
can be used in the regression. It is important that this genotype file is
perfectly sorted to match the outputs of the POST step.

```
usage: cisVar.py geno [-h] [-F PREFIX_NAME] [-r TRIAL_DEPTHS] [-g GENOSFILE]
                      [-i INDIVIDUALSLIST]

Munge genotypes to prepare for regression

optional arguments:
  -h, --help            show this help message and exit

Run Options:
  -F PREFIX_NAME, --SampleName PREFIX_NAME
                        sample/population name (default: cis_var)
  -r TRIAL_DEPTHS, --readDepth TRIAL_DEPTHS
                        minimum read depth per variant (default: 20)

Genotype Options:
  -g GENOSFILE, --genoFile GENOSFILE
                        The genotypes file, (Optional, default is file created
                        in prep)
  -i INDIVIDUALSLIST, --individualsFile INDIVIDUALSLIST
                        list of individuals matching genotype matrix; one indv
                        per line
```

#### qtls

This is the actual regression step, it makes sure all the files are in the right
place and then calls `regression_qtls.R` to do the actual regression.

```
usage: cisVar.py qtls [-h] [-F PREFIX_NAME] [-r TRIAL_DEPTHS] [-n NUMINDV]

Run the regression

optional arguments:
  -h, --help            show this help message and exit

Run Options:
  -F PREFIX_NAME, --SampleName PREFIX_NAME
                        sample/population name (default: cis_var)
  -r TRIAL_DEPTHS, --readDepth TRIAL_DEPTHS
                        minimum read depth per variant (default: 20)

Regression Options:
  -n NUMINDV, --numberIndividuals NUMINDV
                        The number of individuals in the pool (if omitted,
                        calculated from genotype file length)
```

The regression produces z-scores and p-values, and additionally writes
coefficients and some simple summary plots in separate files.

#### tidy

This step calls the open/closed alleles and produces a final integrated file
with all available data as both a tad-delimited file and as a pandas dataframe.

```
usage: cisVar.py tidy [-h] [-F PREFIX_NAME] [-r TRIAL_DEPTHS] [-b BEDFILE]
                      [-t TEXTFILE] [-p PANDASFILE]

Tidy up regression, call open/closed

optional arguments:
  -h, --help            show this help message and exit

Run Options:
  -F PREFIX_NAME, --SampleName PREFIX_NAME
                        sample/population name (default: cis_var)
  -r TRIAL_DEPTHS, --readDepth TRIAL_DEPTHS
                        minimum read depth per variant (default: 20)

inputs:
  -b BEDFILE, --bedfile BEDFILE
                        BED file to extract rsIDs from (optional)

outputs:
  -t TEXTFILE, --textfile TEXTFILE
                        Parsed output
  -p PANDASFILE, --pandasfile PANDASFILE
                        Parsed dataframe
```

#### get_snake

This option just downloads the Snakefile and config files from this repo, for
easy access when code is installed via pip.

```
usage: cisVar.py get_snake [-h] [-x]

Download Snakefile and config to current dir

optional arguments:
  -h, --help   show this help message and exit
  -x, --extra  Get additional sample and cluster configs
```

### combine.py

This script is separate and is in the `scripts` folder. It takes a search string
as an input and produces both combined and merged DataFrames. The combined
dataframe is just all dataframes combined in order with sample data added as a
column and to the index. The merged dataframe is a collapsed dataframe that has
one entry per SNP with p-values combined using Fisher's method and supporting
population data. It also includes information on the level of support for the
open and closed calls.

The search string should match your prefix and depth from the main pipeline. For
example, if you used a name of 'cis_var' plus a sample name (the variable part)
of e.g. CEU and YRI, and a read depth of 20, your search string would be:
`cis_var.{sample}.20.regression.pd`.

The script will write `cis_var.combined.20.regression.pd` and
`cis_var.merged.20.regression.pd`.

```
usage: combine.py [-h] [-c COLUMN_NAME] [--no-merge] search_str

Combine a bunch of cisVar pandas files by sample (e.g. population).

Requires a search string such as prefix.{sample}.regression.pd.

Writes
------
prefix.combined.regression.pd
    A simple merger of all DataFrames
prefix.merged.regression.pd
    A per-snp merger based on p-value

positional arguments:
  search_str            e.g. name.{sample}.regression.pd, used to find files

optional arguments:
  -h, --help            show this help message and exit
  -c COLUMN_NAME, --column-name COLUMN_NAME
                        Name for combine column, e.g. population
  --no-merge            Produce only a combined dataframe, not a merged
                        dataframe. merging can add half an hour over
                        combination, which takes seconds
```

### plot

There is an additional script in `scripts` called `plot_fun.R` that takes a
single argument—the output of the regression step (e.g.
`cis_var.YRI.20.totals.txt`) and creates a simple density pre-freq vs post freq
plot.
