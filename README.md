# README

`cisVar.py` is the pipeline written in python3 and uses `regression_qtls.R`.

`regression_qtls.R` is called inside of `cisVar.py` so you will not run it
directly.

Right now the default minor allele frequency filter is `0.1>MAF<0.99`.
To change these and other constants edit the `regressions_qtls.R` script.

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
```

Requires `samtools v1.9<` and `bedtools v2.26.0<`.

File formats are described in the script help section

## Script help

```
-h	help

-g	The genotypes file is a tab separated file, sorted by individuals(columns) and by snps (rows).

	EX:
		chr   position    ref alt NA18486 NA18489 ...
		chr1  10583   g   a   0   1 ...
		chr1  10611   C   G   0   0 ...
		chr1  13302   c   t   1   0 …

-i	The individualsFile is just a list of the individuals (one per line) in your pool that are to be extracted from a larger list of individuals.

	EX:
		NA18486
		NA18489
		…
-a	The alleles file is a tab-sep file with no header (chr#, position, Ref_allele, Alt_allele). This file should match the length of the pileup file.

	EX:
		chr1	10505	A	T
		chr1	10506	C	G
		chr1	10511	G	A

-f	fasta file with all chromosomes concatenated and sorted.

-p	The mpileup.bed file is a tab-sep standard bed file with no header (where the SNP location is the third column and the second is SNP-1.
    	EX:

	    chr1	10504	10505
	    chr1	10505	10506
	    chr1	10510	10511

-B	sortedBam is the sorted without the -n flag i.e; samtools sort [in.bam]
```

## Additional Scripts

The `scripts` folder includes a few additional scripts that may be useful

- `convert_vcf2bed.sh` — Used bedops to convert vcfs to bed files (which can then be concatenated to make the required bed
- `make_alleles.py` — Makes the alleles file described above from the mpileup file
- `vcf_to_indi_and_geno.py` — Make the individual and genotype files from VCFs
- `post_process_regression.py` — Use pandas to process the output of the qtls step of the pipeline to add open/closed alleles and rsIDs
- `rmdup.py` — Remove duplicates from a bam file, picks duplicate to keep randomly, works with paried end data
