# Global 5mC DNA Methylation Profiling using Whole Genome Bisulfite Sequencing (WGBS)

__Jimmy Breen (jimmy.breen@sahmri.com)__  
Head of Bioinformatics, South Australian Genomics Centre (SAGC)  
Group Leader, Robinson Research Institute  
Senior Research Associate, Adelaide Medical School, University of Adelaide  

## Command-line linux

Today I will be going through WGBS analysis via a Linux command-line interface, using a terminal-style program located on your assigned virtual machine (VM). 
This system might seem unfamiliar to many of you, but it is an essential part of most Bioinformatics techniques.
To avoid the tutorial from being too long, we won't go into the basics of the command-line. 
For those completely new to the system, I would point you towards the [University of Adelaide Bioinformatics Hub's Introduction to Linux course](https://github.com/UofABioinformaticsHub/BASH-Intro/blob/master/notes/1_bash.md), which has detailed outline of all the common linux commands that are used in Bioinformatics work.  

For today we will be using these basic commands:
- `cd`: Change directory
- `ls`: List the files in a directory
- `mv` & `cp`: Move/rename and copy files
- The pipe `|`: A device for tying multiple commands together

## Whole genome bisulfite sequencing (WGBS)

From [Enseqlopedia](http://enseqlopedia.com/wiki-entry/dna-sequencing-methods/epigenetics/bs-seqbisulfite-seqwgbs/):  

BS-Seq/Bisulfite-seq or WGBS is a well-established protocol to detect methylated cytosines in gDNA (Feil et al., 1994). In this method, gDNA is treated with sodium bisulfite and then sequenced, providing single-base resolution of methylated cytosines in the genome. Upon bisulfite treatment, unmethylated cytosines are deaminated to uracils which, upon sequencing, are converted to thymidines. Simultaneously, methylated cytosines resist deamination and are read as cytosines. The location of the methylated cytosines can then be determined by comparing treated and untreated sequences. Bisulfite treatment of DNA converts unmethylated cytosines to thymidines, leading to reduced sequence complexity. Very accurate deep sequencing serves to mitigate this loss of complexity.

_Advantages:_  
- CpG and non-CpG methylation throughout the genome is covered at single-base resolution
- Covers 5mCs in dense and less dense repeat regions

_Disadvantages:_
- Bisulfite converts unmethylated cytosines to thymidines, reducing sequence complexity, which can make it difficult to create alignments
- SNPs where a cytosine is converted to thymidine will be missed upon bisulfite conversion
- Bisulfite conversion does not distinguish between 5mC and 5hmC

## Analysis Workflow



## Available community-driven analysis workflows

Over the past couple of years, there has been a greater push in the Bioinformatics community to develop more community driven analysis standards.
These standards enable new researchers to get straight into 

https://www.encodeproject.org/pipelines/ENCPL985BLO/

https://github.com/nf-core/methylseq

https://github.com/dohlee/snakemake-bismark-methyldackel