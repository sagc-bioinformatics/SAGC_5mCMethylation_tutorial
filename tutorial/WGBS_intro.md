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

For today we will be using these basic commands (care of Fides D Lay: UCLA):

- Where am I?: `pwd`
- Home directory: `~/`
- Change directory: `cd ~/data`
- Move up one level: `cd ..`
- List files in folder: `ls`
- Look at a file: `less fileName`
- Copy a file: `cp ~/data/file ~/otherdir/`
- Delete a file: `rm fileName`
- Delete a directory: `rmdir ~/dirName/`
- Move a file: `mv ~/data/file ~/otherdir/file`
- Secure copy: `scp user@host1:dir/file user@host2:dir/file`
- Compress a file: `gzip –c file > file.gz`
- Uncompress a file: `gunzip file.gz`
- Print a compressed file: `zcat file.gz`
- Make a new folder: `mkdir data2`
- Current directory: `./`
- Count lines in a file: `wc –l fileName`

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

The basic workflow required to go from sequencing to methylation calling is fairly standard across WGBS datasets, with slight variations coming from the different library preparation strategies and/or alignment algorithms

1. Sequencing QC: Basic quality control of raw high-throughput sequencing data
2. Adapter & Quality Trimming: Using information gathered in the QC step, sequencing adapters and other known biases can be identified and removed. Poor quality or biases sequence composition has the ability to impact downstream analyses
3. Sequence mapping: This computational intensive step can take significant amounts of time given the complexities of mapping bisulfite converted sequences. 
4. Post-mapping quality control: After mapping, base composition biases may be identified at ends of fragments, requiring additional investigation. Duplicates may also need to be removed in order to identify accurate methylation calls.
5. Methylation calling: Takes WGBS mapping data and calls a methylation proportion at each site (CpG, CHG or CHH).
6. Differentially Methylated Region (DMR) detection: Much like identifying differential expression from RNA-seq data, this step uses differential statistics to identify regions of the genome that have changed methylation status in accordance with the experimental design context.
7. Data intergation: This final step is where methylation data, whether it be DMRs or global methylation patterns, is combined with other genome annotation in order to give context to results.


## Available community-driven analysis workflows

Over the past couple of years, there has been a greater push in the Bioinformatics community to develop more community driven analysis standards.
These standards and analysis protocols enable new researchers to get straight into the analysis and enable strong reproducibility across important datasets.
Much of the work described in this tutorial is based on community standards in the Epigenetics field, driven by our experience at the SAGC, international sequencing projects such as ENCODE, NIH Epigenomics Roadmap and European Union Blueprint consortium.

If you are new to the field, I would highly recommend implementing standardised datasets, relying on a workflow language system such as `nextflow`, `snakemake` or Common Workflow Language (CWL).
Most of these have clear standard inputs and outputs and are extensively documented.


Below are a number of recommendations for DNA methylation:

- [ENCODE's DNA methylation pipeline](https://www.encodeproject.org/pipelines/ENCPL985BLO/)
- [Nextflow's methylseq nf-core workflow](https://github.com/nf-core/methylseq)
- [Snakemake DNA methylation example workflow](https://github.com/dohlee/snakemake-bismark-methyldackel)

---

NEXT: Lets get straight into it with [Sequencing QC](tutorial/BS_quality_control.md)