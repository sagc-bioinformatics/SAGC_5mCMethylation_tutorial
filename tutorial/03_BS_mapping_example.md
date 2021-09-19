# WGBS Genome Mapping

## Preparing the Reference Genome

Aligning to a reference sequence when you have treated your sample with Sodium bisulfite is no easy task.
Global or "end-to-end" read aligners are generally very susceptible to sequence variants. 
During bisulfite conversion, all cytosines (C) will be converted to a Uracil (U), while 5' methyl cytosine's will be left unconverted as a cytosine (C). 
During amplification the converted Uracil's (U) are read as Thymine's (T).

![bisulfite conversion](https://www.diagenode.com/img/applications/bisulfite.png)  
*www.diagenode.com*

The consequence is that the conversion introduces a significant amount of sequence diversity which impacts genome alignment.
However, there are a number of programs that use existing global read aligners, such as `bwa` and `bowtie(1/2)`, by aligning reads to converted and non-converted reference genomes.

![Bismark alignment strategy](https://d3i71xaburhd42.cloudfront.net/a114acb9c90d29d9611674824b01a007b6b7a115/2-Figure1-1.png)

The mapping procedure and program that we will be using today is the `bismark` WGBS program, which uses the `bowtie2` alignment tool to align reads.
As shown above, reads are aligned to a reference genome that is converted into two different `bowtie2` indexes: 
- a reference sequence that converts all the C's to T's
- a reference sequence that converts all the G's to A's (i.e. C->T on the opposite strand)

With our reduced reference genome of _Arabidopsis_thaliana_, we now need to prepare our genome for alignment.
Move the reference genome sequence (thats in a compressed format "fa.gz") into a new directory and use the `bismark_genome_preparation` command to run all the `bowtie2` indexing.

	# Create directory for bismark WGBS genome
	mkdir Athal

	# Move data file in there
	mv TAIR10_chr1_cp.fa.gz Athal/

	# Run bismark to format our Athaliana genome (TAIR10)
	bismark_genome_preparation --bowtie2 Athal

Now we have all of the index sequences we need to map our WGBS reads to _Arabidopsis_thaliana_.
To demonstrate all the files `bismark_genome_preparation` has made, lets list all of the files found within the new directory.

	$ ls -R Athal

**Questions**  
1. How many `bowtie2` index files are created for each converted sequence?

---

## Mapping

Single end data


	bismark -p 2 --bam --bowtie2 Athalchr1cp SRR534177_colWT_trimmed.fq.gz


## Data output

Lets have a look at our brand new alignment file.
Firstly, lets look at the header:

	$ samtools view -H SRR534177_colWT_trimmed_bismark_bt2.bam
	@HD	VN:1.0	SO:unsorted
	@SQ	SN:1	LN:30427671
	@SQ	SN:chloroplast	LN:154478
	@PG	ID:Bismark	VN:v0.22.1	CL:"bismark -p 8 --bam --bowtie2 Athalchr1cp SRR534177_colWT_trimmed.fq.gz"

Now lets have a look at the rest of the BAM file.
Because there is a lot of lines in this file, lets just look at the top portion:

	$ samtools view  SRR534177_colWT_trimmed_bismark_bt2.bam | head -n 4
	SRR534177.13_SN603_WA034:3:1101:45.50:110.90/1	16	1	10011165	42	40M	*	0	0	AATCGCATATTCACATACATATAATAAAATTTTATAAAAT	GHE?4>HBCHBHBBHFFEEBHAIIIGIIIIHBFBFFHDDD	NM:i:7	MD:Z:0G3A11G6G4G4G1G4	XM:Z:h...............h......h....h....h.h....	XR:Z:CT	XG:Z:GA
	SRR534177.18_SN603_WA034:3:1101:42.80:117.20/1	16	1	17342024	42	40M	*	0	0	TACTCTAAAAATTATACCATTATTTTGTTAAATCCGTACC	FBIIIIJJJJIIGF?9GGEGGFIGGHBGGHGGFHHFD?FD	NM:i:8	MD:Z:1G4G1G6G5G7G0G6G2	XM:Z:.x....x.h......h.....h....H..hh....Z.h..	XR:Z:CT	XG:Z:GA
	SRR534177.34_SN603_WA034:3:1101:51.50:109.90/1	16	1	19411065	42	39M	*	0	0	CAAATATAAATTCCGCACAAATTCCTACAACATTTTCAT	JJJJJJIJJJGGJIIIHJJJJIGIHGDJIIJHHHHHFFF	NM:i:4	MD:Z:9G9G0G8G9	XM:Z:.........h....Z....xh........x.........	XR:Z:CT	XG:Z:GA
	SRR534177.36_SN603_WA034:3:1101:72.10:110.50/1	0	chloroplast	47672	42	37M	*	0	0	ATAGATTTAAGTTATATATTAAAATGATATTGATATT	FFFHHHHHJJJIJJJJJJJJJJJJJJIJJJJJJJIJJ	NM:i:5	MD:Z:5C0C5C16C5C1	XM:Z:.....hh.....h................x.....h.	XR:Z:CT	XG:Z:CT

As you can see, the last fields

	less SRR534177_colWT_trimmed_bismark_bt2_SE_report.txt



## Summary statistics

The amount of data thats produced from the alignment file 

	$ samtools flagstat  SRR534177_colWT_trimmed_bismark_bt2.bam


---

NEXT: We'll now start using these BAM files to [call methylated sites](04_methylation_calling.md)
