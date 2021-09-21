# Post mapping quality control & methylation calling

After we're done with mapping our data, we need to process our BAMs. 
We ned to do the following tasks:

1. Deduplicate the data to remove clonal deduplicates
2. Sort based on reference genome coordinate
3. Make methylation call files

## Deduplication

Deduplication is an important aspect of genomics applications where a precise quanitification is required at base-pair level.
In WGBS analysis, we need to accurately determine the proportion of methylated to unmethylated bases, meaning that any uneven amplification will impact the methylation value identified during calling.
This is a similar situation to variant calling workflows, that aim to identify the proportion of reference and alternate alleles at each individual site (i.e. Allele-Frequency).

Without a unique tag during library preparation, such as a Unique Molecular Identifier (UMI), it is almost impossible to know what is a PCR duplicate.
However the most common definition of a duplicate for sequencing data is a read that start and ends at identical positions. 
Bismark has its own deduplication command which we will use for today. 

In their own documentation: 
	_"This command will deduplicate the Bismark alignment BAM file and remove all reads but one which align to the the very same position and in the same orientation. This step is recommended for whole-genome bisulfite samples, but should not be used for reduced representation libraries such as RRBS, amplicon or target enrichment libraries."_

So lets deduplicate!

	# SRR534177
	deduplicate_bismark --bam SRR534177_colWT_trimmed_bismark_bt2.bam

	# SRR534239
	deduplicate_bismark --bam SRR534239_met1_trimmed_bismark_bt2.bam


Each produces a BAM file and a summary output that we can show which can be viewed to see how many duplicates were removed:

	# List the files in the directory
	$ ls -l

	# View summary report
	$ less SRR534177_colWT_trimmed_bismark_bt2.deduplication_report.txt


To continue onto the methylation calling stage, we will need to resort the BAM file.
You have probably noticed that in their current states, each raw and deduplicated BAM file are ordered by the first column, which is the read name.

	$ samtools view SRR534177_colWT_trimmed_bismark_bt2.bam | head 

However to call each reference base, we will need to sort the BAM file by coordinate, so sequences can be easily scanned and summarised on a per base level.
`samtools` has a sort command which can then be used to create a new coordinate sorted file.
We will also need to index the file after sorting which will be used later on.

	# Sort SRR534177 and index
	$ samtools sort -o SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted.bam \
		SRR534177_colWT_trimmed_bismark_bt2.deduplicated.bam
	$ samtools index SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted.bam
	

	# Sort SRR534239 and index
	$ samtools sort -o SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted.bam \
		SRR534239_met1_trimmed_bismark_bt2.deduplicated.bam
	$ samtools index SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted.bam



**Questions**  
1. How many duplicates were removed in each sample?

---

## Methylation calling

"Methylation calling" is a process of quantifying the proportion of methylation at a specific 5mC context.
At each individual Cytosine site across the reference sequence, converted and un-converted cytosine bases are tallied and quantified to a specific proportion of methylation (i.e. methylated C/(methylated C + non-methylated C)).
Generally most mammalian BSseq sequencing project can get away with only identifying the proportion of methylated cytosines in CpG conetxts, but in other species 

Bismark does come with methylation call extraction scripts, however I have found that the program `MethylDacyl` (formerly `PileOMeth`) is a very quick caller that is multi-threaded and easy to use.

Before we start to call methylated sites, there are a number of considerations to take into account

1. Methylation bias: While the probability of identifying 5mC sites across a sequenced read should be constant, biases can be introduced by library preparation to impact that. While we do not see that in our FastQC report (and we wont be delving into this for our data), `MethylDacyl` does have a nice `MethylDacyl mbias` subcommand that you can use to quanitfy those impacts. These can then be applied to the command to extract accurate methylation calls from the alignment file
2. Coverage considerations: If you're data is low-coverage, what should be your coverage cut-off for calling a methylated base. For example, if you have only 5x coverage over a site, 2 being methylated and 3 being unmethylated, your 40% methylation figure is probably far less accurate than a 20x coverage calculation. We generally call based on 10x coverage, unless looking at single cell or low-coverage inputs
3. Sequencing & alignment quality: The alignment file will carry through information about the individual base-quality (taken from the original FASTQ file). While these are generally trimmed by `trim-galore` at the end of reads, sometimes they can be found within the read that is aligned. The uniqueness of the alignment also can impact the quality of the methylation call, as the more repetitive the read is the less likely it is to be correctly placed.

Ok lets call some methylation!

The format for the extraction command is as follows:

	MethylDackel extract [OPTIONS] <ref.fa> <sorted_alignments.bam>

We need our reference genome sequence `TAIR10_chr1_cp.fa.gz` (NOT the `bismark` formatted genome directory) and our coordinate sorted, deduplicated alignment file `SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted.bam`.
I have also added some commands which we can use to extract the best results based on sequence quality, mapping quality and coverage recommendations:

	$ MethylDackel extract -d 10 -q 10 -p 5 Athal/TAIR10_chr1_cp.fa.gz SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted.bam

	$ MethylDackel extract -d 10 -q 10 -p 5 Athal/TAIR10_chr1_cp.fa.gz SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted.bam

This produces a [bedGraph file or tab-delimiated file](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) for each CpG/CG site in the genome that satisfies the additional options on quality and coverage.
They are commonly used in genome browsers to display interval sets, such as CG coordinates in this case.

Lets have a look at it:

	# View both files
	$ head *_CpG.bedGraph
	==> SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted_CpG.bedGraph <==
	track type="bedGraph" description="SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted CpG methylation levels"
	1	115	116	100	11	0
	1	911	912	90	10	1
	1	948	949	82	14	3
	1	1438	1439	0	0	11
	1	5448	5449	63	7	4
	1	5458	5459	66	10	5
	1	5464	5465	73	11	4
	1	5487	5488	30	3	7
	1	5848	5849	0	0	10

	==> SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted_CpG.bedGraph <==
	track type="bedGraph" description="SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted CpG methylation levels"
	1	4107	4108	0	0	12
	1	5458	5459	0	0	10
	1	8642	8643	0	0	11
	1	8658	8659	0	0	13
	1	8661	8662	0	0	11
	1	12227	12228	10	1	9
	1	12302	12303	0	0	10
	1	12394	12395	0	0	10
	1	12726	12727	0	0	11

As you can see, the first line is a [track definition line](https://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK) containing information relating to what the data contains.
This enables genome browsers, most noteably the [UCSC Genome Browser](https://genome.ucsc.edu/index.html), to show the user 

Ignoring the first line, what do each of these columns represent?
- Column 1 is Chromosome
- Column 2 is the start position of the CG site
- Column 3 is the end position of the CG site
- Column 4 is percentage methylation at each CG site
- Column 5 is number of methylated bases at that CG site
- Column 6 is number of unmethylated bases at that CG site

So the actual proportion of methylation, basing the quality and coverage thresholds outlined within the command, is contained within the fourth column.

We can slo count the number of CpG sites within the samples using `wc -l` (i.e. count the lines in the file -1)

	$ wc -l *_CpG.bedGraph
	92955 SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted_CpG.bedGraph
	43195 SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted_CpG.bedGraph
	136150 total

So the colWT sample has almost double the number of CpG sites than the met1 sample, with a similar number of input reads.
What about per chromosome? 
Lets look at those figures:

	$ cut -f1 SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted_CpG.bedGraph | \
                sort | uniq -c 
	88028 1
	4926 chloroplast
		1 track type="bedGraph" description="SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted CpG methylation levels"

	$ cut -f1 SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted_CpG.bedGraph | \
                sort | uniq -c
	38188 1
	5006 chloroplast
		1 track type="bedGraph" description="SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted CpG methylation levels"

What about CpG sites that contain more than 50% methylation?
We can use the linux command awk to filter those sites, by filtering sites with >50% methylation in column 4

	$ awk '$4>50' SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted_CpG.bedGraph

	$ awk '$4>50' SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted_CpG.bedGraph



Extending this approach, lets output CHG and CHH metrics along with CpG by adding the `--CHG --CHH` flags to the same command.


	$ MethylDackel extract -d 10 -q 10 -p 5 --CHG --CHH Athal/TAIR10_chr1_cp.fa.gz  \
	        SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted.bam

	$ MethylDackel extract -d 10 -q 10 -p 5 --CHG --CHH Athal/TAIR10_chr1_cp.fa.gz \
		SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted.bam

Non-CpG methylation is an [important component to many organisms (see Wikipedia)](https://en.wikipedia.org/wiki/DNA_methylation):

	In plants and other organisms, DNA methylation is found in three different sequence contexts: CG (or CpG), CHG or CHH (where H correspond to A, T or C). In mammals however, DNA methylation is almost exclusively found in CpG dinucleotides, with the cytosines on both strands being usually methylated. Non-CpG methylation can however be observed in embryonic stem cells, and has also been indicated in neural development. Furthermore, non-CpG methylation has also been observed in hematopoietic progenitor cells, and it occurred mainly in a CpApC sequence context.

In _Arabdopsis_ and other plants, different genes are responsibile for the deposition and maintainence of non-CpG DNA methylation types, making them more important biologically compared to mammalian species.

Lets have a quick look at the CHG and CHH files.

	$ head *_CHG.bedGraph
	==> SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted_CHG.bedGraph <==
	track type="bedGraph" description="SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted CHG methylation levels"
	1	116	117	33	4	8
	1	928	929	38	7	11
	1	1432	1433	0	0	12
	1	5435	5436	0	0	10
	1	5450	5451	0	0	12
	1	5478	5479	9	1	10
	1	5793	5794	0	0	18
	1	6082	6083	0	0	10
	1	8444	8445	0	0	18

	==> SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted_CHG.bedGraph <==
	track type="bedGraph" description="SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted CHG methylation levels"
	1	4106	4107	0	0	11
	1	4115	4116	0	0	10
	1	4132	4133	0	0	15
	1	8444	8445	8	1	11
	1	8648	8649	0	0	11
	1	12245	12246	7	1	12
	1	16729	16730	0	0	10
	1	23368	23369	0	0	11
	1	24051	24052	0	0	12


	$ head *_CHH.bedGraph
	==> SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted_CHH.bedGraph <==
	track type="bedGraph" description="SRR534177_colWT_trimmed_bismark_bt2.deduplicated.sorted CHH methylation levels"
	1	238	239	0	0	10
	1	240	241	10	1	9
	1	241	242	10	1	9
	1	1409	1410	0	0	10
	1	1715	1716	0	0	13
	1	1720	1721	0	0	14
	1	1721	1722	0	0	14
	1	4125	4126	0	0	10
	1	4131	4132	0	0	10

	==> SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted_CHH.bedGraph <==
	track type="bedGraph" description="SRR534239_met1_trimmed_bismark_bt2.deduplicated.sorted CHH methylation levels"
	1	4102	4103	0	0	11
	1	4125	4126	0	0	15
	1	4131	4132	0	0	15
	1	4155	4156	0	0	11
	1	5161	5162	0	0	10
	1	5162	5163	0	0	10
	1	8264	8265	0	0	10
	1	8650	8651	0	0	11
	1	12042	12043	0	0	10



## Questions

1. How many >50% CpG sites were identified in colWT and met1 samples?
2. How many >50% CpG sites were identified in chr1 and chloroplast sequences in colWT?
3. Based on the CpG profiles shown for the colWT and met1 samples, what can you say about the genes the biological impact of _met1_ on 5mC DNA methylation? 
4. What is the impact of this gene on CHG and CHH contexts?
5. Based on the CpG profiles shown for colWT, what can you say about the amount of 5mC DNA methylation in organelles vs autosomal sequences? 


---

LASTLY: To finish off we're going to skip the DMR step and get stuck into [data integration](05_data_integration.md)
