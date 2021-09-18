
# Dataset & Quality Control

## Tutorial Data

Our tutorial is only lasting a few hours and because we don't have a significant amount of resources to use, we will be using a subset WGBS sequencing run as today's dataset. 
To further reduce the strains on our resources, we'll be mapping our dataset against the model plant genome _Arabidopsis_thaliana_, which co-incidently was one of the first organisms to be run using whole-genome bisulfite sequencing (WGBS/BS-seq/MethylC-seq) in [Ryan Lister's 2008 _Cell_ paper](http://www.sciencedirect.com/science/article/pii/S0092867408004480). 
_Arabidopsis_thaliana_ is roughly 125Mb in size, making it a easy organism to map to, and there have been some great papers published which sequence the methylomes of a range of DNA methyltransferase (DNMT) mutants.

In this tutorial we're taking two samples from the paper:

    Stroud H, Greenberg MV, Feng S, Bernatavichute YV, Jacobsen SE. Comprehensive analysis of silencing mutants reveals complex regulation of the Arabidopsis methylome. Cell. 2013 Jan 17;152(1-2):352-64. doi: 10.1016/j.cell.2012.10.054. Epub 2013 Jan 11. Erratum in: Cell. 2015 Jun 18;161(7):1697-8. PMID: 23313553; PMCID: PMC3597350.

![](https://els-jbs-prod-cdn.jbs.elsevierhealth.com/cms/attachment/2a889f3d-454d-433a-9ac2-2b020fa86705/fx1.jpg)

For plant researchers, or anyone interested in the mechanisms of DNA methylation maintenance for that matter, this paper is a great resource.
There is a significant amount of sequencing data produced from _Arabidopsis thaliana_ lines that contain silencing mutations in genes involved in DNA methylation-derived gene regulation. 
All these sequencing datasets are publicly available and can be used to study how methyltransferases deposit DNA methylation markers and what genes are regulated by each.
Most of these mutants have multiple replicates, but for today we are going to be taking a subset of two WGBS libraries sequenced from a [cytosine-DNA-methyltransferase gene mutant (_met1_)](https://pubmed.ncbi.nlm.nih.gov/12663548/) and the standard Col or Columbia ecotype (_colWT_).
_MET1_ is the gene responsible for the maintenance of CpG methylation in Plants, so silencing this gene has a major impact on DNA methylation throught the organism.

We will need to rename these FASTQ files to make them a little more easy to 

	# Run fastQC
	fastqc -t 2 SRR534239_met1.fastq.gz SRR534177_colWT.fastq.gz

The two samples are the Columbia wild-type reference and _met1_, a mutation of the CG methyltransferase MET1 resulting in elimination of CG methylation throughout the genome.

The reference genome of _Arabidopsis_thaliana_ is in the /data directory, so lets format that for bisulfite sequence mapping. We need to make sure that bismark recognises that you want to use bowtie2 to do the mapping.


	# Create directory for bismark BS-seq genome
	mkdir Athal

	# Move data file in there
	mv data/TAIR10.fa Athal/

	# Run bismark to format our Athaliana genome (TAIR10)
	bismark_genome_preparation --bowtie2 Athal

Now our genome 

We may need to trim quality and adapters if the fastQC report indicates that there is adapter contamination and/or poor quality at the end of sequence reads.


	# Trim adapters and quality
	trim_galore -o ./ --clip_R1 5 --three_prime_clip_R1 5  SRR534177.fastq.gz

## Clipping WGBS data