

## WGBS Mapping

Because we're dealing with DNA that has been treated with Sodium bisulfite, all cytosines (C) will be converted to a Uracil (U), while 5' methyl cytosine's will be left unconverted as a cytosine (C). During amplification the converted Uracil's (U) are read as Thymine's (T)

![bisulfite conversion](https://www.diagenode.com/img/applications/bisulfite.png)  
*www.diagenode.com*

There are multiple programs that can do this sort of analysis.
Each implements a strategy of using existing global read aligners, such as `bwa` and `bowtie(1/2)`.

## Bismark

## Preparing the Reference Genome

The reference genome of _Arabidopsis_thaliana_ is in the /data directory, so lets format that for bisulfite sequence mapping. We need to make sure that bismark recognises that you want to use bowtie2 to do the mapping.


	# Create directory for bismark BS-seq genome
	mkdir Athal

	# Move data file in there
	mv data/TAIR10.fa Athal/

	# Run bismark to format our Athaliana genome (TAIR10)
	bismark_genome_preparation --bowtie2 Athal

Now our genome 

## Mapping

Single end data

```
bismark --bowtie2 -p 8 --bam --temp_dir $TEMP_DIR --path_to_bowtie $BOWTIE_PATH/ -o $OUTPUT_DIR $GENOME_PATH -U 
bismark -p 2 --bam --bowtie2 Athal colWT_trim.truncated.gz
```

## Data output

## Summary statistics

- How much has mapped


