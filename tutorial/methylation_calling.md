# Post mapping quality control & methylation calling

After we're done with mapping our data, we need to process our BAMs. We ned to do the following tasks:

1. Deduplicate the data to remove clonal deduplicates
2. Sort based on reference genome coordinate
3. Make methylation call files

## Deduplication


## Resort for methylation calling

## Methylation calls

### Best practises

- Methylation bias
- Coverage considerations
- Sequencing quality

Now we have an alignment for each of our samples, we need to extract the respective 5mC contexts from the data. Generally most mammalian BSseq sequencing project can get away with only identifying the proportion of methylated cytosines in CpG conetxts, but with whole-genome data it is always recommended to look at the other contexts to get a wide picture of the methylation profile.

Bismark does come with methylation call extraction scripts, however I have found that the program PileOMeth is very quick and easy to use.

	for i in *.sorted.markdup.bam
 	 do
  		PileOMeth extract --mergeContext -o ${i%%.*}.CpG Athal/TAIR10.fasta $i
  		PileOMeth extract --noCpG -CHG --mergeContext -o ${i%%.*}.CHG Athal/TAIR10.fasta $i
  		PileOMeth extract --noCpG -CHH --mergeContext -o ${i%%.*}.CHH Athal/TAIR10.fasta $i
	done

## Sanity check

- WT vs MET1
- Chr1 vs Chloroplast

