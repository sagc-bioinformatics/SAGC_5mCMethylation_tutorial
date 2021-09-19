# Post mapping quality control & methylation calling

After we're done with mapping our data, we need to process our BAMs. 
We ned to do the following tasks:

1. Deduplicate the data to remove clonal deduplicates
2. Sort based on reference genome coordinate
3. Make methylation call files

## Deduplication


## Methylation calls

"Methylation calling" is a process of quantifying the proportion of methylation at a specific 5mC context.
At each individual Cytosine site across the reference sequence, converted and un-converted cytosine bases are tallied and quantified to a specific proportion of methylation (i.e. methylated C/(methylated C + non-methylated C)).
Generally most mammalian BSseq sequencing project can get away with only identifying the proportion of methylated cytosines in CpG conetxts, but in other species 

Bismark does come with methylation call extraction scripts, however I have found that the program `MethylDacyl` (formerly `PileOMeth`) is a very quick caller that is multi-threaded and easy to use.

Before we start to call methylated sites

- Methylation bias
- Coverage considerations
- Sequencing quality

Now we have an alignment for each of our samples, we need to extract the respective 5mC contexts from the data. 




	for i in *.sorted.markdup.bam
 	 do
  		PileOMeth extract --mergeContext -o ${i%%.*}.CpG Athal/TAIR10.fasta $i
  		PileOMeth extract --noCpG -CHG --mergeContext -o ${i%%.*}.CHG Athal/TAIR10.fasta $i
  		PileOMeth extract --noCpG -CHH --mergeContext -o ${i%%.*}.CHH Athal/TAIR10.fasta $i
	done

## Sanity check

- WT vs MET1
- Chr1 vs Chloroplast

---

LASTLY: To finish off we're going to skip the DMR step and get stuck into [data integration](05_data_integration.md)