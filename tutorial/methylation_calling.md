## Post mapping quality control & methylation calling

After we're done with mapping our data, we need to process our BAMs. We ned to do the following tasks:

1. Sort our bam file
2. Deduplicate the data to remove clonal deduplicates
3. Make methylation call files

Firstly to sort and deduplicate. One of the most widely-used Bioinformatics tools used today is samtools. And it can be frustrating to use because it lacks multi-threading options which would allow scalability to larger systems.

However now we can use sambamba which fills those needs. Lets use it to sort our bismark-made BAMs and mark and remove duplicates.

	sambamba sort -t 2 met1.TAIR10.bam
	sambamba sort -t 2 colWT.TAIR10.bam

	sambamba markdup -p -r -t 2 met1.TAIR10.sorted.bam met1.TAIR10.sorted.markdup.bam
	sambamba markdup -p -r -t 2 colWT.TAIR10.sorted.bam colWT.TAIR10.sorted.markdup.bam

Now we have two cleaned bam files that 

## Methylation calls



Now we have an alignment for each of our samples, we need to extract the respective 5mC contexts from the data. Generally most mammalian BSseq sequencing project can get away with only identifying the proportion of methylated cytosines in CpG conetxts, but with whole-genome data it is always recommended to look at the other contexts to get a wide picture of the methylation profile.

Bismark does come with methylation call extraction scripts, however I have found that the program PileOMeth is very quick and easy to use.

	for i in *.sorted.markdup.bam
 	 do
  		PileOMeth extract --mergeContext -o ${i%%.*}.CpG Athal/TAIR10.fasta $i
  		PileOMeth extract --noCpG -CHG --mergeContext -o ${i%%.*}.CHG Athal/TAIR10.fasta $i
  		PileOMeth extract --noCpG -CHH --mergeContext -o ${i%%.*}.CHH Athal/TAIR10.fasta $i
	done

