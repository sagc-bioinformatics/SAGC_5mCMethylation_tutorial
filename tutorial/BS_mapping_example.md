
Because we're dealing with DNA that has been treated with Sodium bisulfite, all cytosines (C) will be converted to a Uracil (U), while 5' methyl cytosine's will be left unconverted as a cytosine (C). During amplification the converted Uracil's (U) are read as Thymine's (T)

![bisulfite conversion](https://www.diagenode.com/img/applications/bisulfite.png)  
*www.diagenode.com*

```
bismark --bowtie2 -p 8 --bam --temp_dir $TEMP_DIR --path_to_bowtie $BOWTIE_PATH/ -o $OUTPUT_DIR $GENOME_PATH -U 
bismark -p 2 --bam --bowtie2 Athal colWT_trim.truncated.gz
```

## Data output

## Summary statistics


