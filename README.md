
# SAGC 5mC DNA Methylation tutorial

September 21st 2021  
Presenters: Jimmy Breen, Bastien Llamas, Melanie Smith, Nathan Watson-Haigh  

## What you'll need  

- Laptop that can access the command-line
    - Preferrably a Mac or Linux but newer windows versions have the windows linux subsystem capabilities (https://docs.microsoft.com/en-us/windows/wsl/install-win10)

For this tutorial we'll need the following programs:
- `bowtie2`
- [bismark](http://www.bioinformatics.babraham.ac.uk/projects/download.html#bismark)
- `samtools`
- [MethylDacyl](https://github.com/dpryan79/MethylDackel)
- `trim-galore`
- `bedtools`

## Workshop schedule

| Time      | Description                                           |
|-----------|-------------------------------------------------------|
| 900       | Registration                                          |
| 910-920   | Introduction to SAGC                                  |
| 920-945   | DNA methylation intro (Bastien)                       |
| 945-1010  | miRNA-seq intro (Melanie)                             |
| 1010-1035 | Histone Modifications (Jim)                           |
| 1035-1050 |                                       [Coffee Break]  |
| 1050-1115 | Chromatin Accessibility (Natalie)                     |
| 1115-1140 | Chromosome conformation capture (Jimmy)               |
| 1140-1200 | Hands-on bioinformatics tutorial setup (Jimmy/Nathan) |
| 1200-1230 |                                               [Lunch] |
| 1230-1400 | Hands-on bioinformatics tutorial (cont.)              |
| 1400-1430 |                                        [Coffee Break] |
| 1430-1600 | Hands-on bioinformatics tutorial (cont.)              |
| 1600-1700 |                              [Drinks at the West Oak] |

## Presenters

<img align="left" src="https://researchers.adelaide.edu.au/sites/default/files/styles/profile_large/public/profile-images/10840.jpeg" width="180" height="180" />
<u>1. Bastien Llamas</u>
  
Bastien is an Associate Professor of Ancient DNA at the University of Adelaide and investigates a range of genetic and epigenetic mechanisms and host-microbiome interactions that facilitate human adaptation to diverse environmental and cultural stressors. He uses a range of advanced new analytical methods (e.g., long-read sequencing, genome graphs) to integrate past and present Indigenous genetic diversity from populations around the world into a new human pangenome reference.
  
  
<img align="right" src="https://www.adelaide.edu.au/directory/melanie.smith?attr=data;dsn=directory.image;field=image;id=56547;m=view" width="200" height="200" />
<u>2. Melanie Smith</u>
  
Melanie is a Postdoctoral Researcher in Professor Claire Roberts's lab at Flinders University. Her research is focused on investigating the role of miRNAs to regulate gene expression in the placenta during pregnancy. She is focused on determining molecular biomarkers for monitoring placenta and pregnancy help using multiomics and bioinformatics analysis methods.
  

<img align="left" src="https://www.sahmriresearch.org/user_assets/7af8aac6c4da0cfb8ccfb1ba486c6d74b5992988/dimitrios_cakouros_cropped.jpg" width="180" height="200" />
<u>3. Dimitrios (Jim) Cakouros</u>
  
Jim is a Postdoctoral Researcher in Professor Stan Gronthos's lab at SAHMRI where his research is focused on investigating epigenetic mechanisms employed by mesenchymal stem cells (MSC) to regulate stem cell renewal and lineage determination. He uses genome wide chromatin immunoprecipitation (ChIP), drug inhibition and conditional knockout models to investigate the function of epigenetic enzymes in skeletal development and bone related diseases such as osteoporosis, fractures and saethre chotsen syndrome.
  
  

<img align="right" src="https://portal.sahmriresearch.org/files-asset/35193617/Stevens.Natalie_Dr._Precision_Medicine_3_PURE.jpg" width="160" height="200" />
<u>4. Natalie Stevens</u>
  
Natalie is a Postdoctoral Researcher in Immunology in Professor David Lynn's EMBL Computational and Systems Biology Laboratory. She works to understand the consequences of immune challenges on longterm innate immune function and pathogen defence. Natalie is passionate about evidence-based medicine, communicating science to the public and advocating for science literacy.
  

<img align="left" src="https://pbs.twimg.com/profile_images/1232553882804350976/R7_bUSmc_400x400.jpg" width="180" height="200" />
<u>5. Jimmy Breen</u>
  
Jimmy is Bioinformatics leader at the SAGC, specialising in many bioinformatics analysis techniques mainly centred around Functional Genomics. He leads a research group at the Robinson Research Institute (University of Adelaide) developing methods to analyse multiomics datasets, primarily focused on gene regulation in human reproductive systems.


## Tutorial

1. [Whole Genome Bisulfite Sequencing (WGBS) introduction](tutorial/01_WGBS_intro.md)
2. [Data prep & quality control](tutorial/02_BS_quality_control.md)
3. [Reference genome mapping](tutorial/03_BS_mapping_example.md)
4. [Post mapping quality control & methylation calling](tutorial/04_methylation_calling.md)
5. [Data integration](tutorial/05_data_integration.md)

### Not discussed but information available

- [Identification of Differentially Methylated Regions (DMRs)](tutorial/DMR_analysisR.md)

### Full dataset download

- [Download the full dataset for future work](tutorial/fillDataDownload.md)