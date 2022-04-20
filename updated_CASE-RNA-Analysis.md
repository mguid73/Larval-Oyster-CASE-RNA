# CASE RNASeq Analysis : QC, Read Trimming, Read Alignment, and Transcript Assembly

**Author: Megan Guidry**

Original analysis and markdown authored by [Maggie Schedl](https://github.com/meschedl). Maggie's analysis workflow is annotated in the `CASE-RNA-Analysis.md` file also in this directory. 

**This document utilizes the analyses that Maggie compiled in her workflow but with an _improved reference genome_ for the eastern oyster.**

This markdown references Erin Robert's [RNASeq pipeline](https://github.com/erinroberts/apoptosis_data_pipeline/blob/master/Streamlined%20Pipeline%20Tutorial/Apoptosis_Pipeline_Tutorial_with_Reduced_Dataset.md) or Kevin Wong's [project](https://github.com/jpuritz/BIO_594_2018/blob/master/FinalAssignment/KevinWong_FinalAssignment/P.dam_DE_Analysis.md) in [J. Puritz's](https://github.com/jpuritz) Bio594 2018.  If other resources were used they should be linked in this markdown. If anything is not linked properly, not-sourced, or looks missing please contact at meschedl@uri.edu or mguidry@uri.edu.

_All analysis was done on our lab shared server, KITT, made by [J. Puritz](https://github.com/jpuritz)_  

Programs Installed/Needed for this Project:  
- HISAT2 
- StringTie 
- gffcompare 
- fastp, fastQC, multiqc
- samtools

conda create -n CASE-RNA hisat2 stringtie gffcompare fastp fastQC multiqc samtools

File Naming and Information:

location of reference genome on KITT: `/RAID_STORAGE2/Shared_Data/Oyster_Genome/masked/masked.cvir.genome.fasta`

- `masked.cvir.genome.fasta`: Eastern Oyster genome
- `ref_C_virginica-3.0_top_level.gff3`: annotation file for the Eastern Oyster, I also got this from Erin, it has the matching header line convention to work with the above genome
- `prepDE.py`: python script for converting read count information from StringTie into a matrix fo DESeq2, full code [here](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py)
- `CA_J06`: this is an example of the naming convention of the samples, CA refers to the treatment (coastal acidification), and the J number refers to the replicate jar

-------

## Create conda environment & install packages
```
conda create -n CASE-RNA hisat2 stringtie gffcompare fastp fastQC multiqc samtools
conda activate CASE-RNA
```

`conda list` will allow you to look at all of the packages you have installed in your environment.

## Quality Control and Read Trimming
Steps:
1. 

Reads were already de-multiplexed and assigned to each individual sample by [J. Puritz](https://github.com/jpuritz), and linked to a directory called CASE_RNA.

 I made a "working" directory to work in, and then linked in the files to that directory.
```
mkdir Working-CASE-RNA
cd Working-CASE-RNA
ln -s /home/mguidry/CASE_RNA/* .
ls
```
Output:
```
CA_J06.F.fq.gz  CA_J08.F.fq.gz  CA_J11.F.fq.gz  CA_J18.F.fq.gz  CASE_J03.F.fq.gz  CASE_J09.F.fq.gz  CASE_J12.F.fq.gz  CASE_J13.F.fq.gz  CON_J02.F.fq.gz  CON_J05.F.fq.gz  CON_J10.F.fq.gz  SE_J01.F.fq.gz  SE_J04.F.fq.gz  SE_J07.F.fq.gz
CA_J06.R.fq.gz  CA_J08.R.fq.gz  CA_J11.R.fq.gz  CA_J18.R.fq.gz  CASE_J03.R.fq.gz  CASE_J09.R.fq.gz  CASE_J12.R.fq.gz  CASE_J13.R.fq.gz  CON_J02.R.fq.gz  CON_J05.R.fq.gz  CON_J10.R.fq.gz  SE_J01.R.fq.gz  SE_J04.R.fq.gz  SE_J07.R.fq.gz 
```

### 1. Check out read counts
The first thing I did was look at the read counts for each file. I used a code from [this website](http://www.sixthresearcher.com/list-of-helpful-linux-commands-to-process-fastq-files-from-ngs-experiments/) and made it into a for-loop that went through all the files. 

It outputs the filename and the number of reads in that file.

```
for fq in *.fq.gz
> do
>	echo $fq
> zcat $fq | echo $((`wc -l`/4))
> done
```

Output:

Takes roughly ~1-2min per sample
```
CA_J06.F.fq.gz
20445176
CA_J06.R.fq.gz
20445176
CA_J08.F.fq.gz
21746189
CA_J08.R.fq.gz
21746189
CA_J11.F.fq.gz
25550864
CA_J11.R.fq.gz
25550864
CA_J18.F.fq.gz
37263541
CA_J18.R.fq.gz
37263541
CASE_J03.F.fq.gz
26925142
CASE_J03.R.fq.gz
26925142
CASE_J09.F.fq.gz
31720810
CASE_J09.R.fq.gz
31720810
CASE_J12.F.fq.gz
24582739
CASE_J12.R.fq.gz
24582739
CASE_J13.F.fq.gz
36132924
CASE_J13.R.fq.gz
36132924
CON_J02.F.fq.gz
28850301
CON_J02.R.fq.gz
28850301
CON_J05.F.fq.gz
27446573
CON_J05.R.fq.gz
27446573
CON_J10.F.fq.gz
35291136
CON_J10.R.fq.gz
35291136
SE_J01.F.fq.gz
28376966
SE_J01.R.fq.gz
28376966
SE_J04.F.fq.gz
24827716
SE_J04.R.fq.gz
24827716
SE_J07.F.fq.gz
30894132
SE_J07.R.fq.gz
30894132
```

### 2. Raw read quality check
Next, I looked at the quality of the reads in each `.fq.gz` file using [MultiQC](https://multiqc.info/) & [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). I'll go back and run the same on the trimmed data in the next couple of steps too. But first, we want to get an idea of what the data look like first.

```
mkdir fastqc-raw
cd fastqc-raw/
fastqc ../*fq.gz   #started at 4:51pm
mv *fastqc.* fastqc-raw/
cd fastq-raw/
multiqc .
```

