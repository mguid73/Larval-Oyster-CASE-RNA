# CASE RNASeq Analysis : QC, Read Trimming, Read Alignment, and Transcript Assembly

**Author: Megan Guidry**

Original analysis and markdown authored by [Maggie Schedl](https://github.com/meschedl). Maggie's analysis workflow is annotated in the `CASE-RNA-Analysis.md` file also in [this directory](https://github.com/mguid73/Larval-Oyster-CASE-RNA). 

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
2. 
3. 
4. 
5. 


### Set up
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
fastqc ../*fq.gz   #took about 1.5-2 hrs for these data
mv *fastqc.* fastqc-raw/
cd fastq-raw/
ls  #the *fastqc* files should all be in the fastq-raw directory now
multiqc .
```

Now I have a `multiqc_report.html` file which I transfered on to my local computer to look at in a web browser. 
_____________________
There are a couple ways I could have 'transfered' that file over to local. I downloaded it from Visual Studio code's Explorer pane. I could have used the following line to file transfer from KITT to my local computer. This code would be run from my local terminal (not on KITT). Alternatively, I could have used a file transfer software like Cyberduck.
```
ls
scp -P zzzz mguidry@KITT.uri.edu:/home/mguidry/Working-CASE-RNA/fastqc-raw/multiqc_report.html ~/Desktop
```
______________________
Once you have `multiqc_report.html` locally, simply run `open multiqc_report.html` in terminal from the directory that the file is in. This will open the report in your browser and from there you can poke around on it looking through the different module reports.

#### Digging into the raw reads
Here, I will focus on the quality scores of the raw reads. 

![raw mean quality scores](images/raw-fastqc_per_base_sequence_quality_plot.png)
***Mean Phred score at each base in the sequence.*** 

**These quality scores look good! Phred score>30 is an indicator of good quality!** 

Follows the expected trend of a slight decrease in quality at the end of the read. There is a slight dip at the beginning of some of the reads - not sure why? Overall though, this is good quality data.

![raw per sequence quality scores](images/raw-fastqc_per_sequence_quality_scores_plot.png)
***Count of Sequences across quality scores.*** 

**This plot also indicates good quality with a large peak in the number of sequences with a Phred score around 40 and most above 30!** 

However there is still a tail of sequences with lower scores that can be trimmed.
_______________
**There are lots of other visualizations in the multiQC report that have informative data, but for now, I will move on to trimming our raw reads.**
___________

### 3. Trimming raw reads with [`fastp`](https://github.com/OpenGene/fastp)

`fastp` flags:

`-i` = input file, forward

`-I` = input file, reverse

`-o` = output file, forward

`-O` = output file, reverse

`-f` 


In my `Working-CASE_RNA` directory, I used `fastp` with the following code. 
```
mkdir trimmed  #directory for trimmed read files
```

This takes roughly ~2-25 mins per sample. I ran this for each set of paired end reads. I know there is a more elegant way to loop through these though. 
```
fastp -i CA_J06.F.fq.gz -I CA_J06.R.fq.gz -o trimmed/CA_J06.F.trim.fq.gz -O trimmed/CA_J06.R.trim.fq.gz  -f 5 -F 5 -q 15 -u 50 -j trimmed/CA_J06.json -h trimmed/CA_J06.html
```

Check how many files are in the trimmed directory. It should be 56 files total (1 F, 1 R, 1 .html, 1 .json *for each sample*). There are 14 samples.
```
cd trimmed
ls | wc -l  #check total number of files (56)
ls *.json | wc -l   #output: 14
ls *.html | wc -l   #output: 14
ls *.R.trim.fq.gz | wc -l  #output: 14
ls *.F.trim.fq.gz | wc -l  #output: 14
```
Sweet! Now all of the reads have been trimmed with `fastp`. We'll check out the quality scores of the trimmed sequences in the next step. 

### 4. Trimmed read quality check
Running QC like I did for the raw reads. Directing fastqc files and reports to a `fastqc-trimmed` directory. 

```
mkdir fastqc-trimmed
fastqc /trimmed/*fq.gz   #took 2.5-3hrs
cd trimmed
mv *fastqc.* ../fastqc-trimmed/
cd ../fastqc-trimmed/
ls  #the *fastqc* files should all be in the fastqc-trimmed directory now
multiqc .
```

#### Checking out the trimmed reads
Looking at some plots from the MultiQC report we just generated.

Again there is a lot to look at in these reports but here I've just included a few sequence quality plots. 

![trimmed mean quality scores](images/trimmed-fastqc_per_base_sequence_quality_plot-2.png)
***Trimmed data. Mean Phred score at each base in the sequence.*** 


![trimmed per sequence quality scores](images/trimmed-fastqc_per_sequence_quality_scores_plot.png)
***Trimmed data. Count of Sequences across quality scores.*** 

These plots look a little more smoothed out than the pre-QC data!
Next I moved on to mapping those new, trimmed fasta files to the Eastern Oyster Genome.

### 5. Alignment to the *Crassostrea virginica* genome


# Notes to self
Next steps:
* QC on trimmed reads - done 5/10/22
* mapping and aligning 
* fill in the flag meaning for fastp
* add TOC
