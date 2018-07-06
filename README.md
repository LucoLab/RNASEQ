
<p align="center">
<img align="center"   src="/img/TAPAS.jpeg" alt="TAPAS Logo">
</p>


Transcription Analysis Plus Alternative Splicing for RNA-SEQ
=============


[TAPAS](https://github.com/LucoLab/RNASEQ) is the worklow we use in the lab to process RNA-SEQ.

Here we describe the different steps involved in this workflow.

The goal of this tutorial is to help people to reproduce the analysis.

---

1. Quick overview
2. Set up Tools
4. Set up Files
2. Listing of projects to analyse



## Quick overview


Once you installed all the necessary softs and files, you need to set up a _init.json_ file inside a directory called _config_ in the RNA-SEQ directory. 

And it's done.

One command is used to launch all steps :

```shell
    python3 signature.py -l listofRnaSeqProject.tsv
```

- Alignment with STAR to generate bigWigs for visualisation
- Transcript Quantification with SALMON
- Differential Expression with DESEQ2
- Splicing Analysis with WHIPPET
- Splicing Analysis with RMATS
- Custom filters, formats ...

You end with bed files describing exons more splice In or splice Out between two conditions.

Each step can be launched separately if you need to (nevertheless previous steps need to be completed).

_**ALIGNEMENT**_ : [here](https://github.com/LucoLab/RNASEQ/blob/master/ALIGNEMENT.md)
_**DIFF_EXPRESSION**_ : [here](https://github.com/LucoLab/RNASEQ/blob/master/DIFF_EXP.md)


Here we wrote a resume describing the sequence of scripts and what they do.

_**WHOLE STEPS**_ : [here](https://github.com/LucoLab/RNASEQ/blob/master/ALLSCRIPTS.md)


## Init file


First of all, you need to change init.json inside config directory to set up correctly paths and files.

You will need to download some files from Gencode, create a genome Index with STAR, and create some custom files for which scripts are provided to generate them before using the main pipeline.


## Set up Tools


First you need to check that the following tools are installed on server/computer.


Scripts available here are in Python3.  
It's not required but advised to install Conda if python3 is not set up on your computer.   
It will make things easier then for installing tools or switching to older python version if needed.

To be sure python3 is there, you can call it as follows : 

```shell
python3 --version
```

Should return something like  : _Python 3.5.5 :: Anaconda custom (64-bit)_

Then you can install all tierce tools needed by the pipeline. See below.


_**Conda**_ : [here](https://www.continuum.io/downloads)

Follow the link and the look set up tools paragraph for each section.

_**ALIGNEMENT**_ : [here](https://github.com/LucoLab/RNASEQ/blob/master/ALIGNEMENT.md)
_**DIFF_EXPRESSION**_ : [here](https://github.com/LucoLab/RNASEQ/blob/master/DIFF_EXP.md)


## Set up Files


Download  fasta sequence of the genome of interest. (Here we use PRIM_ASSEMBLY from GENCODE R25 based on ENSEMBL 85 for Human and GENCODE M15 based on ENSEMBL 90 for Mouse.  
Create this file with chromosome files.  

```shell
conda install pyfaidx
faidx yourgenome.fasta -i chromsizes > yourgenome.genome.sizes
```

Then configure in init.json the value for _path_to_chrom_length_ref_ with the directory where you create yourgenome.genome.sizes directory

Also you need to create index for STAR :  

```
STAR --runMode genomeGenerate --genomeDir /home/jean-philippe.villemin/mount_archive2/commun.luco/ref/indexes/GRCm38_PRIM_GENCODE_M15/ --genomeFastaFiles /home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genome/GRCm38_PRIM_GENCODE_M15/GRCm38.primary_assembly.genome.fa --runThreadN 30 --sjdbGTFfile /home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCm38_PRIM_GENCODE_M15/gencode.vM15.primary_assembly.annotation.gtf
```

Generate the fasta file index :  
```shell
samtools faidx reference.fa 
```
This creates a file called reference.fa.fai

Generate the sequence genome.2bit :  
```shell
faToTwoBit genome.fa genome.2bit
```


Download the transcripts from gencode for your genome version and then inside a dir called transcriptome for instance,
you will create indexes for  Salmon as follows : 

```shell  
salmon index -t gencode.vM15.transcripts.fa.gz -i gencode.m15.transcripts.index --type quasi -k 31  
```

Also run this R script. It will create two files : one with genomic non-redundant exons coordinates and the other with gene length using sum of exon length. 

- exons.non-redundant.csv
- gene.length.exon.sums.csv


```R  
Rscript geneSize.R -f Yourannotation.gtf
```


## Listing of projects to analyse


The file listing all the RNA-SEQ Projects you want to analysis should be as follows :

![alt text](https://github.com/LucoLab/RNASEQ/blob/master/img/listing.png "Listing")

Column names has to be the sames as showed in the picture. The order of the columns can be changed.


- **STUDY :** Can be anything you want.
- **RUNACCESSION :** Can be anything you want.
- **LIBRARY_LAYOUT :** KeyWords to use are "PAIRED" or "SINGLE"
- **FASTQ :** Absolute path to your dataset or ftp url to download automatically, if paired-end it can be data_1_fastq.gz;data_2_fastq.gz  or data_R1_fastq.gz;data_R2_- fastq.gz 
- **TREATMENT :** Anything you want
- **CELL_LINE :** Anything you want
- **CONDITION :** Keywords to use are "TEST" and "CONTROL"
- **REP_NUMBER :** should be 1,2,3 or rep1, rep2, rep3
- **TREATMENT_DAY :** A number.
- **KMER :** Keywords to use are "normal" or "short". Short is used when read are less than 50 and you want to give a try with whippet with lowers index value for kmersize.




