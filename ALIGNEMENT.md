
<p align="center">
<img align="center"   src="/img/TAPAS.jpeg" alt="TAPAS Logo">
</p>


Transcription Analysis Plus Alternative Splicing for RNA-SEQ
=============


[TAPAS](https://github.com/LucoLab/RNASEQ) is the worklow we use in the lab to process RNA-SEQ.

---

Here we describe the different steps involved in this workflow.

The goal of this tutorial is to help people to reproduce the analysis.


## Quick overview

---

First, be sure you can call python3 as follows : 

```shell
python3 --version
```

Should return something like  : _Python 3.5.5 :: Anaconda custom (64-bit)_

Then you can install all tierce tools needed by the pipeline. See _Set up Tools_ paragraph.

Once you installed all the necessary softs, you need to set up a _init.json_ file inside a directory called _config_ in the RNA-SEQ directory. 

And it's done.

One command is used to launch all steps :

```shell
    python3 signature.py -l listofRnaSeqProject.tsv
```

- Alignment with STAR to generate bigWigs for visualisation
- Transcript Quantification with SALMON
- Differential Expression with DESEQ2
- Splicing Analysis with Whippet
- Splicing Analysis with Rmats
- Merge and Filters

You end with bed files describing exons more splice In or splice Out between two conditions.

Each step can be launched separately if you need to.


## Listing of projects to analyse

---

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


---
1. Quick overview
2. Listing of projects to analyse
3. Set up Tools
4. Set up Files
5. Set up ConfigFile for Alignment.py
6. Launch Alignment
7. Outputs

## Set up Tools

---

First you need to check that the following tools are installed on server/computer.


Scripts available here are in Python3.  
It's not required but advised to install Conda if python3 is not set up on your computer.   
It will make things easier then for installing tools or switching to older python version if needed.

_**Conda**_ : [here](https://www.continuum.io/downloads)

ALIGNMENT : 

_**STAR**_ : Aligner [here](https://github.com/alexdobin/STAR)

_**Salmon**_ : Compute TPM values from fastq [here](https://github.com/COMBINE-lab/salmon)

_**Samtools**_ : Bam handler [here](http://www.htslib.org/download/)

_**wigToBigWig**_ : Include in KentTools suite [here](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)

## Set up Files

Download  fasta sequence of the genome of interest. (Here we use PRIM_ASSEMBLY from GENCODE R25 based on ENSEMBL 85 for Human and GENCODE M15 based on ENSEMBL 90 for Mouse.  
Create this file with chromosome files.  

```shell
conda install pyfaidx
faidx yourgenome.fasta -i chromsizes > yourgenome.genome.sizes
```

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


Last step  : 

Download the transcripts from gencode for your genome version and then inside a dir called transcriptome for instance,
you will create indexes for  Salmon as follows : 

```shell  
salmon index -t gencode.vM15.transcripts.fa.gz -i gencode.m15.transcripts.index --type quasi -k 31  
```
## Set up ConfigFile.json

---

We use a json file to create a configuration file before doing your alignment. 

The json file is a *key:value* listing which defines all parameters for the pipeline.
- number of cores used
- path to output/input directory
- name of analyse
- ...


See config directory. You will find an example called paired.set1_align.json for a test dataset.


## Launch Alignment

---

```shell
	python3 pathTo/alignment.py -c pathToConfigFile/condition.json
```

You can test the pipeline on a  short dataset provided in test directory.

It will launch an alignment with STAR.
It will creates a directory (path defined in your configuration file) with inside :

1. Bam file.
2. Bigwig files per strand normalized by CPM (count per millions) for visualization purpose in UCSC.
3. Reads count per gene - Sample_ReadsPerGene.tab (used for the gene expression step).
4. Intermediate files used in the following steps of the pipeline.


You need to use this on each sample/replicate.

Finally you get the following directories as output : 


![alt text](https://github.com/LucoLab/RNASEQ/blob/master/img/output_alignment.png "Outputs")

![alt text](https://github.com/LucoLab/RNASEQ/blob/master/img/output_alignment_open.png "Outputs")




