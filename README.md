
<p align="center">
<img align="center"   src="/img/TAPAS.jpeg" alt="TAPAS Logo">
</p>


Transcription Analysis Plus Alternative Splicing for RNA-SEQ
=============


[TAPAS](https://github.com/LucoLab/RNASEQ) is the workflow we use in the lab to process RNA-SEQ.


---

<p align="center">1. Quick overview</p>
<p align="center">2. Create Listing File</p>
<p align="center">3. Create Init File</p>
<p align="center">4. Set up Tools</p>
<p align="center">5. Set up Annotation Files</p>

## Quick overview

One command is used to launch all the steps :

- Alignment with STAR to generate bigWigs for visualisation
- Differential Expression with DESEQ2
- Splicing Analysis with WHIPPET 10.4
- Custom filters, formats ...

```shell
    python3 signature.py -l listofRnaSeqProject.tsv -i init_disciplus.json
```

You end with bed files/excel describing spliced exons between two conditions.

Two files need to be created before.


## Create listing File 

![alt text](https://github.com/LucoLab/RNASEQ/blob/master/img/listing.png "Listing")

- **STUDY :** Can be anything you want.
- **RUNACCESSION :** Can be anything you want.
- **LIBRARY_LAYOUT :** KeyWords to use are "PAIRED" or "SINGLE"
- **FASTQ :** Absolute path to your dataset or ftp url to download automatically, if paired-end it can be data_1_fastq.gz;data_2_fastq.gz  or data_R1_fastq.gz;data_R2_fastq.gz 
- **TREATMENT :** Anything you want
- **CELL_LINE :** Anything you want
- **CONDITION :** Keywords to use are "TEST" and "CONTROL"
- **REP_NUMBER :** should be 1,2,3 or rep1, rep2, rep3
- **TREATMENT_DAY :** A number.
- **KMER :** Keywords to use are "normal" or "short". Use "normal" everytime. Short is used when read are less than 50 and you want to give a try with whippet with lowers index value for kmersize.

Column names has to be the sames as showed in the picture. The order of the columns can be changed.


__Very Important__  

1. You create a directory where you put your fastq files. The results will be in this directory. "Output" path in init.json must be the same as in listing.tsv in column fastq. Input and Output directories sould be the same. 

2. You also have to create a directory called FINAL in this directory.


## Create Init File 


You need to set up a _init.json_ . You will set up path to softwares and other annotation files.

First of all, you need to change init.json inside config directory to set up correctly paths and files.

---


## Under the Hood / Go Fast for Splicing


Whippet works in two commands. First you create the .psi files with _whippet-quant.jl_

Then you do a differential analysis using _whippet-delta.jl_ to create a .diff file.

Finally, you apply fancy filters to keep reliable events. (P>=0.95, only CE , |dpsi| >= 0.1)

A Rscript (annotSymbol.R) is used to add the gene symbol.

```shell

/home/jean-philippe.villemin/bin/julia-6/julia /home/jean-philippe.villemin/bin/Whippet.0.10.4/bin/whippet-quant.jl -x /home/jean-philippe.v
illemin/data/data/index_whippet/julia-6/human/9kmer/index_whippet  path2/A_1.fastq.gz path2/A_2.fastq.gz -o path2/A

/home/jean-philippe.villemin/bin/julia-6/julia /home/jean-philippe.villemin/bin/Whippet.0.10.4/bin/whippet-delta.jl -a A.psi.gz,B.psi.gz,C.psi.gz -b D.psi.gz,E.psi.gz,F.psi.gz

PATH2=/somewhereYouMustConfig/
FILE=output
TYPE="CE"

/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="CE" && sqrt($8*$8) >= sqrt(0.1*0.1) && $9 >= 0.95) print  ary[1],$3,$4,$5,$6,$7,$8,$9,$10,$11 ; }' ${PATH2}${FILE}.diff > ${PATH2}${FILE}.clean.${TYPE}.diff
/bin/sed -i $'1 i\\\ngene\tcoordinates\tstrand\tevent\tpsiA\tpsiB\tdpsi\tprobability\tcomplexity\tentropy'  ${PATH2}${FILE}.clean.${TYPE}.diff
Rscript /home/luco/work/annotSymbol.R --organism=human --file=${PATH2}${FILE}.clean.CE.diff
/usr/bin/gawk -F ","  'BEGIN {OFS="_";}  {  print $2,$1,$4; }'  ${PATH2}${FILE}.clean.${TYPE}.diff.annoted.csv > ${PATH2}${FILE}.final.txt
```

If you have discrepency with this go fast method and the pipeline it's because pipeline, remove some events that are on weakly expressed gene.

---

## Notes For Expression  

Genes differentially expressed are with a | foldchange | > 1.5 and p-value-adjusted < 0.05

You will find them in directory final.

DESEQ_res_annotated_DE.csv (genes up and down significant)  

DESEQ_all_res_annoted_sorted_pvalAdj.csv contains everything without filtering.

In the outpout of STAR mapper, you will also have a file for TPM computed by TPMcalulator. gene_TPM.txt.  

I used the column TPMexon.



---

The next part explain what you need to do if you start from scratch (softs to install, files to download or create)

You will need to download some files from Gencode, create a genome Index with STAR (explain in Alignment part), and create some custom files for which scripts are provided to generate them before using the main pipeline.

An example of the init.json is showed in config directory.  
In the next section, I explain all the files you need.

---
---
---

## Set up Annotation Files

Much more of the annotation files are given in the annot directory.

Download  from ensembl biomart a text file with _Ensembl Gene ID_ &	_Associated Gene Name_
Should look as follows (you can also the one provided) : 

```
Ensembl Gene ID	Associated Gene Name
ENSG00000252760	RNA5SP54
ENSG00000252830	5S_rRNA
```

Now you can set up _genes_biomart_ensembl_ value in init.json  

Download from gencode the gtf file to set up _path_to_gtf_.  

Download  fasta sequence of the genome of interest. (Here we use PRIM_ASSEMBLY from GENCODE R25 based on ENSEMBL 85 for Human and GENCODE M15 based on ENSEMBL 90 for Mouse.  

Create this file with chromosome files.  

```shell
conda install pyfaidx
faidx yourgenome.fasta -i chromsizes > yourgenome.genome.sizes
```

Then configure in init.json the value for _path_to_chrom_length_ref_ & _genomeDir_ with the directory where you create yourgenome.genome.sizes directory as showed in the provided example

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

This gives you the value for _transcriptome_index_ in init.json

Also run this R script. It will create two files : one with genomic non-redundant exons coordinates and the other with gene length using sum of exon length. In annot, you have two zipped examples for the human Gencode Annotation version 25 (Ensembl 85) for these two files. These files are also provided in annotation dir.

- exons.non-redundant.csv
- gene.length.exon.sums.csv


```R  
Rscript geneSize.R -f Yourannotation.gtf
```
 
Now you can complete _postitionGenomicExon_ and _gene_length_ values in init.json.


NB : Andrew pointed out you need to covert exons.non-redundant.csv to a bed format.

--
---
---


## Set up Tools


Check that the following tools are installed on server/computer.


Scripts available here are in Python3.  
It's  advised to install Conda if python3 is not set up on your computer.   
It will make things easier then for installing tools or switching to older python version if needed.

To be sure python3 is there, you can call it as follows : 

```shell
python3 --version
```

Should return something like  : _Python 3.5.5 :: Anaconda custom (64-bit)_

_**Conda**_ : [here](https://www.continuum.io/downloads)



Then you can install all tierce tools needed by the pipeline. See below.


Follow the links below and look the _Set up tools_ paragraph for each section.

-ALIGNEMENT : [here](https://github.com/LucoLab/RNASEQ/blob/master/ALIGNEMENT.md)  

-DIFF_EXPRESSION : [here](https://github.com/LucoLab/RNASEQ/blob/master/DIFF_EXP.md)  

__NB__ : RMysql package need to be installed in R also.

-SPLICING :

_**JULIA 0.6**_ : [here](https://julialang.org/downloads/)  

 ```shell
 	wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.3-linux-x86_64.tar.gz
 	tar -zxvf julia-0.6.3-linux-x86_64.tar.gz
 ```

_**WHIPPET 10.4**_ : [here](https://github.com/timbitz/Whippet.jl)  

In julia :  

```shell
	Pkg.update()
	Pkg.add("Whippet")
	using Whippet
```


_**bedtools**_ : Bed handler [here](http://bedtools.readthedocs.io/en/latest/index.html)

