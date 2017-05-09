# RNA-SEQ

---

Here we describe the different steps involved in the analysis of RNASEQ in the laboratory.  
This script has been first developed for stranded pair-end RNASEQ.  

This tutorial goal is to make people understand how the raw data has been process in the team (tools and workflow) and to help them to reproduce analysis.


---

1. Set up Environment
2. Set up ConfigFile for alignment only
3. Launch Alignment
4. Compute Statistics
5. Gene Expression Analysis
5. Gene Ontology Analysis
7. Alternative Splicing


## Set up Environment

---

First you need to check that the following tools are installed on server/computer.

_**STAR**_ : Aligner [here](https://github.com/alexdobin/STAR)

_**Samtools**_ : Bam handler [here](http://www.htslib.org/download/)

_**wigToBigWig**_ : Include in KentTools suite [here](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)

_**rMATS**_ : For Alternative Splicing Analysis [here](http://rnaseq-mats.sourceforge.net/index.html)

Scripts available here are in Python3.  
It's not required but advised to install Conda if python3 is not set up on your computer.   
It will make things easier then for installing tools or switching to older python version if needed.

_**Conda**_ : [here](https://www.continuum.io/downloads)

_**Salmon**_ : Compute TPM values from fastq [here](https://github.com/COMBINE-lab/salmon)

_**LiftOver**_ : Include in KentTools suite. Lift Over is used to transform coordinates between several genome build [here](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)

_**WiggleTools**_ : _Optionnal_ , can be useful when you want to merge bigwig files and apply some basic numerical operations [here](https://github.com/Ensembl/WiggleTools)

## Set up ConfigFile.json

---

We use a json file to create a configuration file before doing your alignment. 

The json file is a *key:value* listing which defines all parameters for the pipeline.
- number of cores used
- path to output/input directory
- name of analyse
- ...


## Launch Alignment

---

```shell
	python3 pathTo/align_V1.py -c pathToConfigFile/condition.json
```

It will launch an alignment with STAR.
It will creates a directory (path defined in your configuration file) with inside :

1. Bam file.
2. Bigwig files per strand normalized by CPM (count per millions) for visualization purpose in UCSC.
3. Reads count per gene - Sample_ReadsPerGene.tab (used for the gene expression step).
4. Intermediate files used in the following steps of the pipeline.


You need to use this on each sample/replicate.

_**Manually create your reads count per gene matrice :**_



You have to manually grap each *read count per Gene* column in each _Sample_ReadsPerGene.tab_ to make only one file.

Sort each file by geneName column before joining.  
Depending on your library type , you should pick first or second column. 
Select the output according to the strandedness of your data.

Here kit used for library preparation is of type directional first strand that means R2 is always forward.  
The second read (read 2) is from the original RNA strand/template, first read (read 1) is from the opposite strand. 

Tcheck [here](https://www.youtube.com/watch?v=n2XEsw7EJLw&feature=youtu.be)  to get more details on the library preparation.
[here](http://chipster.csc.fi/manual/library-type-summary.html)  
Be careful, different vocabulary is used depending of the softwares : ISR(Salmon),fr-secondstrand(RMATS)

Select second column which name is *GoodPairs(F2R1+/F1R2-)*  
At the end, your file should looks like what follows :

| GeneName  | Condition1Replicat1    | Condition1Replicat2  | Condition2Replicat1  | ... |
| ---   |:---:| ---: | ----:| -----:|
| Gene1 | 100 | 140  | ...  |
| Gene2 | 0   |  0   | ...  |
| Gene3 | ...  | ... | ...  |

## Compute Statistics

---

This will create a file named *statistics* in you current directory with :

1. NameOfYourSample
2. Input Reads
3. Uniquely Mapped (%)
4. Multimapping (%)

```shell
	python3 pathTo/computeStatistics.py -d pathToDirOfSample -n nameOfYourSample
```
You need to use this on each sample/replicate.


## Gene Expression Analysis

First, you need to create a csv file describing your experiments.
Column name can be changed but then you need to make changes directly inside retrieved code.  
*Sample_id* and *condition* columns are mandatory.  
It should looks as follows :

| sample_id  | condition    | group  |
| ---   |:---:| ---: | 
| Cond1rep1 | Cond1 | 2016  | 
| Cond1rep2 | Cond1   |  2016   | 
| Cond1rep1 | Cond1  |  2015  |

This will be your first input.
Then you will pass also the matrice containing all your reads count per gene and per sample.You did that the step before.  
You have to give the names of the condition you want to compare between each others to make the script works correctly.

NB  : diff_expRscript need some R packages to work. After first launch if packages are not installed it will complain and tell you the package missing.

```R
	Rscript diff_expRscript ${PATH_TO_SCRIPT}/diff_exp.R  --dir ${PATH_TO_DATA}/[DIR_NAME] --cond1 [COND1]  --cond2 [COND2]  ${PATH_TO_DATA}/[DESIGN.csv] ${PATH_TO_DATA}/[GENE_READ_COUNT.csv] > ${PATH_TO_DATA}/logR/[ANY_NAME].out
```
This script will give you list of differentially expressed genes with fold changes and associated pvalues. Graphics ( PCA plot -see above-, Heatmap, Histogram , ...) that describes quality of your data will also be created.
![PCA](https://github.com/ZheFrenchKitchen/pics/blob/master/pca2015.png)



One file with all pvalue and fc will be created per comparison.

You can merge all this files into one. You need first to sort them by gene ensembl id column to be sure the fusion will be correct.

> _Note1 :_  Look in logR directory, the file finishing by 'out' contains number of genes up & down.

> _Note2 :_  There is also others text files created which already apply filter on the foldchange(log2FC > 1.5) .You can change that directly in the code if you want.

> _**Note3 :**_ Important, you need to modify inside the code the path to the files because they are hard-coded.  

```R
	Rscript merged_all_FC_files ${PATH_TO_SCRIPT}/merged_all_FC_files.R  --dir ${PATH_TO_DATA}/[DIR_OUTPUT]  
```
At the end, it will create a file (*FC_collapse_with_pval.csv*) in the directory you passed as parameter.

## Gene Ontology Analysis

Plot are automatically created in the ${PATH_TO_DATA}/[DIR_NAME]/pathway . It gives you enriched ontology terms for the genes which are differentially expressed using a cutoff of | FC | > 1.5.  

![REACTOME](https://github.com/ZheFrenchKitchen/pics/blob/master/DESEQ_REACTOME_DOTPLOT.png)

You also can separately use enrichment.R script to pass a list of EnsemblGeneID to plot the same graphics.  

```R
	Rscript enrichment.R --file ${PATH_TO_MY_GENE_LIST}  
```
## Computing TPM with Salmon

TPM - Transcripts Per Million - computed here let you compare gene expression inside a condition.
The following script will launch Salmon tools to compute this values for all samples using the json configuration file you set up at the beginning of the tutorial.  
It assuming that you didn't change files location of course.

```shell
	python3 pathTo/salmonGoToTheRiver.py -c pathToConfigFile/condition.json
```
All R1 and R2 fastq files are concatenated here separately and are then used by Salmon to create a file quant.genes.sf per condition as follows : 

| Name |	Length |	EffectiveLength |	TPM |	NumReads |
| ---   | --- | --- | --- | --- | 
| ENSG00000210678.1	| 68 |	9.78894 |	13564.3 |	9021 |
| ENSG00000210196.2	| 66 |	9.45298 |	308.302 |	198 |
| ENSG00000210194.1	| 71 |	10.2649 |	3903.13 |	2722 |

These files will be used to filter out some splicing events by gene TPM value when expression is very low.

## Disjoin Bam by Read Strand for stranded RNASEQ  (optional)

For several purpose, you could need to separate bam file per Read Strand (specific to stranded RNASEQ).  
This can be useful when looking at visualization in UCSC genome browser for designing primers ...

```shell
rnaseq_dispatch_stranded.sh pathToBam NameOfYourChoiceForOutput PathToOutputDir
```

You will end with two bam files for each strand.

## Trim (optional)

Rmats needs to have reads of the same length.
It give a script to trim your fastq.[here](http://rnaseq-mats.sourceforge.net/user_guide.htm)
To trim the poor quality 3' end of reads, use the trimFastq.py script found in the bin directory of RMATS.  
This script is in python2. (see also in Alternative Splicing section to switch to python2 env)

```shell
	source activate python2
	python bin/trimFastq.py pathToFastQ pathToTrimmedOutputFastQtrimmed.fastq DesiredLength
```

## Alternative Splicing Analysis

* _**Extract all fastq.gz.**_  

Remember that Salmon concatenated and create big R1/R2 fastq gunzipped files.  
You need to unzip these files to work with RMATS.

```shell
	gunzip pathToFastq.gz
```
It will uncrompress fastq.gz file in the same directory.

* _**Run Rmats**_ [here](http://rnaseq-mats.sourceforge.net/index.html)

Rmats works only with python2.
If your system is in python3 , you still can easily create with conda a python2 environment.

Do this once when you are connected to your server : 
```shell
	conda create -n python2 python=2.7 anaconda
```

Then switch to you python2 environment before launching Rmats scripts: 
```shell
	source activate python2
	python RNASeq-MATS.py -s1 rep1_1[:rep1_2][,rep2_1[:rep2_2]]* -s2 rep1_1[:rep1_2][,rep2_1[:rep2_2]]* -gtf gtfFile -bi STARindexFolder -o outDir -t paired -len readLength -c 0.1 -libType r-secondstrand -analysis P
```

* _**Get  final excel file with alternative events :**_ 

You need to create another json configuration file. 
The file should descrive the following keys : 

1. _**path_to_output**_: clearly enough  
2. _**list_files_splicing**_ : Values can be A3SS,A5SS,SE,RI,MXE . You can use all at the same time but we advice to separate each event in separate configuration json file.

3. _**analysis**_ :  Name of comparison between two conditions   Cond1_vs_Cond2   Always Test_vs_Control   
3.1.1. _**pathToRMATSOutput**_: splicing_dirpath path to MATS_output for comparison   
3.1.2. _**expression_file_path**_ : path DESEQ_all_res_annotated_sorted_pvalAdj_cond1_vs_cond2.csv   
3.1.3. _**sample_control_for_quantification**_ : Name of the condition Cond2   
3.1.4. _**sample_control_for_quantification**_ : Name of the condition Cond1

4. _**samples_for_quantification**_ :  
4.1. _**Name of each condition**_  ex: Cond1
4.1.1 _**replicat_id & path to quant_genes.sf**_ file

5. _**tabs**_ :  Describe the different tabs you want in the final excel.

First you want the analysis you have done between all possible tuples. Here you will have 3 tabs in your file.
> "Cond1 vs Cond2": { "analysisTocheck": ["Cond1_vs_Cond2"],"0":"nevermind"	}, #Comparison1   
> "Cond1 vs Cond2": { "analysisTocheck": ["Cond1_vs_Cond2"],"0":"nevermind"	}, #Comparison2   
> "Cond2 vs Cond3": { "analysisTocheck": ["Cond2_vs_Cond3"],"0":"nevermind"	},

Then you want also , tabs for intersection of events that are present in Comparison 1 and Comparison2.  
You want to tcheck if psi is increasing, decreasing or anyway. In this case you will have 3 tabs in your file.

> "Cond1 vs Cond2 Intersect Cond2 vs Cond3 Decrease": { "analysisTocheck": ["Cond1_vs_Cond2","Cond2_vs_Cond3"],"1":"upper"	},  
> "Cond1 vs Cond2 Intersect Cond2 vs Cond3 Increase": { "analysisTocheck": ["Cond1_vs_Cond2","Cond2_vs_Cond3"],"1":"lower"	},  
> "Cond1 vs Cond2 Intersect Cond2 vs Cond3": { "analysisTocheck": ["Cond1_vs_Cond2","Cond2_vs_Cond3"],"1":"idem"	}

At the end , you will have 6 tabs in your excel but it could be more or less depending what you define in the json config file. Look at the json example in configs.

```shell
  python3 ./mergeFinal.py -c /pathTo/AnyNameSplicing.json 
```

You should end up with an excel file looking like this : 

![Excel](https://github.com/ZheFrenchKitchen/pics/blob/master/fileExcel.png)



