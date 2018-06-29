
## Quick overview

---

Here we describe the **differential expression** part of the workflow.

An analyse should not last more than 10 minutes.


## Set up Tools

---

First you need to check that the following tools are installed on server/computer.


Scripts available here are in Python3.  
It's not required but advised to install Conda if python3 is not set up on your computer.   
It will make things easier then for installing tools or switching to older python version if needed.

_**Conda**_ : [here](https://www.continuum.io/downloads)


**DIFFERENTIAL EXPRESSION :**

_**R**_ : R [here](https://pbil.univ-lyon1.fr/CRAN/)


Install following packages for R : 

```r
source("https://bioconductor.org/biocLite.R")
	biocLite("data.table")
	biocLite("reshape2")
	biocLite("edgeR")
	biocLite("DESeq2")
	biocLite("limma")
	biocLite("RColorBrewer")
	biocLite("gplots")
	biocLite("heatmap3")
	biocLite("grDevices")
	biocLite("genefilter")
	biocLite("ggplot2")
	biocLite("GenomicFeatures")
	biocLite("AnnotationDbi")
	biocLite("biomaRt")
	biocLite("stringr")
	biocLite("org.Hs.eg.db")
	biocLite("vsn")
	biocLite("plyr")
	biocLite("pheatmap")
	biocLite("PoiClaClu")
	biocLite("gtools")
```

## Create config file in json format

---

You need an _init.json_ and _diff_exp.json_ to launch this script.

_init.json_ is called automatically. Create the file in in _configs_ directory.

Only scriptDir variable need to be set up in your _init.json_ :

		"scriptDir"                      : "/home/jean-philippe.villemin/code/RNA-SEQ/",

To get an overview of the json, look into configs directory.

Here we show an example for the diff_exp.json :

![alt text](https://github.com/LucoLab/RNASEQ/blob/master/img/dif_exp_json.png "Json for Differential Expression")


## Launch differential expression

---

```shell
	python3 pathTo/diffGeneExp.py -c pathToConfigFile/diff_exp.json -p TestConditionName_vs_NormalConditionName
```

This script is a wrapper calling a Rscript called _diff_exp.R_. _diff_exp.R_ will use _Design.csv_ & _Raw_read_counts.csv_ created by the python wrapper using json configuration file.

_Design.csv_ & _Raw_read_counts.csv_ should be in _$path_to_output/output/$project_name/_ directory.

If you already have _Design.csv_ & _Raw_read_counts.csv_ , you can execute directly the Rscript as follows :

```R
Rscript ${PATH_TO_SCRIPT}/diff_exp.R  --dir ${PATH_TO_DATA}/[DIR_NAME] --cond1 [COND1]  --cond2 [COND2]  ${PATH_TO_DATA}/[DESIGN.csv] ${PATH_TO_DATA}/[GENE_READ_COUNT.csv] 
```

This is how _Design.csv_ should be :

![alt text](https://github.com/LucoLab/RNASEQ/blob/master/img/design.png "Design")

When you call the script, the p parameter need TestConditionName_vs_NormalConditionName to be set. It should be set in accordance with what you wrote in design.csv.

__Note :__ No need of last column.

This is how _Raw_read_counts.csv_ should be :

![alt text](https://github.com/LucoLab/RNASEQ/blob/master/img/raw_read.png "Raw_read")


Finally you get the following directories as output : 


![alt text](https://github.com/LucoLab/RNASEQ/blob/master/img/output_diffexp.png "Outputs")
)

