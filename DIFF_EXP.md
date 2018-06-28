
## Quick overview

---

Here we describe the **differential expression** part of the workflow.



## Set up Tools

---

First you need to check that the following tools are installed on server/computer.


Scripts available here are in Python3.  
It's not required but advised to install Conda if python3 is not set up on your computer.   
It will make things easier then for installing tools or switching to older python version if needed.

_**Conda**_ : [here](https://www.continuum.io/downloads)

DIFFERENTIAL EXPRESSION : 

_**R**_ : R [here](https://pbil.univ-lyon1.fr/CRAN/)

Packages : 


## Set up Files


```shell
conda install pyfaidx
faidx yourgenome.fasta -i chromsizes > yourgenome.genome.sizes
```

Also you need to create index for STAR :  

See config directory. You will find an example called paired.set1_align.json for a test dataset.


## Launch differential expression

---

```shell
	python3 pathTo/diffGeneExp.py -c pathToConfigFile/diff_exp.json -p TestConditionName_vs_NormalConditionName
```



You need to use this on each sample/replicate.

Finally you get the following directories as output : 


![alt text](https://github.com/LucoLab/RNASEQ/blob/master/img/output_alignment.png "Outputs")
)

