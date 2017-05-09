library(optparse)
library(data.table)
library(reshape2)
library(gplots)
library(heatmap3)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

#################################################################
#
# date: Ocotober 22, 2016
# platform: Ubuntu 16.04
# R.version : 3.2.2
# author: Villemin Jean-Philippe
# team: Epigenetic Component of Alternative Splicing - IGH
#
# merged_all_FC_files
# Usage : 
# merged_all_FC_files ${PATH_TO_SCRIPT}/merged_all_FC_files.R  --dir ${PATH_TO_DATA}/[DIR_OUTPUT]  
# 
# Merge in the same file all the fold changes obtained for the different comparisons.
# Path to files are hard coded in the code. Change it manually.
#
# Note :
# Need to sort each file by ensembl geneName column before executing this one
#################################################################

option_list = list(

  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="dir output", metavar="character")
); 

parser = OptionParser(usage = "%prog [options] file ",option_list=option_list);

arguments = parse_args(parser, positional_arguments = 0);

opt <- arguments$options
print(opt$dir)

# 2015/2016
#/mount_archive2/commun.luco/EMT_2015_SEQ_data/EMT_RNASEQ_FILES/RNASEQ_2016_RESULTS/MET.2015.2016/Expression14.11.16
#late_vs_early<- read.table(paste0(opt$dir,"DESEQ_all_res_annotated_sorted_pvalAdj_late_vs_early.csv",collapse="/"),sep=",", header=FALSE)
#early_vs_unT <- read.table(paste0(opt$dir,"DESEQ_all_res_annotated_sorted_pvalAdj_early_vs_unT.csv",collapse="/"),sep=",", header=FALSE)
#mant_vs_late<- read.table(paste0(opt$dir,"DESEQ_all_res_annotated_sorted_pvalAdj_mant_vs_late.csv",collapse="/"),sep=",", header=FALSE)
#late_vs_unT<- read.table(paste0(opt$dir,"DESEQ_all_res_annotated_sorted_pvalAdj_late_vs_unT.csv",collapse="/"),sep=",", header=FALSE)
#mant_vs_unT<- read.table(paste0(opt$dir,"DESEQ_all_res_annotated_sorted_pvalAdj_mant_vs_unT.csv",collapse="/"),sep=",", header=FALSE)
# 2017
#/mount_archive2/commun.luco/EMT_2015_SEQ_data/EMT_RNASEQ_FILES/RNASEQ_2017_RESULTS/MET.2017/Expression17.03.17
early_vs_unT  <- read.table(paste0(opt$dir,"DESEQ_all_res_annotated_sorted_pvalAdj_early_vs_unT.csv",collapse="/"),sep=",", header=FALSE)
late_vs_unT   <- read.table(paste0(opt$dir,"DESEQ_all_res_annotated_sorted_pvalAdj_late_vs_unT.csv",collapse="/"),sep=",", header=FALSE)
late_vs_early <- read.table(paste0(opt$dir,"DESEQ_all_res_annotated_sorted_pvalAdj_late_vs_early.csv",collapse="/"),sep=",", header=FALSE)

fc_dt_with_pval <- data.table(
            GENE=early_vs_unT[,1],
            CHROM=early_vs_unT[,2],
            START=early_vs_unT[,3],
            END=early_vs_unT[,4],
            STRAND=early_vs_unT[,5],
            BIOTYPE=early_vs_unT[,6],
           
           
            early_vs_unT=early_vs_unT[,13],
            early_vs_unT_padj=early_vs_unT[,12],
            late_vs_unT=late_vs_unT[,13],
            late_vs_unT_padj=late_vs_unT[,12],
            #mant_vs_unT=  mant_vs_unT[,13],  
            #mant_vs_unT_padj=  mant_vs_unT[,12],
            late_vs_early=late_vs_early[,13],
            late_vs_early_padj=late_vs_early[,12]
            #mant_vs_late =mant_vs_late[,13], 
            #mant_vs_late_padj =mant_vs_late[,12] 
           )

fc_dt_with_pval<- fc_dt_with_pval[rowSums(is.na(fc_dt_with_pval)) != ncol(fc_dt_with_pval)-6, ]

write.csv(fc_dt_with_pval,row.names=FALSE,file=paste0(c(opt$dir,"FC_collapse_with_pval.csv"),collapse="/"))

