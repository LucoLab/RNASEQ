library(optparse)
library(data.table)
library(readr)

#################################################################
#
# date: May 10, 2016
# platform: Ubuntu 16.04
# R.version : 3.2.2
# author: Villemin Jean-Philippe
# team: Epigenetic Component of Alternative Splicing - IGH
#
# filter
# Usage : 
# RScript ${PATH_TO_SCRIPT}/filer  --file1 ${PATH_TO_FILE}  --file2 ${PATH_TO_FILE} -p PREFIX
# 
# FIlter out from DEG  list genes of interest 
# 
#################################################################
#Rscript filter.R --file1 /home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/ZinFinger.tsv --file2 /home/jean-philippe.villemin/RNASEQ_2017_RESULTS/9_05_2017_early_treated_vs_unT/final/DESEQ_res_annotated_DE.csv -p ZNF
#Rscript filter.R --file1 /home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/ZinFinger.tsv  --file2 /home/jean-philippe.villemin/RNASEQ_2017_RESULTS/9_05_2017_late_treated_vs_unT/final/DESEQ_res_annotated_DE.csv -p ZNF
#Rscript filter.R --file1 /home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/ZinFinger.tsv --file2 /home/jean-philippe.villemin/RNASEQ_2016_RESULTS/9_05_2016_early_treated_vs_unT/final/DESEQ_res_annotated_DE.csv -p ZNF
#Rscript filter.R --file1 /home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/ZinFinger.tsv  --file2 /home/jean-philippe.villemin/RNASEQ_2016_RESULTS/9_05_2016_late_treated_vs_unT/final/DESEQ_res_annotated_DE.csv -p ZNF
#

option_list = list(
  
  make_option(c("-f1", "--file1"), type="character", default=NULL, 
              help="list file", metavar="character"),
  make_option(c("-f2", "--file2"), type="character", default=NULL, 
              help="input file you want to extract occurencies of your list", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix your ouput file", metavar="character")
); 

parser = OptionParser(usage = "%prog [options] file ",option_list=option_list);

arguments = parse_args(parser, positional_arguments = 0);

opt <- arguments$options
print(opt$file1)
print(opt$file2)

dir_output <- dirname(opt$file2)


data_list <- fread(opt$file1)
names(data_list)[2] <- "SYMBOL"
gene_of_interest <- data_list$SYMBOL

data_input <- fread(opt$file2)
gene <- c(data_input$hgnc_symbol)

gene_found <- intersect(gene_of_interest,gene)
gene_found_not_found <- setdiff(gene,gene_of_interest)

data_list_filtered<- data_input[which(data_input$hgnc_symbol %in% gene_found==TRUE)]
head(gene_found_not_found)
head(data_list_filtered)

name <- paste(opt$prefix,basename(opt$file2))
name_notFOUND <- paste("NOTFOUND_",name)

write.csv(data_list_filtered, paste0(c(dir_output,name),collapse="/"),row.names=F)
write.csv(gene_found_not_found, paste0(c(dir_output,name_notFOUND),collapse="/"),row.names=F)
