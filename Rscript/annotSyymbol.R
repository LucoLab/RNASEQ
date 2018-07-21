#################################################################
#
# date: Ocotober 22, 2016
# platform: Ubuntu 16.04
# R.version : 3.2.2
# author: Villemin Jean-Philippe
# team: Epigenetic Component of Alternative Splicing - IGH
#
# annotFileWithSymbol
# Usage : 
# 
# Annotation of whippet output files.
# Pass a file called kikou
# Output is kikou.annotated.csv
#
#################################################################

####################################################################################
#####################################  Package Loading  ############################
####################################################################################
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(GenomicFeatures)))
suppressWarnings(suppressMessages(library(AnnotationDbi)))
suppressWarnings(suppressMessages(library(biomaRt)))
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(optparse)))
####################################################################################
######################### Parameters  ##############################################
####################################################################################
#Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/annotSyymbol.R --file=/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCm38_PRIM_GENCODE_M15/gencode.v15.gene.length.exon.sums.csv --organism=human

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="Absolute file input path", metavar="character"),
  make_option(c("-o", "--organism"), type="character", default=NULL, help="organism", metavar="character")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

######################################################
#########         MAIN                     ###########
######################################################

#' Trick to handle if script is launched directly or by generateDoc
#+ opt-param, include=FALSE
if (length(opt) != 3 ) { 
  opt$file <- params$file
  opt$organism <- params$organism

}


# COUNT DATA & SAMPLE DESIGN
sampleTable <- try(read.table(opt$file ,sep="\t",header=TRUE), silent = TRUE)

if (!inherits(sampleTable, 'try-error')){ 

organism = ""
host     = ""
symnol_description=''

if (opt$organism=="mouse") {
	print("Organism : Mouse")
	organism ="mmusculus_gene_ensembl"
	host="aug2017.archive.ensembl.org"
	symbol_description='mgi_symbol'

}
if (opt$organism=="human") {
	print("Organism : Human")

	organism = "hsapiens_gene_ensembl"
	host="jul2016.archive.ensembl.org"
	symbol_description='hgnc_symbol'

}

# BIOMART OBJECT
edb = useMart("ENSEMBL_MART_ENSEMBL", dataset=organism,host=host)
colnames(sampleTable)[colnames(sampleTable) == 'gene'] <- 'ensembl_gene_id'

# Retrieve gene infos And entrezeneId needed fr KEGGPATHWAY
gene_infos = getBM(attributes=c('ensembl_gene_id',symbol_description,'gene_biotype'),values=sampleTable$ensembl_gene_id,filters='ensembl_gene_id',mart=edb)

res_annotated <- join(gene_infos, sampleTable, by='ensembl_gene_id', type='left', match='all')

if (opt$organism=="human") {
	res_annotated$hgnc_symbol <- ifelse(res_annotated$hgnc_symbol == "", res_annotated$ensembl_gene_id, res_annotated$hgnc_symbol )
}
if (opt$organism=="mouse") {
	res_annotated$mgi_symbol <- ifelse(res_annotated$mgi_symbol == "", res_annotated$ensembl_gene_id, res_annotated$mgi_symbol )
}

filewithoutExtension= unlist(strsplit(opt$file, split='.csv', fixed=TRUE))[1]

output=paste(c(filewithoutExtension,"annoted.csv"),collapse="")
output

write.csv(res_annotated,row.names=FALSE,file=output,quote = FALSE,sep="\t")

 }