library(GenomicFeatures)
library(optparse)
library(rtracklayer)
library(stringr)
library(reshape2)
library(plyr) 

# Based on an idea of https://www.biostars.org/p/83901/
#Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/geneSize.R --file /home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/gencode.v25.primary_assembly.annotation.gtf

option_list = list(
  
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="file input", metavar="character")
); 

parser = OptionParser(usage = "%prog [options] file ",option_list=option_list);

arguments = parse_args(parser, positional_arguments = 0);

opt <- arguments$options
print(opt$file)

# First, import the GTF-file that you have also used as input for htseq-count
txdb <-makeTxDbFromGFF(opt$file,format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})

exon.genomics <- lapply(exons.list.per.gene,function(x){reduce(x)})

head(exon.genomics)

df <- do.call("rbind", lapply(exon.genomics, as.data.frame)) 

write.csv(df, file = "exons.non-redundant.csv")
write.csv(exonic.gene.sizes , file = "gene.length.exon.sums.csv")

# need to be transform in bed 
# remove things like this KI270728.1
# need to be sort 

