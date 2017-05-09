#################################################################
#
# date: April 14, 2017
# platform: Ubuntu 16.04
# R.version : 3.2.2
# author: Villemin Jean-Philippe
# team: Epigenetic Component of Alternative Splicing - IGH
#
# enrichment
# Usage : 
# Rscript ${PATH_TO_SCRIPT}/enrichment.R  --file ${PATH_TO_DATA_FILE}
# Plot enrichment for Go KEGG & INTERACTOME.
#
# Note :
# Work only on my laptop, on the labserver there is a bug with png lib...and I'm not root to upgrade anything.
#
#################################################################

library(optparse)
library(data.table)
library(reshape2)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(GenomicFeatures)
library(AnnotationDbi)
library(biomaRt)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(plyr)

option_list = list(
  make_option(c("-i", "--file"), type="character", default=NULL, help="dataset file name", metavar="character")
 
); 

parser = OptionParser(usage = "%prog [options] file ",option_list=option_list);

arguments = parse_args(parser, positional_arguments = 0);

opt <- arguments$options

input <- opt$file
print(input)
# hardcoded
dir_pathway <- dirname(input)


data<- read.table(input, sep=",", header=TRUE)

edb = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="jul2016.archive.ensembl.org")

gene_infos = getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype','chromosome_name','start_position','end_position','strand','entrezgene'),values=data[,1],filters='ensembl_gene_id',mart=edb)

# SET INPUT LIST
entrez_id  <- gene_infos$entrezgene
ensembl_id <- gene_infos$ensembl_gene_id

print(head(entrez_id))
print(head(ensembl_id))

####################################################################################
##########################      ENRICHMENT    ####################################
####################################################################################

# ##########################  KEGG ENRICHMENT #########################################

kegg <- enrichKEGG(entrez_id, organism="human", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2,use_internal_data=FALSE)

write.csv(as.data.frame (kegg),file=paste0(c(dir_pathway,"DESEQ_KEGG_ENRICHMENT.csv"),collapse="/"))

png(paste0(c(dir_pathway,"DESEQ_KEGG_DOTPLOT.png"),collapse="/"),width=1000,height=1000)
barplot(kegg, showCategory=30)
dev.off()

png(paste0(c(dir_pathway,"DESEQ_KEGG_ENRICHMAP.png"),collapse="/"),width=1000,height=1000)
try(enrichMap(kegg, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai))
dev.off()

# ##########################  GO ENRICHMENT #######################################
# 
# 

go_mf <- enrichGO(gene=entrez_id,OrgDb = org.Hs.eg.db,ont = "MF",pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

write.csv(as.data.frame (go_mf),file=paste0(c(dir_pathway,"DESEQ_GOMF_ENRICHMENT.csv"),collapse="/"))


png(paste0(c(dir_pathway,"DESEQ_GOMF_DOTPLOT.png"),collapse="/"),width=1000,height=1000)
barplot(go_mf, showCategory=30)
dev.off()


pdf(paste0(c(dir_pathway,"DESEQ_GOMF_GOGRAPH.pdf"),collapse="/"))
try(plotGOgraph(go_mf))
dev.off()

png(paste0(c(dir_pathway,"DESEQ_GOMF_ENRICHMAP.png"),collapse="/"),width=1000,height=1000)
try(enrichMap(go_mf, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai))
dev.off()

#bp2 <- simplify(go_mf, cutoff=1.2, by="p.adjust", select_fun=min)
#png("./results/DESEQ_GOMF_ENRICHMAP_SYMPLIFY.png",width=1000,height=1000)
#enrichMap(bp2,vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
#dev.off()

#go_mf <- setReadable(go_mf, OrgDb = org.Hs.eg.db,keytype="auto")
#cnetplot(go_mf,foldChange=log2_fold)


go_bp <- enrichGO(gene=entrez_id, OrgDb = org.Hs.eg.db,  ont = "BP", pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

write.csv(as.data.frame (go_bp),file=paste0(c(dir_pathway,"DESEQ_GO_BP_ENRICHMENT.csv"),collapse="/"))

png(paste0(c(dir_pathway,"DESEQ_GOBP_DOTPLOT.png"),collapse="/"),width=1000,height=1000)
barplot(go_bp, showCategory=30)
dev.off()

png(paste0(c(dir_pathway,"DESEQ_GOBP_ENRICHMAP.png"),collapse="/"),width=1000,height=1000)
try(enrichMap(go_bp, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai))
dev.off()

pdf(paste0(c(dir_pathway,"DESEQ_GOBP_GOGRAPH.pdf"),collapse="/"))
try(plotGOgraph(go_bp))
dev.off()

#go_bp <- setReadable(go_bp, OrgDb = org.Hs.eg.db,keytype="auto")
#png(paste0(c(dir_pathway,"DESEQ_GOBP_CNETPLOT.png"),collapse="/"))
#cnetplot(go_bp,foldChange=log2_fold)
#dev.off()


# ##########################  REACTOME ENRICHMENT #########################################

reactome <- enrichPathway(gene=entrez_id,organism = "human",pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2,readable=TRUE)

write.csv(as.data.frame (reactome),file=paste0(c(dir_pathway,"DESEQ_REACTOME_ENRICHMENT.csv"),collapse="/"))

png(paste0(c(dir_pathway,"DESEQ_REACTOME_DOTPLOT.png"),collapse="/"),width=1000,height=1000)
barplot(reactome, showCategory=30)
dev.off()

png(paste0(c(dir_pathway,"DESEQ_REACTOME_ENRICHMAP.png"),collapse="/"),width=1000,height=1000)
try(enrichMap(reactome, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai))
dev.off()

#reactome <- setReadable(reactome, OrgDb = org.Hs.eg.db, keytype="ENTREZID")



