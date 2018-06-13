#################################################################
#
# date: Ocotober 22, 2016
# platform: Ubuntu 16.04
# R.version : 3.2.2
# author: Villemin Jean-Philippe
# team: Epigenetic Component of Alternative Splicing - IGH
#
# diff_exp
# Usage : 
# diff_expRscript ${PATH_TO_SCRIPT}/diff_exp.R  --dir ${PATH_TO_DATA}/[DIR_NAME] --cond1 [COND1]  --cond2 [COND2]  ${PATH_TO_DATA}/[DESIGN.csv] ${PATH_TO_DATA}/[GENE_READ_COUNT.csv] > ${PATH_TO_DATA}/logR/[ANY_NAME].out
# 
# Find Differentially expressed gene.
# Output list of genes with Fold change values.
# Also create graphics (PCA,histogram,volcanoPlot) to verify dataset is reliable.
#
# Note :
# Last part which is doing GeneEnrichment has been removed because of missing library in the server used.
# Can be uncommented if you are root on your machine and want to generate gene enrichment graphics.
#
#################################################################

library(optparse)
####################################################################################
######################### Parameters  ##############################################
####################################################################################

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default="output", 
              help="directory name of output [default= %default]", metavar="character"),
  make_option(c("-cd1", "--cond1"), type="character", 
              help="first condition to use for comparison ", metavar="character"),
  make_option(c("-cd2", "--cond2"), type="character", 
              help="second condition to use for comparison ", metavar="character"),
  make_option(c("-e", "--exploration"),  action="store_true" , default = FALSE,
              help="Only Explore data...default False")
); 

#

parser = OptionParser(usage = "%prog [options] design_file matrice_file",option_list=option_list);

arguments = parse_args(parser, positional_arguments = 2);

opt <- arguments$options
file_design_sample  <- arguments$args[1]
file_matrice_sample <- arguments$args[2]

print(opt$dir)
print(opt$cond1)
print(opt$cond2)
print(file_design_sample)
print(file_matrice_sample)

####################################################################################
######################### Create Output Dir  #######################################
####################################################################################

if( !file.exists(file_design_sample) ) {stop(sprintf("Specified file ( %s ) does not exist", file_design_sample))} 
if( !file.exists(file_matrice_sample) ) {stop(sprintf("Specified file ( %s ) does not exist", file_matrice_sample))} 

dir_exploration = paste(opt$dir,"/exploration",sep="")
dir_pathway     = paste(opt$dir,"/pathway",sep="")
dir_final       = paste(opt$dir,"/final",sep="")
dir_scatter       = paste(dir_exploration,"/scatter",sep="")



#createDir <- ifelse(!dir.exists(file.path(opt$dir)), dir.create(file.path( dir_exploration),recursive = TRUE),FALSE)
createDir <- ifelse(!dir.exists(opt$dir), dir.create(file.path(dir_exploration),recursive = TRUE),FALSE)
if(!createDir){stop ("Dir already exist")}

dir.create(file.path(dir_scatter),recursive = TRUE)
dir.create(file.path(dir_pathway),recursive = TRUE)
dir.create(file.path(dir_final),recursive = TRUE)

####################################################################################
#####################################  Package Loading  ############################
####################################################################################

library(data.table)
library(reshape2)
library(edgeR)
library(DESeq2)
library(limma)
library(RColorBrewer)
library(gplots)
library(heatmap3)
library(grDevices)
library(genefilter)
library(ggplot2)
library(GenomicFeatures)
library(AnnotationDbi)
library(biomaRt)
library(stringr)
library(org.Hs.eg.db)
library(vsn)
library(plyr)
library(pheatmap)
library(PoiClaClu)
library(gtools)
library(clusterProfiler)
library(ReactomePA)

####################################################################################
#####################################  Package Loading  ############################
####################################################################################


#library(clusterProfiler)
#library(ReactomePA)
# COUNT DATA & SAMPLE DESIGN
counts_dt <- read.delim(file_matrice_sample,sep=",", stringsAsFactors = FALSE,header=TRUE)
sampleTable <- read.table(file_design_sample ,sep=",",header=TRUE)

counT = as.vector(colnames(counts_dt)[-1])
sampleT = as.vector(sampleTable$sample_id)
print("SAmple : ")
sampleT
print("Matrice : ")
counT
if(identical(sampleT,counT) == FALSE){stop("WARNING : TCHECK FOR SAMPLE ORDER IN MATRICE AND DESIGN SAMPLE FILES")}
print("Sample and Matrice Columns seems to be ok ! Go on... ")

# DATATABLE To DATAFRAME
setDF(sampleTable)

countdata <- counts_dt[,-(1)]
rownames(countdata) <- counts_dt[,1]
####################################################################################
#####################################  EdgeR  ######################################
####################################################################################

dge <- DGEList(counts=countdata)
dge <- calcNormFactors(dge)

png(paste0(c(dir_exploration,"EDGER_barplot_library_size.png"),collapse="/"))
par(mar=c( 8.1 ,4.1, 4.1 ,2.1))
barplot(dge$samples$lib.size,names=colnames(dge),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()


# that remove genes that have no counts, or only a single count across all samples:
# dge <- dge[rowSums(dge$counts) > 1,]
# 
# paste0(c("NUMBER OF GENES WITH AT LEAST ONE READ : ",nrow(dge)),collapse=NULL)

####################### Retrieve Gene Length For RPKM ##########################

# BIOMART OBJECT
# edb = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="jul2016.archive.ensembl.org")
# 
# dirty_ensembl_id <-rownames(dge)
# clean_ensembl_id <- str_match(dirty_ensembl_id,"^(\\w+)\\.[0-9]+")
# 
# gene_infos = getBM(attributes=c('ensembl_gene_id','start_position','end_position'),values=clean_ensembl_id[,2],filters='ensembl_gene_id',mart=edb)
# 
# gene_infos$start_position <- as.numeric(gene_infos$start_position)
# gene_infos$end_position   <- as.numeric(gene_infos$end_position)
# gene_infos$gene_length    <- gene_infos$end_position - gene_infos$start_position

####################### TRANSFORMATION ###############################
# 
# cpm   <- cpm(dge) 
# lrpkm <- rpkm(dge,gene.length=gene_infos[,3],log=TRUE)
# lcpm  <- cpm(dge, log=TRUE)
# 
# paste0(c("DGE LIBRARY SIZES: ",dge$lib.sizes),collapse=NULL)
# nb_samples <- length(colnames(dge))
# 
# # For each gene , set TRUE OR FALSE VECTOR
# keep.exprs <- rowSums(cpm>1)>= nb_samples
# 
# dge_filtered <-dge[keep.exprs, keep.lib.sizes=FALSE]
# 
# lcpm_filtered <- cpm(dge_filtered, log=TRUE)
# paste0(c("NUMBER OF GENES WITH CPM FILTRATION : ",dim(dge_filtered)),collapse=NULL)
# 
# ########################### BOX PLOT #############################
# 
# png(paste0(c(dir_exploration,"EDGER_boxplot_lcpm.png"),collapse="/"))
# par(mar=c( 8.1 ,4.1, 4.1 ,2.1))
# # Check distributions of samples using boxplots
# boxplot(lcpm, xlab="", ylab="Log2 counts per million",las=2)
# # Let's add a blue horizontal line that corresponds to the median logCPM
# abline(h=median(lcpm),col="blue")
# title("Boxplots of logCPMs (unnormalised)")
# dev.off()
# 
# 
# png(paste0(c(dir_exploration,"EDGER_boxlot_lrpkm.png"),collapse="/"))
# par(mar=c( 8.1 ,4.1, 4.1 ,2.1))
# # Check distributions of samples using boxplots
# boxplot(lrpkm, xlab="", ylab="Log2  reads per kilobase per million (log-RPKM)",las=2)
# # Let's add a blue horizontal line that corresponds to the median logCPM
# abline(h=median(lrpkm),col="blue")
# title("Boxplots of logRPKM (unnormalised)")
# dev.off()
# 
# ########################### DENSITY PLOT #############################
# 
# png(paste0(c(dir_exploration,"EDGER_DensityCurve_CPM.png"),collapse="/"))
# par(mar=c( 8.1 ,4.1, 4.1 ,2.1))
# nsamples <- ncol(dge)
# col <- brewer.pal(nsamples, "Paired")
# par(mfrow=c(1,3))
# 
# plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,main="", xlab="")
# title(main="B pre filtered", xlab="Log-cpm")
# abline(v=0, lty=3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col=col[i], lwd=2)
# }
# legend("topright", colnames(dge), text.col=col, bty="n")
# 
# 
# plot(density(lrpkm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
# title(main="B pre fltered", xlab="log rpkm")
# abline(v=0, lty=3)
# for (i in 2:nsamples){
#   den <- density(lrpkm[,i])
#   lines(den$x, den$y, col=col[i], lwd=2)
# }
# legend("topright", colnames(dge), text.col=col, bty="n")
# 
# plot(density(lcpm_filtered[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
#      main="", xlab="")
# title(main="C post fltered", xlab="log-cpm")
# abline(v=0, lty=3)
# for (i in 2:nsamples){
#   den <- density(lcpm_filtered[,i])
#   lines(den$x, den$y, col=col[i], lwd=2)
# }
# legend("topright", colnames(dge_filtered), text.col=col, bty="n")
# 
# 
# dev.off()

####################################################################################
############################# DESEQ ################################################
####################################################################################

#By removing the low count genes from
#the input to the FDR procedure, we can find more genes to be significant among those that we keep,
#and so improved the power of our test. This approach is known as independent filtering.

#Data transformation can be applied in different ways to compensate for the dispersion at low and high 
#expression ranges. These methods were developed after the first wave of RNASeq analysis 
#to try compensate for effects associated with this type of data.
#Remember that weakly expressed genes show much higher log fold changes than strongly expressed genes.
#This phenomenon, seen in most RNASeq datasets, is a direct consequence of dealing with count data, in which ratios are inherently noisier when counts are low.
#Therefore, it is useful to transform data to diminish this effect. 
#If data are used on the original count scale, the result will be dominated by highly expressed, highly variable genes; if logarithm-transformed data are used, undue weight will be given to weakly expressed genes, which show exaggerated LFCs, as discussed above. Therefore, DESeq2 implements a regularized logarithm transformation (rlog), which behaves similarly to a log2 transformation for genes with high counts, while shrinking the values for genes with low counts.
# Control for emt (replicate)

################################################################################################
# 2017
dds <- DESeqDataSetFromMatrix(countData=countdata,colData=sampleTable, design =~ condition)

# 2016/2015
#dds <- DESeqDataSetFromMatrix(countData=countdata,colData=sampleTable, design =~ emt + condition)
################################################################################################

#Note: it is prefered in R that the first level of a factor be the reference level (e.g. control, or untreated samples), so we can relevel the dex factor like so:
#dds$emt
dds$condition
#nrow(dds)
# Filter...could be more stringent
dds <- dds[rowSums(counts(dds)) > 1,]
print("Taill DDS")
dim(dds)
print("Taill DDS after  FITRAGE")
dim(dds)

dds <- DESeq(dds)
print("SIZEFACTORS: ")
sizeFactors(dds)

# When TRUE Unsupervised regularized-logarithm transformation (account for seqencing depth)
rld <- rlog(dds, blind=FALSE) 
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
fpm_sf <- fpm(dds, robust = TRUE)
fpm_rc <- fpm(dds, robust = FALSE)
print("assay fpm")
head(fpm_rc,5)
print("count fpm")
head(count(fpm_rc),5)
#write.csv(mcols(dds,use.names=TRUE),file=paste0(c(dir_final,"DEseq_OBJECT.csv"),collapse="/"))

#head(counts(dds))
#head(assay(rld))

###################  DESeq        TEST              #######################################
# alpha=.05 change after vennDiagramm series
res <- results(dds, contrast=c("condition",opt$cond1,opt$cond2), alpha=0.05)
print("Summary res < 0.05")
head(summary(res),10)
print("ResultsNames")
resultsNames(res)
print("MCols")
mcols(res, use.names = T)
print("TotGene<0.05:")
sum(res$pvalue < 0.05, na.rm=TRUE)
print("isNotNa:")
sum(!is.na(res$pvalue))
print("gene < padj 0.05")
sum(res$padj < 0.05, na.rm=TRUE)
################################### FILTERING #######################################

resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.05)
#FC 1.5
#
resEffectSig <- subset(resOrdered, ( abs(resOrdered$log2FoldChange) >= 1.5 & padj < 0.05 ) )
#resEffectSig <- subset(resOrdered, ( abs(logratio2foldchange(resOrdered$log2FoldChange, base=2) ) >= 1.5 & padj < 0.05 ) )
#Reorder to be sure
resEffectSig <- resEffectSig[order(resEffectSig$padj),]
#table(resEffectSig$padj < 0.05)
#SIMPLE COUNT TRUE FALSE
print("summary(resEffectSig)")
head(summary(resEffectSig),10)

################################### DENSITY PLOT #######################################
nsamples <- ncol((assay(rld)))
png(paste0(c(dir_exploration,"DESEQ_density.png"),collapse="/"),width=2000)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(log2(counts(dds, normalized=FALSE)+ 1)[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
title(main="A", xlab="Log2(x+1)")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(log2(counts(dds, normalized=FALSE)+ 1)[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames((assay(dds))), text.col=col, bty="n")

#plot(density(log2(counts(dds, normalized=TRUE)+ 1)[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
#title(main="B", xlab="Log2 normalised")
#abline(v=0, lty=3)
#for (i in 2:nsamples){
#  den <- density(log2(counts(dds, normalized=TRUE)+ 1)[,i])
#  lines(den$x, den$y, col=col[i], lwd=2)
#}
#legend("topright", colnames((assay(dds))), text.col=col, bty="n")

plot(density(assay(rld)[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
title(main="C", xlab="RLOG transformation")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(assay(rld)[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames((assay(rld))), text.col=col, bty="n")


#plot(density(log2(fpm_rc[,1]+1)), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
#title(main="C", xlab="log2(FPM+1) normalized by RC")
#abline(v=0, lty=3)
#for (i in 2:nsamples){
#  den <- density(log2(fpm_rc[,1]+1))
#  lines(den$x, den$y, col=col[i], lwd=2)
#}
#legend("topright", colnames(fpm_rc), text.col=col, bty="n")


#plot(density(log2(fpm_sf[,1]+1)), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
#title(main="C", xlab="log2(FPM+1) normalized by mean geometric * size factor")
#abline(v=0, lty=3)
#for (i in 2:nsamples){
#  den <- density(log2(fpm_sf[,i]+1))
#  lines(den$x, den$y, col=col[i], lwd=2)
#}
#legend("topright", colnames((fpm_sf)), text.col=col, bty="n")



# TODO
#dds <- dds[rowSums(counts(dds)) > 1,]

#log2(counts(dds, normalized=TRUE)+ 1)[,i]
# Filter...could be more stringent
#dds <- dds[rowSums(counts(dds)) > 1,]

# Unsupervised regularized-logarithm transformation (account for seqencing depth)
#rld <- rlog(dds, blind=FALSE) 
dev.off()

########################### BOXPLOT  GRAPHIC  ##################################
png(paste0(c(dir_exploration,"DESEQ_boxplot.png"),collapse="/"),width=1500)
par( mfrow = c( 1, 3 ),mar=c( 8.1 ,4.1, 4.1 ,2.1) )

boxplot(log2(counts(dds, normalized=FALSE)+ 1),las=2,col = c('red','red','sienna','sienna','palevioletred1','palevioletred1','orange','orange','yellow','yellow','darkolivegreen1','darkolivegreen4','forestgreen'))
                                                        
title("LOG2(x+1) ")

boxplot(log2(counts(dds, normalized=TRUE)+ 1),las=2,col = c('red','red','sienna','sienna','palevioletred1','palevioletred1','orange','orange','yellow','yellow','darkolivegreen1','darkolivegreen4','forestgreen'))
title("LOG2(x+1)  normalized")

#boxplot(log2(fpm_rc+1),las=2) # normalized
#title("log2(FPM+1) - Raw Count normalisation")

#boxplot(log2(fpm_sf+1),las=2) # normalized
#title("log2(FPM+1) -geometric mean x size factor normalisation")

boxplot(assay(rld),las=2,col = c('red','red','sienna','sienna','palevioletred1','palevioletred1','orange','orange','yellow','yellow','darkolivegreen1','darkolivegreen4','forestgreen')) # normalized
title("RLOG transformation")



dev.off()

########################### MEAN GRAPHIC  ##################################
png(paste0(c(dir_exploration,"DESEQ_disp_vs_mean.png"),collapse="/"))
plotDispEsts(dds)
dev.off()

pdf(paste0(c(dir_exploration,"DESEQ_sd_vs_mean.pdf"),collapse="/"))

notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds, normalized=FALSE)+ 1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))

dev.off()

#he log2 fold change for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor
#Genes with an adjusted p value below a threshold (here 0.1, the default
#This plot demonstrates that only genes with a large average normalized count contain sufficient information to yield a significant call
png(paste0(c(dir_final,"DESEQ_maplot.png"),collapse="/"))
plotMA(res, ylim=c(-5,5))
abline(h=1,col="dodgerblue",lwd=2)
abline(h=-1,col="dodgerblue",lwd=2)
dev.off()


########################### SCATTER PLOT SAMPLE VS SAMPLE  ##################################

# 1 - >mantrep1  2 -> mantrep2   3 ->  t1rep1    4 -> t1rep2    5 - >t6rep1   6 ->  t6rep2 7 -> unt6rep1 8-> unt6rep2
# 9 ->unt1rep1 10-> unt1rep2 11 -> t0rep1_2015 12 -> t0rep2_2015 13 -> t6rep1_2015 14 -> t6rep2_2015 15 -> t1rep1_2015 1 ->t1rep2_2015
plot_transformation <- function(index_1,index_2) {
  nameFile <- paste0(c(colnames(dds)[index_1],colnames(dds)[index_2],".png"),collapse="_")
  png(paste0(c(dir_scatter,nameFile),collapse="/"),width=5000)
  
  par( mfrow = c( 1, 3 ) )
  plot(counts(dds, normalized=FALSE)[,index_1:index_2] ,pch=16,ylim=c(0,25000),xlim=c(0,25000), cex=0.3)
  title("Raw Count")
  plot(log2(counts(dds, normalized=FALSE)[,index_1:index_2] + 1),pch=16, cex=0.3)
  title("LOG2(X+1)")
  #dds <- estimateSizeFactors(dds)
  #plot(log2(counts(dds, normalized=TRUE)[,index_1:index_2] + 1),pch=16, cex=0.3)
  #title("R-LOG2(X+1) Normalised")
  plot(assay(rld)[,index_1:index_2],pch=16, cex=0.3)
  title("RLOG(X) Transformation")
  dev.off()
}
#Scatterplot of transformed counts from two samples.
i <-1
nmax <- length(colnames(dds)) - 1
for (n in c(1:nmax) ){
  plot_transformation(n,n+1)
}

sampleDists <- dist( t( assay(rld) ) )

########################### HEATMAP   ##################################

# Heatmap of sample-to-sample distances using the rlog-transformed values.
png(paste0(c(dir_exploration,"HeatMap_Distance_euclidian.png"),collapse="/"))

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$sample_id )
colnames(sampleDistMatrix) <- paste(  rld$sample_id )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#Heatmap of sample-to-sample distances using the Poisson Distance
png(paste0(c(dir_exploration,"HeatMap_Distance_poisson.png"),collapse="/"))

poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$sample_id )
colnames(samplePoisDistMatrix) <- paste( rld$sample_id )
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)
dev.off()

########################### PCA  ##################################
pdf(paste0(c(dir_exploration,"PCA_MDS_plot.pdf"),collapse="/"))

#plotPCA(rld,  intgroup = c("emt"))
plotPCA(rld,  intgroup = c("condition"))
plotPCA(rld,  intgroup = c("sample_id"))

#data <- plotPCA(rld, intgroup = c( "emt", "condition"), returnData=TRUE)
data <- plotPCA(rld, intgroup = c( "condition"), returnData=TRUE)

percentVar <- round(100 * attr(data, "percentVar"))
#ggplot(data, aes(PC1, PC2, color=emt, shape=condition)) +   geom_point(size=3) +
ggplot(data, aes(PC1, PC2,  shape=condition)) +   geom_point(size=3) +

  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

########################### MDS  ##################################

#mdsData <- data.frame(cmdscale(sampleDistMatrix))
#mds <- cbind(mdsData, as.data.frame(colData(rld)))
#euclidianMDS <- ggplot(mds, aes(X1,X2,color=emt,shape=condition)) + geom_point(size=3)
#euclidianMDS <- ggplot(mds, aes(X1,X2,shape=condition)) + geom_point(size=3)

#euclidianMDS + ggtitle("MDS - Euclidian Distrance from rlog-transformed values")

#Creating the same plot for the PoissonDistance :
#mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
#mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds)))
#poissonMDS <- ggplot(mdsPois, aes(X1,X2,color=emt,shape=condition)) + geom_point(size=3)
#poissonMDS <- ggplot(mdsPois, aes(X1,X2,shape=condition)) + geom_point(size=3)
#poissonMDS + ggtitle("MDS - PoissonDistance values")

#dev.off()


#we test for genes that show significant effects of treatment 
#on gene counts more than doubling or less than halving, because 21 = 2.

if (opt$exploration){stop("EXLORATION TERMINATED !")}


########################### HIST and VOLCANO GRAPHIC  ##################################

#The Volcano plot arranges the DE genes along dimensions of biological and statistical significance. 
#The first (horizontal) dimension is the log fold change between the two groups, 
#and the second (vertical) axis represents the adjusted p-value 
#(on an inverse scale so smaller p-values appear higher up).
#The first axis indicates biological impact of the change
png(paste0(c(dir_final,"DESEQ_histogram_FC_padj_inf_005.png"),collapse="/"),width=1500)

# plot distribution of LFC in this subset
hist(resSig$log2FoldChange,
     xlim=c(-6,6),
     breaks=20,col="red",
     main="LFC distribution for genes with padj < 0.05")
dev.off()

#png(paste0(c(dir_final,"DESEQ_volcano_FC_padj_inf_005.png"),collapse="/"),width=1500)
#plot(resSig$log2FoldChange, 1-resSig$padj,
#     xlim=c(-6,6),
#     main="volcano plot for genes with padj < 0.05")
#dev.off()

########################### HIST and VOLCANO GRAPHIC  ##################################


#png(paste0(c(dir_final,"DESEQ_histogram_FC_padj_inf_005_fc_sup_1_5.png"),collapse="/"),width=1500)
#hist(resEffectSig$log2FoldChange,
#     xlim=c(-6,6),
#     breaks=20,
#     main="LFC distribution for genes with |LFC|>=1.5 & padj<0.05")
#dev.off()

#png(paste0(c(dir_final,"DESEQ_volcano_FC_padj_inf_005_fc_sup_1_5.png"),collapse="/"),width=1500)
#plot(resEffectSig$log2FoldChange, 1-resEffectSig$padj,
#     xlim=c(-6,6),ylim=c(0.90,1),
#     main="volcano plot for genes with |LFC|>=1.5 & padj<0.05")
#dev.off()

####################################################################################
##########################     ANNOTATION    #######################################
####################################################################################

# BIOMART OBJECT
edb = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="jul2016.archive.ensembl.org")

dirty_ensembl_id <-rownames(res)

clean_ensembl_id <- str_match(dirty_ensembl_id,"^(\\w+)\\.([0-9]+)")

# Retrieve gene infos And entrezeneId needed fr KEGGPATHWAY
gene_infos = getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype','chromosome_name','start_position','end_position','strand','entrezgene'),values=clean_ensembl_id[,2],filters='ensembl_gene_id',mart=edb)
res_data_frame <- as.data.frame(res)

# Convert rownames to column
tibble::rownames_to_column(res_data_frame)
# Set name for column gene
colnames(res_data_frame)[1] <- "ensembl_gene_id"
res_data_frame[,1] <- clean_ensembl_id[,2]

res_annotated <- join(gene_infos, res_data_frame, by='ensembl_gene_id', type='left', match='all')


######################## FILTER OUR RESULTS ##################################

res_annotated <- res_annotated[!duplicated(res_annotated$ensembl_gene_id),]

# Sorted by Adjusted
res_annotated <- res_annotated[order(res_annotated$padj),]

res_annotated_df <- as.data.frame(res_annotated)
rownames(res_annotated_df) <- NULL

res_annotated_df$foldchange <- logratio2foldchange(res_annotated_df$log2FoldChange, base=2)

res_annotated_df$classifier <- -log10(res_annotated_df$padj) * (res_annotated_df$log2FoldChange/abs(res_annotated_df$log2FoldChange))
res_annotated_df$rank <-  rank(res_annotated_df$classifier)

#data.frame(append(res_annotated_df, list(FC=Foldchange), after=match("log2FoldChange", names(res_annotated_df))))
res_annotated_df_sorted <- res_annotated_df[order(res_annotated_df$classifier ),decreasing = TRUE),]
write.csv(res_annotated_df_sorted ,row.names=FALSE,file=paste0(c(dir_final,"DESEQ_all_res_annotated_sorted_pvalAdj.csv"),collapse="/"))

#write.csv(res_annotated_df ,row.names=FALSE,file=paste0(c(dir_final,"DESEQ_all_res_annotated_sorted_pvalAdj.csv"),collapse="/"))

res_annotated_filtered <- subset(res_annotated_df, ( abs(res_annotated_df$log2FoldChange)>= 1.5 & padj < 0.05 ))
res_annotated_filtered <- res_annotated_filtered[order(abs(res_annotated_filtered$log2FoldChange),decreasing = TRUE),]
write.csv(res_annotated_filtered,row.names=FALSE,file=paste0(c(dir_final,"DESEQ_res_annotated_DE.csv"),collapse="/"))

res_annotated_filtered_up <- subset(res_annotated_df, ( res_annotated_df$log2FoldChange >= 1.5 & padj < 0.05 ))
res_annotated_filtered_up <- res_annotated_filtered_up[order(abs(res_annotated_filtered_up$log2FoldChange),decreasing = TRUE),]
write.csv(res_annotated_filtered_up,row.names=FALSE,file=paste0(c(dir_final,"DESEQ_res_annotated_DE_UP.csv"),collapse="/"))

res_annotated_filtered_down <- subset(res_annotated_df, ( res_annotated_df$log2FoldChange <= -1.5 & padj < 0.05 ))
res_annotated_filtered_down <- res_annotated_filtered_down[order(abs(res_annotated_filtered_down$log2FoldChange),decreasing = TRUE),]
write.csv(res_annotated_filtered_down,row.names=FALSE,file=paste0(c(dir_final,"DESEQ_res_annotated_DE_DOWN.csv"),collapse="/"))
