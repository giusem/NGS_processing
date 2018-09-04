

list.of.packages <- c("DESeq2", "gplots", "ggplot2", "RColorBrewer", "xlsx", "SPIA", "geneplotter","ggrepel")
toInstall <- list.of.packages[which(!list.of.packages %in% rownames(installed.packages()))]

if(all( list.of.packages %in% rownames(installed.packages()) )) {
	cat(sprintf("\n All packages are already installed! \n\n"))
} else {
	cat(sprintf("\n Some required packages are missing. \n Installing them now \n\n"))
	source("http://bioconductor.org/biocLite.R")
	biocLite(toInstall)
}

library("DESeq2")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("xlsx")
library("SPIA")
library("geneplotter")
library("ggrepel") #extends ggplot labels to prettier version


## input

setwd("/data/akhtar/group/Giuseppe/smallRNAseq/03-20180319-Dis3-TF2S-kd/02-MERGEDreads/downstream_analysis/04-DESeq2/02-DESeq2/rnaseq_from_ttseq/")

sampleInfoFilePath = ("/data/akhtar/group/Giuseppe/smallRNAseq/03-20180319-Dis3-TF2S-kd/02-MERGEDreads/downstream_analysis/04-DESeq2/02-DESeq2/rnaseq_from_ttseq/sampleInfo.txt")
countFilePath = ("/data/akhtar/group/Giuseppe/smallRNAseq/03-20180319-Dis3-TF2S-kd/02-MERGEDreads/downstream_analysis/04-DESeq2/01-featureCounts/TTseq_fC/counts.txt")
geneNamesFilePath = ("/data/akhtar/group/Giuseppe/supplementary/biomart_gene_name/dm6_ensembl_release79_biomart.tsv")


fdr = as.numeric("0.05") # type a number between the brackets (no " ") to set the false discovery rate (default = 0.05) 
if ( is.na(fdr) ) fdr = 0.05  # default FDR

topN =as.numeric("100") # type a number between the brackets (no " ") to set the number of DEgenes to plot in the heatmap (default = 50)
if ( is.na(topN) ) topN = 100  # use topN genes for plot



## read input
sampleInfo = read.table(sampleInfoFilePath, header = TRUE)
sampleInfo = data.frame(as.data.frame(unclass(sampleInfo)))
sampleInfo = sampleInfo[order(sampleInfo$name, decreasing=F),]  # order by sample name

countdata = read.table(countFilePath, header = TRUE)
countdata = data.frame(countdata)

#check: if names of the setup table are subset of the count matrix column names
if ( ! all( is.element(sort(sampleInfo[,1]), sort(colnames(countdata))) ) ) {
	cat("Error! Count table column names and setup table names do NOT match!\n")
	print(as.character(sort(sampleInfo[,1])))
	print(sort(colnames(countdata)))
	break
	#quit(save = "no", status = 1, runLast = FALSE)   # Exit 1
}

#check: if names of the setup table are subset of the count matrix column names
if ( ! all( is.element(sort(sampleInfo[,1]), sort(colnames(countdata))) ) ) {
	cat("Error! Count table column names and setup table names do NOT match!\n")
	print(as.character(sampleInfo[,1]))
	print(colnames(countdata))
	break
	#quit(save = "no", status = 1, runLast = FALSE)   # Exit 1
}

countdata = countdata[,as.character(sampleInfo[,1])]

###########################################################################################
## Extra data collecting some measures for every sample (e.g. total
## counts, scaling factors,...)
info <- data.frame(row.names=sampleInfo$name)


###########################################################################################

###optional ERCC spike normalization (opt-in)

#ERCC spike normalization 

#spikecountFilePath = ("/data/akhtar/group/gaub/04_diffExpression/DESeq2_ERCC_spike_normalized/spike_counts.tsv")  # in this case the counts on the spikes

#sinfo <- read.delim(sampleInfoFilePath)
#spikecounts <- read.delim(spikecountFilePath)
#rnames<-spikecounts[,1]
#rownames(spikecounts) <- rnames
#spikecounts <- spikecounts[ , 2:ncol(spikecounts)]

#spikedata <- data.matrix(spikecounts)
#spikedata <- as.data.frame(spikedata)


#spikedata <- spikedata[,as.character(sampleInfo[,1])]
#head(spikedata)
#colnames(spikedata)


#dds_spike <- DESeqDataSetFromMatrix(
#  countData = spikedata,
#  colData = sinfo,
#  design = ~ condition)

#dds_spike <- DESeq2::estimateSizeFactors(dds_spike)
#dds_spike <- DESeq2::estimateDispersions(dds_spike)




## Initiate the counts DESeq dataset 
#dds = DESeqDataSetFromMatrix(
#  countData = countdata,
#  colData = sampleInfo,
#  design = ~ condition)
#design = ~ batch + condition)
#dds

#sizeFactors(dds) <- sizeFactors(dds_spike)

## end of spike normalization

###########################################################################################


#create DESeq dataset
dds = DESeqDataSetFromMatrix(
	countData = countdata,
	colData = sampleInfo,
	design = ~ condition)

## add counts per sample to info
info$total_counts = apply(assay(dds), 2, sum)  

## save count table to file
write.table(assay(dds),"DESeq2.counts.tsv", sep="\t", quote=FALSE, col.names=NA)

## DE analysis
assign("last.warning", NULL, envir = baseenv())
dds = DESeq(dds)
warnings()
sink("DESeq2.WARNINGS.txt"); warnings(); sink() # save warnings to file

## add size factors to info
info$size_factors = sizeFactors(dds) 

## DE analysis results
res = results(dds)

###########################################################################################


# save normalized counts to file
write.table(counts(dds, normalized=T),"DESeq2.counts_normalized.tsv", sep="\t", quote=FALSE, col.names=NA)

## dispersion plot
pdf("Fig01-dispersion_plot.pdf")
plotDispEsts(dds, main="Dispersion plot")
dev.off()

## count densities (for genes with non-zero count)
### get the genes that have counts > 0 in all samples
GeneCounts <- counts(dds)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
### plot count densities
pdf("Fig02-AllSamplesDensityPlot.pdf", width=6, height=6)
multidensity( counts(dds, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
dev.off()

###########################################################################################

## tranform gene names into a human readable version (if file available)


if (file.exists(geneNamesFilePath)) { 
	cat(paste("Gene names file found\n")) 
	#geneNames = read.csv(geneNamesFilePath, sep="\t", header=F, row.names=1, stringsAsFactors=FALSE)
	geneNames = read.csv(geneNamesFilePath, sep="\t", header=F, stringsAsFactors=FALSE)
	geneNames = geneNames[!duplicated(geneNames[,1]),]
	rownames(geneNames) = geneNames[,1]
	geneNames[,1] = NULL
	head(geneNames)
	
	if (length( intersect( gsub("\\..*", "", res@rownames), rownames(geneNames) ) ) > 0) {
		cat(paste("Names matching to IDs found\n")) 
		
		## make a dictionary
		gene_names_dic = geneNames[[1]]
		names(gene_names_dic) = rownames(geneNames)
		##gene_names_dic["ENSMUSG00000025332"]
	}
}

## generate a dataframe from ids
id_to_gene_name = function(ids) {
	d = data.frame(IDs=gsub("\\..*", "", ids), gene_names=NA)
	d$gene_names = gene_names_dic[ as.character(d$IDs) ]
	head(d)
	
	# some might be NAs; replace those by original ID
	d[which(is.na(d$gene_names)),]$gene_names = as.character(d[which(is.na(d$gene_names)),]$IDs)
	head(d)
	return(d$gene_name)
}


if ( exists("gene_names_dic") ) {
	cat("Gene names are available\n")
	# update res with gene names
	res@listData$gene_names = id_to_gene_name(rownames(res))
}
###########################################################################################

## write DE results into files

write.table(res,"DESeq2.results.tsv", sep="\t", quote=FALSE, col.names=NA)

de_total = res[which(res$padj < fdr),]
length(de_total[,1])
de_total = de_total[order(de_total$padj, decreasing = F),]
write.table(de_total,"DESeq2.de_all.tsv", sep="\t", quote=FALSE, col.names=NA)
write.xlsx2(x = de_total, file = "DESeq2.de_all.xlsx", row.names = TRUE, col.names = TRUE)

de_up = de_total[which(de_total$log2FoldChange>0),]
de_up = de_up[order(de_up$padj, decreasing=F),]   # order by adjusted p-value
length(de_up[,1])
write.table(de_up,"DESeq2.de_up.tsv", sep="\t", quote=FALSE, col.names=NA)
write.xlsx2(x = de_up, file = "DESeq2.de_up.xlsx", row.names = TRUE, col.names = TRUE)

de_down = de_total[which(de_total$log2FoldChange<0),]
de_down = de_down[order(de_down$padj, decreasing=F),]           # order by adjusted p-value
length(de_down[,1])
write.table(de_down,"DESeq2.de_down.tsv", sep="\t", quote=FALSE, col.names=NA)
write.xlsx2(x = de_down, file = "DESeq2.de_down.xlsx", row.names = TRUE, col.names = TRUE)

###########################################################################################

## plotting

## MA
pdf("Fig03-MA_plot.pdf", width=6, height=6)
par(mfrow=c(1,1))
plotMA(res, alpha=0.1, ylim=c(-3,3), 
			 main=sprintf("MA-plot\n(FDR: %.2f, up: %d, down: %d)",fdr,length(de_up[,1]),length(de_down[,1])),
			 ylab="log2 fold change")
dev.off()

## Volcano
plotVolcano <- function(res_obj, data=plot) {
	xlim = c(-3,3)
	ylim = c(0,20)
	cex=c(0.3,0.5)
	plotdata = data.frame(log2FoldChange=res_obj$log2FoldChange, padj=res_obj$padj )
	plotdata = plotdata[!is.na(plotdata),]
	plotdata$cex = cex[[1]]
	plotdata$pch = 16
	plotdata$col = "#525252"
	plotdata$col[plotdata$padj<=fdr] = "#cd0000"
	
	plotdata$pch[plotdata$log2FoldChange<xlim[[1]]] = 5
	plotdata$cex[plotdata$log2FoldChange<xlim[[1]]] = cex[[2]]
	plotdata$log2FoldChange[plotdata$log2FoldChange<xlim[[1]]] = xlim[[1]]
	
	plotdata$pch[plotdata$log2FoldChange>xlim[[2]]] = 5
	plotdata$cex[plotdata$log2FoldChange>xlim[[2]]] = cex[[2]]
	plotdata$log2FoldChange[plotdata$log2FoldChange>xlim[[2]]] = xlim[[2]]
	
	plotdata$pch[-log10(plotdata$padj) > ylim[[2]]] = 2
	plotdata$cex[-log10(plotdata$padj) > ylim[[2]]] = cex[[2]]
	plotdata$padj[-log10(plotdata$padj) > ylim[[2]]] = 10^-ylim[[2]]
	
	plot(plotdata$log2FoldChange, -log10(plotdata$padj),
			 main=sprintf("Volcano plot\n(FDR: %.2f, up: %d, down: %d)",fdr,length(de_up[,1]),length(de_down[,1])),
			 xlab="log2-fold change",
			 ylab="-log10 q-value",
			 xlim=xlim,
			 ylim=ylim,
			 cex=plotdata$cex, pch=plotdata$pch,
			 col=plotdata$col)
	abline(h=-log10(fdr), col=rgb(0,0,1,0.5), lwd=4)
	abline(v=0, col=rgb(1,0,0,0.5), lwd=4)
}


pdf("Fig04-Volcano_plot.pdf", width=6, height=6)
plotVolcano(res)
dev.off()

## padj histogram
pdf("Fig05-Padj_histogram.pdf")
hist(res$padj, breaks=20, col="grey", main="Histogram of adjusted p-values", xlab="padj")
abline(v=fdr, col="red", lwd=1)
dev.off()

## rlog transform counts
rld = rlog(dds)

## save rlog tranformed counts to file
write.table(assay(rld),"DESeq2.counts_rlog.tsv", sep="\t", quote=FALSE, col.names=NA)

## Sample distances
sampleDists = dist(t(assay(rld)))
## Euclidean sample distance heatmap
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = sprintf("%s", colnames(rld)) #paste(colnames(rld), rld$condition, sep="-")
colnames(sampleDistMatrix) = sprintf("%s", colnames(rld)) #paste(colnames(rld), rld$condition, sep="-")

colours = colorRampPalette(rev(brewer.pal(9, "GnBu")))(255)
pdf("Fig06-Distance_heatmap.pdf", width=6, height=6)
heatmap.2(sampleDistMatrix,trace="none",col=colours,
					main="Heatmap\n(Euclidean distances)",
					key = FALSE,
					notecex=0.2,
					cexRow=0.5, cexCol=0.5, margins=c(10,10),
					cellnote=round(sampleDistMatrix,1),
					notecol="black")
dev.off()

##PCA
data <- plotPCA(rld, intgroup=c("name", "condition"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=name, shape=condition)) +
	geom_hline(aes(yintercept=0, xintercept=0), colour="grey") +
	geom_vline(aes(vintercept=0, xintercept=0), colour="grey") +
	geom_point(size=3) +
	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
	theme_bw(base_size = 14) +
	theme(legend.position = "none") +
	ggtitle("PCA\n") +
	geom_text_repel(aes(label=name)) 
ggsave(file=sprintf("Fig07-PCA.pdf"), width=7, height=6)

##Top N genes plot
### topN genes by pvalue
d = data.frame(id=rownames(de_total), padj=de_total$padj)
if ( length(rownames(d)) < topN ) topN = length(rownames(d))

d_topx_padj = d[order(d$padj, decreasing=F),][1:topN,]

plotdata = assay(rld)[as.character(d_topx_padj$id),] 

### test
setdiff(as.character(d_topx_padj$id), rownames(plotdata))

if ( exists("gene_names_dic") ) rownames(plotdata) = id_to_gene_name(rownames(plotdata))  # exchange ids by gene names

#my_palette <- colorRampPalette(c("#0571b0", "#92c5de", "#fdb863", "#e66101"))(n = 20)
#colors = c(seq(-2,-1,length=100),seq(-0.999,0,length=100),seq(0.001,1,length=100), seq(1.001,2,length=100))
#pdf(sprintf("Fig08-gene_clustering_top%i_DE_genes.pdf",topN), pointsize = 9)
#heatmap.2(plotdata, scale="row", trace="none", dendrogram="column",
#         col=my_palette,
#        main=sprintf("Top %d DE genes (by p-value)", topN), keysize=1,
#       margins = c(8,10),
#      cexRow=0.7, cexCol=0.9, density.info="none")
#dev.off()

pdf(sprintf("Fig08-gene_clustering_top%i_DE_genes.pdf",topN), pointsize = 9)
heatmap.2(plotdata, scale="row", trace="none", dendrogram="both",
					col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
					main=sprintf("Top %d DE genes (by p-value)", topN),
					margins = c(8,10), keysize = 1,
					cexRow=0.7, cexCol=0.9, density.info="none")
dev.off()

pdf(sprintf("Fig09-gene_clustering_top%i_DE_genes_NOscale.pdf",topN), pointsize = 9)
heatmap.2(plotdata, trace="none", dendrogram="both",
					col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
					main=sprintf("Top %d DE genes (by p-value); scale = none", topN),
					margins = c(8,10),
					cexRow=0.7, cexCol=0.9, density.info="none")
dev.off()

###########################################################################################

write.table(info,"DESeq2.info.txt", sep="\t", quote=FALSE, col.names=NA)
sink("DESeq2.session_info.txt")
sessionInfo()
sink()

###########################################################################################
## counts of your favorite gene
msl3_counts = plotCounts(dds, gene="ENSG00000005302.13", intgroup =c("name", "condition"), returnData=TRUE)
ggplot(msl3_counts, aes(x=name, y=count, color=condition)) +
	geom_point(position=position_jitter(width=.1,height=0), size=2) +
	#geom_bar(stat="identity", aes(fill=condition)) +
	geom_text_repel(aes(label=round(count))) +
	geom_hline(aes(yintercept=0, xintercept=0), colour="grey") +
	geom_vline(aes(yintercept=0, xintercept=0), colour="grey") +
	ggtitle("Msl3 counts") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file=sprintf("Fig10-Msl3_counts.pdf"), width=13, height=6)

###########################################################################################
## SPIA

entrez = read.table("/data/akhtar/group/Giuseppe/14-Felicia_RNAseq_HDACinhibitors_Msl3patients/03-DESeq/ensemble_entrez.tsv", header = TRUE, sep = "\t", fill = TRUE)
entrez = entrez[which(!entrez$EntrezGeneID == "NA"),]
entrez = entrez[!duplicated(entrez$EntrezGeneID),]

df <- read.table("DESeq2.de_all.tsv", sep= "\t", header = TRUE, quote = "" )
df_rownames = gsub("\\..*", "", df$X)
df = data.frame(row.names = df_rownames, df[2:ncol(df)])

df_entrez = merge(x = df, y = entrez, by.x = "row.names", by.y ="EnsemblGeneID")
df_entrez = df_entrez[order(df_entrez$log2FoldChange, decreasing = TRUE), ]
df.map = df_entrez$log2FoldChange
names(df.map) = as.vector(df_entrez$EntrezGeneID)

countdata_colnames = as.data.frame(gsub("\\..*", "", row.names(countdata)))
colnames(countdata_colnames) = c("EnsemblGeneID")
allgenes_entrez= merge(x = entrez, y = countdata_colnames, by.x = "EnsemblGeneID", by.y = "EnsemblGeneID")
allgenes.map = as.character(allgenes_entrez$EntrezGeneID)
#source("http://bioconductor.org/biocLite.R")
#biocLite("SPIA")
#runSPIA

#A data frame containing the ranked pathways and various statistics: pSize is the number of genes
#on the pathway; NDE is the number of DE genes per pathway; tA is the observed total preturbation
#accumulation in the pathway; pNDE is the probability to observe at least NDE genes on the pathway
#using a hypergeometric model; pPERT is the probability to observe a total accumulation more extreme
#than tA only by chance; pG is the p-value obtained by combining pNDE and pPERT; pGFdr and
#pGFWER are the False Discovery Rate and respectively Bonferroni adjusted global p-values; and the
#Status gives the direction in which the pathway is perturbed (activated or inhibited). KEGGLINK
#gives a web link to the KEGG website that displays the pathway image with the differentially expressed
#genes highlighted in red.

spia.degenes = spia(df.map, allgenes.map, organism = "hsa", nB = 2000)
spia.degenes$Name = substr(spia.degenes$Name,1,20)

spia.degenes = spia.degenes[order(spia.degenes$pGFWER),]
plotP(spia.degenes,threshold=0.4)

all.spia = spia.degenes[c(1:5,7,9,11)]
#top10 pathways
top.spia = spia.degenes[1:10, c(1:5,7,9,11)]
colnames(top.spia) <- c("Name",colnames(top.spia[2:8]))

top.spia$pGFdr <- -log10(top.spia$pGFdr)


ggplot(top.spia,aes(Name,pGFdr,fill = Status, size = pSize, label = Name)) +
	geom_point(alpha = 0.7, shape = 21) +
	scale_size_area(max_size = 15) + theme_bw(base_size = 15) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
	scale_fill_manual(values = c("forestgreen","Red")) +
	labs(x = "Pathways", y = "-log10(p-value)",fill = "Status", title = "SPIA - Top 10 pathways")

ggsave(file=sprintf("Fig11-SPIA.pdf"), width=13, height=6)
