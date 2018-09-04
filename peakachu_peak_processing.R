#this script processes the output of the PEAKachu peakcaller

#this works with 2 replicate, not implemented for 1 replicate yet

######################
#check for packages
if("ggplot2" %in% rownames(installed.packages())) {
	cat(sprintf("\n ggplot2 installed \n\n"))
} else {
	cat(sprintf("\n ggplot2 package is missing. \n Installing it now \n\n"))
	source("http://bioconductor.org/biocLite.R")
	biocLite("ggplot2")
}

#if("xlsx" %in% rownames(installed.packages())) {
#	cat(sprintf("\n xlsx installed \n\n"))
#} else {
#	cat(sprintf("\n xlsx package is missing. \n Installing it now \n\n"))
#	install.packages("xlsx")
#}


library("ggplot2")
#library("xlsx")

######################
#set the FDR
fdr = as.numeric("0.05") # type a number between the brackets (no " ") to set the false discovery rate (default = 0.05) 
if ( is.na(fdr) ) fdr = 0.05  # default FDR

######################

setwd(".")

inputdir = file.path(getwd(), "peak_tables/")
outputdir = file.path(getwd(), "processed_peaks/")

if(!dir.exists(file.path(outputdir))) {
	dir.create(file.path(outputdir))
} else{
	cat("directory processed_peaks already exists...stopping here")
	break
}

csvfiles = list.files(inputdir, pattern = "*.csv")

csv.list = list()
#read the csv files into a list
for (i in 1:length(csvfiles)){
	csv.list[[i]]<-read.csv(file.path(inputdir, csvfiles[i]), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

csv.list = lapply(csv.list, function(x) {x[,c(1,3:11,15,2)]})

csv.merged = do.call(rbind,	csv.list)

#keep only those peaks from the main chromosomes
main_chr = c("2L", "2R", "3L", "3R", "X")
csv.merged = csv.merged[which(csv.merged[,1] %in% main_chr),]

write.table(csv.merged, file.path(outputdir, "PEAKachu_peaks_all.tsv"), sep="\t", quote=FALSE, col.names=NA)
#write.xlsx2(x = csv.merged, file = file.path(outputdir, "PEAKachu_peaks_all.xlsx"), row.names = FALSE, col.names = TRUE)

bed_all = csv.merged[,c(1:3,12,10,4)]
write.table(bed_all, file.path(outputdir, "PEAKachu_peaks_all.bed"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#padj plot
pdf(file.path(outputdir, "Padj_histogram.pdf"))
hist(csv.merged$total_G_padj_value, breaks=20, col="grey", main="Histogram of adjusted p-values", xlab="padj")
abline(v=fdr, col="red", lwd=1)
dev.off()

##MAplot
M = log2(rowMeans(csv.merged[,c(5,6)]) / rowMeans(csv.merged[,c(7,8)]))

A = 0.5 * (log2(csv.merged$base_means))

df = data.frame(A,M) #the rowname is the peakID

#just in case some values might be NAs or ifinite, remove them for plotting
df = df[which(!is.na(df$M)),]
df = df[which(!is.na(df$A)),]

df = df[which(is.finite(df$M)),]
df = df[which(is.finite(df$A)),]

colnames(df) = c("log2BaseMean", "log2FoldChange")
df$highlight = ifelse(csv.merged$total_G_padj_value < 0.05, TRUE, FALSE)
#rox2 = subset(df, rownames(df) == "FBgn0019660")
#rox1 = subset(df, rownames(df) == "FBgn0019661")

ggplot(df, ggplot2::aes(x = log2BaseMean, y = log2FoldChange, color = highlight)) +
	scale_color_hue(l=40, c=35) +
	geom_point(size = 0.5, alpha = 1) +
	geom_hline(color = "grey", yintercept = 0) +
	geom_vline(color = "grey", xintercept = 0) +
	theme_bw()
ggsave(file.path(outputdir, "MAplot.pdf"), width=7, height=6)
#dev.off()


padj.csv.merged = csv.merged[which(csv.merged$total_G_padj_value < fdr),]

padj.csv.merged = padj.csv.merged[order(padj.csv.merged$total_G_padj_value, decreasing = FALSE),]

write.table(padj.csv.merged, file.path(outputdir, "PEAKachu_peaks_padj.tsv"), sep="\t", quote=FALSE, col.names=NA)
#write.xlsx2(x = padj.csv.merged, file = file.path(outputdir, "PEAKachu_peaks_padj.xlsx"), row.names = FALSE, col.names = TRUE)

bed_padj = padj.csv.merged[,c(1:3,12,10,4)]
write.table(bed_padj, file.path(outputdir, "PEAKachu_peaks_padj.bed"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


