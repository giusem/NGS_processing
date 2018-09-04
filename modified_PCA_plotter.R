#modified plotPCA function

plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
	rv <- rowVars(assay(object))
	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
																										 length(rv)))]
	pca <- prcomp(t(assay(object)[select, ]))
	percentVar <- pca$sdev^2/sum(pca$sdev^2)
	if (!all(intgroup %in% names(colData(object)))) {
		stop("the argument 'intgroup' should specify columns of colData(dds)")
	}
	intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
	group <- if (length(intgroup) > 1) {
		factor(apply(intgroup.df, 1, paste, collapse = " : "))
	}
	else {
		colData(object)[[intgroup]]
	}
	d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
									intgroup.df, name = colData(rld)[,1])
	if (returnData) {
		attr(d, "percentVar") <- percentVar[2:3]
		return(d)
	}
	ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
	
}


data <- plotPCA.san(rld, intgroup=c("name", "condition"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC2, PC3, color=name, shape=condition)) +
	geom_hline(aes(yintercept=0, xintercept=0), colour="grey") +
	geom_vline(aes(vintercept=0, xintercept=0), colour="grey") +
	geom_point(size=3) +
	xlab(paste0("PC2: ", percentVar[1], "% variance")) +
	ylab(paste0("PC3: ", percentVar[2], "% variance")) +
	theme_bw(base_size = 14) +
	theme(legend.position = "none") +
	ggtitle("PCA\n") +
	geom_text_repel(aes(label=name)) 

ggsave(file=sprintf("Fig07-PCA.pdf"), width=7, height=6)