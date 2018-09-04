# plotting the coverage of multiBigwig summary 



mBWs_plotter = function(multiBigwigSummary_output, samplenames){
	mbws_raw = read.table(multiBigwigSummary_output, sep = "\t", header = FALSE)
	colnames(mbws_raw) = c("chr", "start", "end", samplenames)
	
	keep_chrs = c("2L", "2R", "3L", "3R", "X")
	
	mbws_raw = mbws_raw[which(mbws_raw$chr %in% keep_chrs),]
	
	mbws_raw = mbws_raw[,c(1,4:ncol(mbws_raw))]
	
	plotdata.list = list()
	for (i in 1:length(samplenames)){
		plotdata = data.frame(chr = mbws_raw$chr, score = mbws_raw[i+1])
		plotdata = reshape2::melt(plotdata)

		plotdata.list[[i]] = ggplot2::ggplot(plotdata, ggplot2::aes(x=chr, y=value, fill = chr)) +
			ggplot2::geom_boxplot(notch = TRUE) +
			ggplot2::xlab("Chromosome") +
			ggplot2::ylab("Coverage") +
			ggplot2::ggtitle(paste0(" ", samplenames[i])) + 
			ggplot2::theme_classic()
	}
	
	#par(mfrow = c(round(length(samplenames)/2), round(length(samplenames)/2)))
	do.call(gridExtra::grid.arrange, c(plotdata.list, ncol = round(length(samplenames)/2)))
}

