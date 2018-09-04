#various R functions 


import_many_files.function = function(files, path, separator, header, skip){
	peak.list = list()
	for (i in 1:length(files)){
		peak.list[[i]]<-read.csv(file.path(path, files[i]), sep = separator, header = header, skip = skip)
	}
	names(peak.list) = gsub("\\..*", "", files)
	return(peak.list)
}

gene_name_id_exchanger = function(ids, path) {
	ortho = read.table(path, header = TRUE, sep = "\t", na.strings=c("","NA")) 
	ortho = ortho[!duplicated(ortho[,1]),]
	rownames(ortho) = ortho[,1]
	ortho[,1] = NULL
	
	ortho_gene_names_dic = ortho[[1]]
	names(ortho_gene_names_dic) = rownames(ortho)
	
	d = data.frame(IDs=gsub("\\..*", "", ids), gene_names=NA)
	d$gene_names = ortho_gene_names_dic[ as.character(d$IDs) ]
}
