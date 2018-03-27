

validation_data_silhoutte_width_analysis <- function(my_file, classification_data, info, rda_dir ){
  
  my_data <- read.delim(my_file, sep = "\t", header = TRUE, row.names=1)
  colnames(my_data) <- gsub("\\_.+|\\.CEL", "", colnames(my_data) )
  
  common_samples <- intersect( classification_data[,1], colnames( my_data) )
  
  classification_data <- classification_data[ classification_data[,1] %in% common_samples, ]
  classification_data <- classification_data[ order(classification_data[,3]), ]
  
  org_info <- unique( classification_data[,2])
  # selecting the samples
  gene_exp_data_num_sel <- subset( my_data, select = classification_data[,1])
  
  # now performing the sam analysis
  # performing the silhoutte analysis using the euclidean distance
  si=silhouette(as.numeric(factor(classification_data[,3])), dist(t(gene_exp_data_num_sel), "euclidean"))
  
  plot_file <- paste0(info,"_", org_info, "_silhoutte_width_plot.pdf")
  
  # converting the silhoutte result into matrix
  si_1_1=matrix(as.numeric(si), nrow=(ncol(gene_exp_data_num_sel)), ncol=3)
  si_1=data.frame( classification_data, si_1_1)
  
  # now saving the image for the silloute width
  save(si_1, file = paste0(rda_dir, "/",info,"_", org_info, "_silhoutte_width.rda") )
  
  sil2plot <- si_1[, c(1,3,5,6,7)]
  sil2plot[,2] <- gsub("\\-score|\\.score", "", sil2plot[,2] )
  sil2plot[,2] <- gsub("\\.", "\\-", sil2plot[,2] )
  
  sil_summary <- plot_silhoutte_results( silhoutte_data=sil2plot, plot_file=plot_file)
    
}