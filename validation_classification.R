# this function classifies the validation samples
validation_classification <- function(my_files, pam.centroids, output_dir, info, sample_info){
  
  # setting the output directory
  setwd(output_dir)
  
  lapply( my_files, function(x, pam_cent= pam.centroids){
    
    sg.nano.dat <- read.table(x, header = TRUE, sep = "\t" )
    sg.nano.dat.f <- apply(sg.nano.dat[, -1], 2, as.numeric )
    row.med.sg <- apply( sg.nano.dat.f, 1, median )
    sg.nano.dat.final <- sg.nano.dat.f - row.med.sg
    data2use <- data.frame( sg.nano.dat[, 1], sg.nano.dat.final)
    
    my_classification_result <- classify_validation_cohort( ucsc_pancreatic = data2use, my_pam_centroids = pam_cent)
    
    rownames(my_classification_result) <- gsub( "\\_.+|\\.CEL", "", rownames(my_classification_result) )
    
    x1 <- gsub(".+\\/", "", x)
    f2 <- gsub(".txt", "sd0_retraining_classifier_classification_results.txt", x1)
    
    write.table(my_classification_result, f2, sep = "\t", quote = FALSE )
    
  })
  
  
  ag_class_files <- list.files( path=output_dir, pattern = "classifier_classification_results.txt", full.names = TRUE, recursive = TRUE)
  agilent_samples_classification <- do.call("rbind", lapply( ag_class_files, read.table, header=TRUE ) )
  rownames(agilent_samples_classification) <- gsub( "\\_.+", "", rownames(agilent_samples_classification) )
  anno.Dat.all <- merge( sample_info, agilent_samples_classification[, c(7,11)], by.x = colnames(sample_info)[1], by.y = "row.names")
  anno.Dat <-  anno.Dat.all[, c(1,2,4)]
  # writing a rdata for the classification
  save( anno.Dat, file = paste0( rda_dir,"/",info ,"_sample_classification.rda") )
  
  # writing a rdata for the sample annotation
  # creating a proporpotion plot
  require(RColorBrewer)
  
  # considering the mixed samples
  
  create_subtype_prop_plot(df = anno.Dat, prop_pdf_file =paste0( output_dir,"/",info ,"_sample_classification_proportions.pdf")  )
  
  # writing the table
  write.table( table(anno.Dat$label, anno.Dat$final_sub), file = paste0( output_dir,"/",info ,"_sample_classification_proportions.txt"), sep = "\t", quote = FALSE) 
  
  # creating a pie chart
  produce_pie_chart_per_dataset(data_table = anno.Dat, plot_file = paste0( output_dir,"/",info ,"_sample_classification_piechart.pdf") )
  
  
    # removing the mixed samples the mixed samples
  
  anno.Dat.f <- anno.Dat[ !(anno.Dat$final_sub %in% c("ND")) , ]
  
  create_subtype_prop_plot(df = anno.Dat.f, prop_pdf_file =paste0( output_dir,"/",info ,"_removingND_sample_classification_proportions.pdf")  )
  
  # writing the table
  write.table( table(anno.Dat.f$label, anno.Dat.f$final_sub), file = paste0( output_dir,"/",info ,"_removingND_sample_classification_proportions.txt"), sep = "\t", quote = FALSE) 
  
  # creating a pie chart
  produce_pie_chart_per_dataset(data_table = anno.Dat.f, plot_file = paste0( output_dir,"/",info ,"_removingND_sample_classification_piechart.pdf") )
  
  
  
}

