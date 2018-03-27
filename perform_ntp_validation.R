# function to perform the kras dependency
perform_ntp_validation <- function( gene.exp.file, output.fol){
  
  load("/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Collection/2018/Multiorgan/rdata/sample_info.rda")
  
  f2 <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Collection/2018/Multiorgan/compare_metagenes/MOT_SD0.8_Pancan_SD1.3_metagenes_comparision_commongenes_ntp_marker.txt"  
  
  ff <-  gsub( ".txt","" , basename(gene.exp.file) )
  
  
  gene_dat <- read.delim2(gene.exp.file, sep = "\t", header = TRUE)
  
  
  
  changeNameTCGA <- function(sample.names){
    syn_dat_ID <- do.call( "rbind", strsplit( sample.names, "\\."))[,1:3]
    # creating the modified TCGA name
    xx <- paste(syn_dat_ID[,1], syn_dat_ID[,2],syn_dat_ID[,3],  sep=".")
    return(xx) 
    
  }
  
  
  # for the ICGC data
  
  ind_icgc <- grep("ICGC", colnames(gene_dat))
  
  if(length(ind_icgc) ==0){
    
  
    ind_tcga <- grep("TCGA", colnames(gene_dat))
    if(length(ind_tcga) >=1) colnames(gene_dat) <- changeNameTCGA(colnames(gene_dat) )
  
  
    rm(ind_tcga)
  
    ind_tcga <- grep("TCGA", sample_info$geo_accession)
  
    if(length(ind_tcga) >=1) sample_info$geo_accession[ind_tcga] <- changeNameTCGA(sample_info$geo_accession[ind_tcga])
  
    rm(ind_tcga)
    ind_tcga <- grep("TCGA", sample_info$geo_accession)
  
    if(length(ind_tcga) >=1) sample_info$geo_accession[ind_tcga] <- changeNameTCGA(sample_info$geo_accession[ind_tcga])
  
    # doing the same for the geo datasets
    ind_gsm <- grep("GSM|\\_", colnames(gene_dat))
    if(length(ind_gsm) >=1) colnames(gene_dat) <- gsub("\\_.+|\\.CEL", "", colnames(gene_dat) , ignore.case = TRUE)
    
    # doing the same for the geo datasets
    ind_gsm <- grep("GSM|\\_", sample_info$geo_accession)
    if(length(ind_gsm) >=1) sample_info$geo_accession[ind_gsm] <- gsub("\\_.+|\\.CEL", "", sample_info$geo_accession[ind_gsm], ignore.case = TRUE )
    
    
  }
  
  
  sample_info <- sample_info[sample_info$geo_accession %in% colnames(gene_dat),  ]
  
  
  org <- unique( sample_info[,2] )
  
  out_fol <- paste0(output.fol,"/" ,org)
  
  dir.create(out_fol, recursive = TRUE, showWarnings = FALSE)
  #f3 <- paste0(out_fol, "/", ff, "_mixed_samples_removed.txt")
  
  setwd(out_fol)
  # running the NTP on the mixed samples removed data file
  NTPez(gene.exp.file, f2)
  
  out_data_f <- list.files( path=out_fol, pattern="NTP_prediction_result.xls", full.names = TRUE)
  out_data <- read.delim( out_data_f, sep = "\t", header = TRUE)
  
  rm(ind_gsm)
  
  ind_icgc <- grep("ICGC", out_data[,1])
  
  if(length(ind_icgc) ==0){
    
    
    ind_tcga <- grep("TCGA", out_data[,1] )
    if(length(ind_tcga) >=1) out_data[,1] <- changeNameTCGA(out_data[,1] )
    
    
    rm(ind_tcga)
    
  # doing the same for the geo datasets
  ind_gsm <- grep("GSM|\\_", out_data$sample.names )
  if(length(ind_gsm) >=1) out_data$sample.names[ind_gsm] <- gsub("\\_.+|\\.CEL", "", out_data$sample.names[ind_gsm], ignore.case = TRUE )
  
  }
  
  if( length( unique(out_data$predict.label)) ==6 ){ out_data$subtypes <- factor( out_data$predict.label, labels  = c("Basal", "Classical",  "Inflammatory", "SP", "Stem-like", "TA"))
  out_data$subtypes <- as.character(out_data$subtypes)
  }else{
    
    out_data$subtypes <- out_data$predict.label
    out_data$subtypes <- gsub("1", "Basal", out_data$predict.label)
    out_data$subtypes <- gsub("2", "Classical", out_data$predict.label)
    out_data$subtypes <- gsub("3", "Inflammatory", out_data$predict.label)
    out_data$subtypes <- gsub("4", "SP", out_data$predict.label)
    out_data$subtypes <- gsub("5", "Stem-like", out_data$predict.label)
    out_data$subtypes <- gsub("6", "TA", out_data$predict.label)
  }
  ind_sel <- which( out_data$BH.FDR > 0.2)
  out_data$predict.label[ind_sel] <- "ND"
  out_data$subtypes[ind_sel] <- "ND"
  
  b <- table( out_data$predict.label, out_data$subtypes )
  rownames(b) <- paste0("NTP_labels_", rownames(b))
  
  tag <- gsub( ".+\\/|\\.txt", "", gene.exp.file )
  
  write.table(b, paste0( tag, "-NTP-prediction_cutoff_FDR0.2.counts.txt"), sep = "\t")
  
  
  data2plot_pie <- merge( out_data, sample_info, by.x = colnames(out_data)[1], by.y = colnames(sample_info)[1] )
  
  
  classification_file <- paste0( tag, "-NTP-prediction_cutoff_FDR0.2.classification.txt")
  
  write.table(data2plot_pie, classification_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
 # f6 <- paste0( tag, "-NTP-prediction_cutoff_FDR0.2.proportions.pdf")
#  produce_pie_chart_per_dataset(data_table = data2plot_pie[, c(1,9,8)], plot_file = f6)
  
  
}
