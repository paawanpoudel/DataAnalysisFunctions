# changing the name for the samples in TCGA
changeNameTCGA_miR <- function(sample.names){
  
  sample.names <- gsub("\\-","\\.", sample.names )
  syn_dat_ID <- do.call( "rbind", strsplit( sample.names, "\\."))[,1:5]
  # creating the modified TCGA name
  xx <- paste(syn_dat_ID[,1], syn_dat_ID[,2],syn_dat_ID[,3],syn_dat_ID[,4], syn_dat_ID[,5], sep=".")
  return(xx) 
  
}



changeNameTCGA_rppa <- function(sample.names){
  
  sample.names <- gsub("\\-","\\.", sample.names )
  syn_dat_ID <- do.call( "rbind", strsplit( sample.names, "\\."))[,1:4]
  # creating the modified TCGA name
  xx <- paste(syn_dat_ID[,1], syn_dat_ID[,2],syn_dat_ID[,3],syn_dat_ID[,4], sep=".")
  return(xx) 
  
}




changeNameTCGA_mutation <- function(sample.names){
  
  sample.names <- gsub("\\-","\\.", sample.names )
  syn_dat_ID <- do.call( "rbind", strsplit( sample.names, "\\."))[,1:3]
  # creating the modified TCGA name
  xx <- paste(syn_dat_ID[,1], syn_dat_ID[,2],syn_dat_ID[,3], sep=".")
  return(xx) 
  
}
