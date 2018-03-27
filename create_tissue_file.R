create_tissue_file <- function(){
  
  require(reshape2)
  
  tissue_file1 <- "/Users/ppoudel/Dropbox (SPM)//Pawan/ANALYSIS/Internal/Collection/2018/Multiorgan/tiger/modified_files//hs2tissue-Table_mod.txt"
  tissue_file2 <- "//Users/ppoudel/Dropbox (SPM)//Pawan/ANALYSIS/Internal/Collection/2018/Multiorgan/tiger/modified_files//symbol2hs-Table_mod.txt"
  
  tsg1 <- read.delim(tissue_file1, sep = "\t", header = FALSE)[-1, ]
  tsg2 <- read.delim(tissue_file2, sep = "\t", header = FALSE)[-1, ]
  
  
  tsg1.1 <-  melt(tsg1, id.vars = c("V1"))
  colnames(tsg1.1) <- c("unigene", "var", "organ")
  
  tsg1.2 <-  melt(tsg2, id.vars = c("V1"))
  colnames(tsg1.2) <- c("hgnc", "var", "unigene")
  
  # now merging the 2 columns
  my_tiger_genes <- merge( tsg1.1[,c(1,3)], tsg1.2[, c(1,3)], by.x="unigene", by.y="unigene")
  
  # remove the columns with no organs
  my_tiger_genes <- my_tiger_genes[ which(my_tiger_genes$organ !=""), ]
  
  return(my_tiger_genes)
  
}

# creating the histogram for the sd value files for each cancer type
create_sd_histogram <- function(file) {
  
  hist_file <- gsub(".txt", "SD_histogram.pdf", file)
  gene_sd_file <- gsub(".txt", "gene_SD.txt", file)
  
  sd_data <- read.delim2(file, sep="\t", header=TRUE)
  sd_data_1 <- apply(sd_data[,-1], 2, as.numeric)
  my_sd <- apply(sd_data_1, 1, sd)
  
  # genes and sd values
  my_gene_sd <- data.frame(sd_data[,1], my_sd)
  
  df <- melt(my_sd)
  # using the ggplot to creat the histogram
  hh <- ggplot(data=df, aes(value)) +
    geom_histogram(breaks=seq(min(my_sd), max(my_sd), by = 0.1),
                   fill="grey",
                   col="black",
                   alpha = .2) +
    labs(title="Histogram for Standard deviation") +
    labs(x="Standard deviation", y="Number of genes") + scale_fill_Publication()+ theme_Publication()
  
  # writing the
  ggsave(hh, file=hist_file)
  
  org<- gsub(".+\\/","", file)
  org<- gsub(".+\\/\\/","", file)
  org<- gsub("\\.txt","", file)
  
  # write the gene list to the file
  colnames(my_gene_sd) <- c("Genes", org)
  write.table( my_gene_sd,gene_sd_file, sep="\t", quote = FALSE)
  
}

get_sd_genes <- function(sd_value, sd_data){
  
  my_sd <- apply(sd_data[,-1], 1, min)
  # selecting the sd and comparing with the tissue specific genes
  ind_sd <- which(my_sd >= sd_value)
  # get the genes in sd 0.8
  genes_sd <- as.character(sd_data[ind_sd,1] )
  
  return(genes_sd)
  # after sd selections
  
  
}

# comparing the SD genes with the tissue specific genes
compare_genes <- function(sd_genes, tissue_specific_genes){
  
  m1 <- match(as.character( tissue_specific_genes[,1]), as.character(sd_genes ))
  w1 <- which(!is.na(m1))
  gene_sel <- sd_genes[m1[w1]]
  
  return(gene_sel)
}


# this function creates the plot for the SD versus tissue specific genes
create_tissue_sd_plots <- function(sd0_files, output_dir, info){
  
  lapply(sd0_files, create_sd_histogram)
  
  if(info=="mot") gene_sd_file <- list.files(path=dirname(sd0_files)[1], pattern="gene_SD", full.names = TRUE)
  if(info == "pancancer") gene_sd_file <- list.files(path=gsub("all_genes.+", "all_genes\\/" ,dirname(sd0_files)[1] ), pattern="gene_SD", full.names = TRUE, recursive = TRUE)
  
  merged_data <-  lapply(gene_sd_file, read.table, sep="\t" )
  
  # combining the sd values from each organ-tumor type to the genes
  sd_merged_matrix <- Reduce(function(x,y) {merge(x,y)}, merged_data)
  
  # sd values for the analysis
  my_sd_values<- c(0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5)
  my_sd_genes <- lapply(as.list(my_sd_values), get_sd_genes, sd_merged_matrix)
  
  if(info=="mot")  mot_tissues <- c("colon", "lung", "mammary_gland", "pancreas", "ovary")
  if(info=="pancancer") mot_tissues <- c("colon", "lung", "mammary_gland", "kidney", "ovary", "bladder", "uterus")
  
  # just considering mot organs and extracting the tissue specific genes 
  mot_tissues_specific_genes <- tissue_specific_genes[ tissue_specific_genes[,2] %in% mot_tissues, c(3,2)]
  compare_sd_genes <- lapply(my_sd_genes, compare_genes, mot_tissues_specific_genes)
  
  num_genes <-unlist(lapply(my_sd_genes, length))
  y <- unlist(lapply(compare_sd_genes, length))
  df <- data.frame(my_sd_values,y, num_genes)
  colnames(df) <- c("SD", "Tissue_Specific_Genes", "Features")

  output_file <- paste0(output_dir, "/", info, "SD_versus_TSG_Analysis_Plots.pdf" )
  # now producing the second plot

  g2<- ggplot(df, aes(x = SD, y = Tissue_Specific_Genes, label=Features)) +
    geom_point(size=5) +
    geom_text(aes(size=14, hjust=0.5, vjust=-1))+scale_x_continuous(breaks = round(seq(min(df$SD), max(df$SD), by = 0.1),2) ) +
    scale_y_continuous(breaks = round(seq(min(df$Tissue_Specific_Genes), max(df$Tissue_Specific_Genes), by = 20),1)) + theme_bw()+
    xlab(colnames(df)[1]) +
    ylab("Tissue specific genes (TIGER)") +
    ggtitle( paste0("Total number of tissue specific genes for ", info," is ", length( unique(mot_tissues_specific_genes[,1]))  ) ) + theme_Publication()

  # writing to the file  
  ggsave(g2, file=output_file)
  
  sd_file <- gsub(".pdf","_Comparision_with_Tiger.txt", output_file )
  write.table(df , sd_file, sep="\t")
  
  return(gene_sd_file)
  
}



# this function creates the plot for the SD versus tissue specific genes for every organ
create_tissue_sd_plots_by_organs <- function(gene_sd_file, output_dir, info, org){
  
  my_sd_data <- read.table(gene_sd_file, sep = "\t", header = TRUE)
  
  # sd values for the analysis
  my_sd_values<- c(0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5)
  
  get_sd_genes_org <- function(sd_value, sd_data){
    
    # selecting the sd and comparing with the tissue specific genes
    ind_sd <- which(sd_data[,-1] >= sd_value)
    # get the genes in sd 0.8
    genes_sd <- as.character(sd_data[ind_sd,1] )
    
    return(genes_sd)
    # after sd selections
    
    
  }
  
  my_sd_genes <- lapply(as.list(my_sd_values), get_sd_genes_org, my_sd_data)
  
  # just considering mot organs and extracting the tissue specific genes 
  mot_tissues_specific_genes <- tissue_specific_genes[ tissue_specific_genes[,2] %in% org, c(3,2)]
  compare_sd_genes <- lapply(my_sd_genes, compare_genes, mot_tissues_specific_genes)
  
  num_genes <-unlist(lapply(my_sd_genes, length))
  y <- unlist(lapply(compare_sd_genes, length))
  df <- data.frame(my_sd_values,y, num_genes)
  colnames(df) <- c("SD", "Tissue_Specific_Genes", "Features")
  
  org.1 <- paste(org, collapse = "_")
  output_file <- paste0(output_dir, "/", info, "_",org.1,"_SD_versus_TSG_Analysis_Plots.pdf" )
  # now producing the second plot
  
  
  g2<- ggplot(df, aes(x = SD, y = Tissue_Specific_Genes, label=Features)) +
    geom_point(size=5) +
    geom_text(aes(size=14, hjust=0.5, vjust=-1))+scale_x_continuous(breaks = round(seq(min(df$SD), max(df$SD), by = 0.1),2) ) +
    scale_y_continuous(breaks = round(seq(min(df$Tissue_Specific_Genes), max(df$Tissue_Specific_Genes), by = 10),1)) + theme_bw()+
    xlab(colnames(df)[1]) +
    ylab("Tissue specific genes (TIGER)") +
    ggtitle( paste0("Total number of tissue specific genes for ", info," ",org.1," is ", length( unique(mot_tissues_specific_genes[,1]))  ) ) + theme_Publication()
  
  # writing to the file  
  ggsave(g2, file=output_file)
  
  sd_file <- gsub(".pdf","_Comparision_with_Tiger.txt", output_file )
  write.table(df , sd_file, sep="\t")
  

  
}

