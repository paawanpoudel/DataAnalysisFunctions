# this function creates the PCA plot for the samples and then it colours the samples by organs

makePcaPlot <-
function(x, sample_info, title = "", plot_file) {
    
    # loading the required packages
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(reshape))
    suppressPackageStartupMessages(require(RColorBrewer))
    
    #arranging the data
    sample_info <- sample_info[order(sample_info[,2]), ]
    
    # loading the color data
    col_data <- create_organ_colour()
    
    # resetting the coloumn names for the sample info file
    colnames(sample_info) <- c("sample", "organ")
    # matching the color data with the cancer types
    #sample_info$organ <- gsub("\\d", "",sample_info$organ )
    #sample_info$organ <- gsub("Rectum", "Colon",sample_info$organ, ignore.case = TRUE )
    # removing the extra spaces from the color data
    col_data[,2] <- gsub(" ","", col_data[,2])

    # removing the extra spaces from the organs
    sample_info$organ <- gsub(" ","", sample_info$organ)
    # get the unique organs from the sample info file
    org_info <- unique(sample_info$organ)
    
    # matching the organ data with the sample info file
    mc <- match( toupper(gsub(" ","",org_info) ) ,toupper(col_data[,1]) )
    wc <- which(!is.na(mc))
    
    # if no organs match then assign new colors to the organ variables using the RColorBrewer
    if(length(wc)==0){
        
        n <- length(org_info)
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        org.col <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:length(org_info)]
        
    }else{
        # if there is a matching organ information then assing the color for the create_organ_colour() function
        org.col <- col_data[mc[wc],2]
    }
    
    # creating the colors for the samples
    my_col <- factor(sample_info$organ, labels=org.col)
    
    # subsetting the data based on the common samples between the gene expression matrix and the sample info file
    data.1 <- subset(x, select=as.character(sample_info[,1]))
    # convering the data to the numberic matrix
    data.1 <- apply(data.1, 2, as.numeric)
    
    # get the rows where the variance is very small
    ind <- apply( data.1, 1, var)
    ind_sel <- which(ind ==0)
    # ignoring the rows with the 0 variance
    if(length(ind_sel) >=1){
    	
    	dx <- data.1[-which(ind ==0), ]
    }else{
    
    	dx <- data.1
    }
    
    
    # transposing the matrix
    mydata <- t(dx)
    # now the genes are in the columns and the samples are in the rows
    rownames(mydata) <- colnames(data.1)
    # running the PCA analysis. Please note that this scaling is actually scaling the genes not the sampels
    mydata.pca <- prcomp(mydata, scale.=TRUE)
    
    # get the scores for the PCA plots
    scores <- mydata.pca$x
    df <- data.frame(rownames(scores), scores[,c(1,2)], sample_info[,2])
    colnames(df) <- c("samples", "PCA1", "PCA2", "Groups")
    
    p <- ggplot(df, aes(PCA1, PCA2, colour = factor(Groups))) + geom_point(size = I(2)) +
    labs(x="PCA-1", y="PCA-2", colour="Groups")+
    scale_colour_manual(values = levels(my_col))
    
    p <- p + labs(title = "PCA plot") + theme_bw()
    
    score_file <- gsub(".pdf", "pca.scores.txt",plot_file)
    
    write.table(scores, score_file, sep="\t")
    
    ggsave(p, filename=plot_file)

}
