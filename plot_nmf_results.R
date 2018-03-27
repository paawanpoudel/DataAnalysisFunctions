# this function plots the results from NMF
plot_nmf_results <-
function(datadir,details, init=2, final=10){
	
    # loading the required packages
    suppressPackageStartupMessages(require(NMF))
    suppressPackageStartupMessages(require(ggplot2))
    
	#setting the data dir
	setwd(datadir)
    # printing the working directory
    cat(datadir)
    	#getting the RData object produced from NMF
	file <- list.files(path=datadir, pattern=".RData", full.names=TRUE)
    
    #loading the RData object
    load(file)
    
    res <- results_nmf
    #setting the value of init and finally, generally the objects are stored in n-1 fashion, so need to do -1
	file_con <- gsub("\\.RData", "_cophenetic_coefficient.pdf", file)
	
	file_con_2 <- gsub("\\.RData", "_cophenetic_coefficient.txt", file)
	write.table(res$measures, file_con_2, sep="\t")
	
	cat("The cophenetic correlation coefficient plot is present in the file- \n ",file_con, "\n")
	
    # saving the nmf plots
	ggsave(plot(res), filename=file_con)

    #resetting the init and final because they are name in n-1 way
	init <- init -1
	final <- final -1
	
    #iterating throught the first and the last object
	for (i in init:final) 
	{
		j=i+1
       
        #accessing the res$fit objects to find the sample classification
		nmf_class=predict(res$fit[[i]], 'samples', prob=TRUE)
        
        #creating the sample file for the classification
		classification_nmf=data.frame(nmf_class[1])
        
		cran_gct_file=data.frame(rownames(classification_nmf), classification_nmf)
		rownames(cran_gct_file)=NULL
        
        #setting the coloumn names
		colnames(cran_gct_file)=c("Names", "membership")
		g_file=paste0(details,"_samples.k.", j, ".txt")
        
        #writing the sample classification to a file
		write.table(cran_gct_file, g_file, row.names=FALSE, sep="\t")
		
        #getting the metagenes information
        my_metagenes= predict(res$fit[[i]], "features", prob = TRUE, dmatrix = TRUE)
        
        #creating the data from of metagenes
        pf=data.frame(my_metagenes)
        m_file=paste0(details, "_metagenes.k.", j, ".txt")
        
        #writing the metagenes information to a file
        write.table(pf, m_file, sep="\t")
        
        }
}
