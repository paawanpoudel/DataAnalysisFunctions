# this function corrects for the batch. it assigns the batch based on the sample information file

run_combat <-
function(df, sample_f){
    # loading the packages
    suppressPackageStartupMessages(require(sva))
    # renaming the colnames from the sample info file
    
    colnames(sample_f) <- c("samples", "batch")
    # assigning the batch for the sample info file
    
    # making the gene expression matrix to numberic (just in case)
    df_1 <- apply(df[,-1], 2, as.numeric)

    m1 <- match(sample_f[,1], colnames(df_1))
    w1 <- which(!is.na(m1))
    
    # creating the batch variable
    batch <- sample_f$batch[w1]
    
    # creating the model matrix by fitting an intercept. No adjustment variable has been included in creating the model matrix. Please note that the adjustment variable will be used later by the combat function
    mod <- model.matrix(~1, data=sample_f[w1,])
    
    
    # now running the combat with tumors/folder as a batch
    data_b_c_1 <- ComBat(dat=df_1[ , m1[w1] ], batch=batch, mod=mod)
    
    # calculating the median
    #rowMed <-  apply(data_b_c_1,1,median)
    
    # median centering the genes
    #df_f <- data_b_c_1 # - rowMed
    
    # creating the batch corrected data matrix
    dbc <- data.frame(df[,1], data_b_c_1 )
    colnames(dbc)[1] <- "Genes"
    
    return(dbc)
        
}
