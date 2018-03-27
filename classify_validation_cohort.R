
# performing the correlation assigners
classify_validation_cohort <- function(ucsc_pancreatic, my_pam_centroids){
  
  my_pam_centroids_1 <- apply( my_pam_centroids[, -1], 2, as.numeric)
  rownames(my_pam_centroids_1) <- my_pam_centroids[,1]
  
  ucscPDAdata <- ucsc_pancreatic[ ucsc_pancreatic[,1] %in% rownames(my_pam_centroids_1), ]
  ucscPDAdata_1 <- apply( ucscPDAdata[,-1], 2, as.numeric )
  rownames(ucscPDAdata_1) <- ucscPDAdata[,1]
  
  genes <- intersect(rownames(ucscPDAdata_1), rownames(my_pam_centroids_1))
  
  # only considering the intersecting genes
  centroids <- my_pam_centroids_1[genes, ]
  Exp <- ucscPDAdata_1[genes, ]
  obs_correlation <- as.data.frame(cor(Exp[genes,], centroids[genes,], use = "pair", method="pearson"))
  obs_correlation$lables <-apply( obs_correlation, 1, function(x) { ind=which.max(x); names(x)[ ind ]   }  )
  
  
  obs_correlation$maxCor <- apply( obs_correlation[,c(1:ncol(centroids))], 1, max)
  
  # get the second highest correlation
  
  obs_correlation$second.high <- apply( obs_correlation[, c(1:ncol(centroids))], 1, function(x){ ind=which.max(x); new_dat <- x[-ind]; ind_2 <- which.max(new_dat); return( new_dat[ind_2] )})
  obs_correlation$delta <- obs_correlation$maxCor - obs_correlation$second.high
  
  # now only selecting those which satifies the required condition
  obs_correlation$final_sub <- "ND"
  
  ind2sel <- which( obs_correlation$maxCor > 0.15 & obs_correlation$delta > 0.06)
  obs_correlation$final_sub[ind2sel] <- obs_correlation$lables[ind2sel]
  
  return(obs_correlation)
  
}

