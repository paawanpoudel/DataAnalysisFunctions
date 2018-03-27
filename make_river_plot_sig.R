make_river_plot_sig <-
function(b, org, plot_file){
    
    # getting the color information for each subtypes
    my_col <- create_subtype_col()
    
    # matching the color information with the 6 subtypes and assigning colours accordingly
    m1 <- match(rownames(b), my_col[,1])
    w1 <- which(!is.na(m1))
    
    subclass.col <- my_col[m1[w1],2]
    
    colnames(b) <- paste0(org,"-", colnames(b))
    
    pp <- round(hypergeometric_test(b), digits = 3)
    xs <- melt(pp)
    colnames(xs) <- c("N1", "N2","adj_pval")
    
    ind_sel <- which(xs$adj_pval <=0.05)
    xs_sel <- xs[ind_sel,]
    
    nodes <- data.frame(ID=unique(c(colnames(b), rownames(b) ) ),
    x=c(rep(1, length(colnames(b))), rep(2, length(rownames(b))) ),
    col= c( rep(NA, length(colnames(b))), subclass.col ),
    labels=unique(c(colnames(b), rownames(b))), stringsAsFactors=FALSE)
    
    edges <- melt(b)
    colnames(edges) <- c("N1", "N2", "Value")
    
    # combing the significant results
    edges_f <- merge(edges, xs_sel, by = c("N1", "N2"))
    
    r <- makeRiver( nodes, edges_f[,-ncol(edges_f)] )
    pdf(plot_file)
    riverplot( r ,node_margin = 0.3, nodewidth = 0.8, plot_area = 0.75, srt=0,yscale = 0.6)
    dev.off()
    
}
