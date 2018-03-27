produce_heatmap <-
function(pval, plot_file){
  
  require(RColorBrewer)
  require(reshape)
    
  df1<- pval
  info <- "Enrichment analysis using hypergeometric test"
  colnames(df1)[1:3] <- c("X1","X2","pvalue")
  
  
 #df1 <- df1[ order(df1$pvalue, decreasing = FALSE), ]


  #pp <- ggplot(df1, aes( x=reorder(X1, pvalue), y=X2, fill=pvalue)) +
  pp <- ggplot(df1, aes( x=X1, y=X2, fill=pvalue)) +
    geom_raster()+geom_tile(aes(fill = pvalue)) + 
    scale_fill_gradientn(colours=brewer.pal(7, 'RdYlBu'),na.value = "transparent",
                         breaks=seq(0,1, 0.1),labels=seq(0,1, 0.1),
                         limits=c(0,1))+theme(
                           axis.text.x  = element_text(colour="black",angle=90, vjust=0.5, size=14),
                           axis.text.y  = element_text(colour="black",angle=0, vjust=0.5, size=14),
                           panel.background = element_rect(fill = "white"))+ labs(title=info) + scale_colour_Publication()+ theme_Publication()
  
  ggsave(filename = plot_file, plot = pp)
  
  pval_file <- gsub(".pdf", ".txt", plot_file)
  write.table(df1, pval_file, sep="\t")
  
}
