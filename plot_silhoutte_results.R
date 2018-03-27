


plot_silhoutte_results <- function(silhoutte_data, plot_file){
  
  df <- silhoutte_data[,c(1,2,5)]
  colnames(df) <- c("samples", "membership", "sil_width")
  
  df <- df[order(df$membership, df$sil_width), ]
  df$name <- factor(rownames(df), levels = rownames(df))
  
  
  df$membership <- as.factor(df$membership)
  df$sil_width <- as.numeric(df$sil_width)
  
  mapping <- aes_string(x = "name", y = "sil_width", fill="membership")
  ave <- tapply(df$sil_width, df$membership, mean)
  n <- tapply(df$membership, df$membership, length)
  sil.sum <- data.frame(cluster = names(ave), size = n, ave.sil.width = round(ave, 2))
  
  my_col <- c( "Basal"="#fdbb84", "Classical"="#e34a33","SP"="#addd8e", "Inflammatory"="#67000d", "Stem-like"="#08519c",  "TA"="#df65b0")
  
  p <- ggplot(df, mapping) + geom_bar(stat = "identity") +  labs(y = "Silhouette width Si", x = "", title = paste0("Clusters silhouette plot ", "\n Average silhouette width: ", round(mean(df$sil_width),  2))) #+ ggplot2::ylim(c(NA, 1))
  q <-   p + scale_fill_manual(values = my_col ) + scale_colour_Publication() + theme_Publication() +theme( axis.text.x=element_blank() )
  
  ggsave(q, file=plot_file)
  return(sil.sum)
  
}