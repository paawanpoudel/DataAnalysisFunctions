create_subtype_prop_plot_hypermutated <-
function(df, plot_file){
    

    colnames(df) <- c("samples", "organs","membership")
    
    df$membership <- gsub("\\.score", "",df$membership)
    df$membership <- gsub("\\.", "\\-",df$membership)
    
    
    # getting the color information for each subtypes
    my_col <- create_subtype_col()
    
    # matching the color information with the 6 subtypes and assigning colours accordingly
    m1 <- match( my_col[,1],df$membership)
    w1 <- which(!is.na(m1))
    
    sub.col <- my_col[w1,2]
    
    if(length(w1)==0){
        
        n <- length( unique(df$membership) )
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        sub.col <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:length(unique(df$membership))]
        
    }
    
       # creating the bar plot
    p <- ggplot(df, aes(x= organs, fill=membership))+ geom_bar(position="fill",aes(fill = membership))
    p <- p + scale_fill_manual(values = sub.col)+ scale_colour_Publication()+ theme_Publication()+ labs(title="Proportions of samples in each subtype",x="Organs", y = "Proportions")
    
    
    ggsave(p, file=plot_file)
    
}
