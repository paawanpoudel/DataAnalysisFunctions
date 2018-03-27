create_subtype_prop_plot <-
function(df, prop_pdf_file){
    

    colnames(df) <- c("samples", "organs","membership")
    
    df$membership <- gsub("\\.score|\\-score", "",df$membership)
    df$membership <- gsub("\\.", "\\-",df$membership)
    
    
    # getting the color information for each subtypes
    my_col <- create_subtype_col()
    
        
  	data2plot_fig1 <- df[ order( df$organs, df$membership), ]
  	# creating the order of the plot so that mixed subtype in on the top
  	
  	
  	
  
  	# calculating the proportion of mixed subtypes and then sorting it in ascending order for plotting
  	xx <- table(data2plot_fig1$organs, data2plot_fig1$membership)
  	
  	
  	ind_sel_nd <- grep( "ND", df$membership)
  	
  	if(length(ind_sel_nd)>=1){
  		
  		my_order <- c("ND",  "Stem-like", "SP","Inflammatory","Basal", "TA", "Classical")
  		my_proportions <- prop.table(xx, margin = 1)[,"ND"]
  	}else{
  		
  		my_order <- unique(df$membership)
  		#my_order <- c("Inflammatory","Basal", "TA", "Stem-like", "SP","Classical")
  		#my_proportions <- prop.table(xx, margin = 1)[,"Basal"]
  	}
  	
  	
  	data2plot_fig1$order <-factor(data2plot_fig1$membership, levels=my_order )
  	  	
  	my_organ_order <- names(sort(my_proportions))
  	data2plot_fig1$my_organ_order <-factor(data2plot_fig1$organs, levels=my_organ_order)
  
  	# matching the color information with the 6 subtypes and assigning colours accordingly
    
  	m1 <- match( levels(data2plot_fig1$order),my_col[,1] )
    w1 <- which(!is.na(m1))
    
    sub.col <- my_col[m1[w1],2]
    
    if(length(w1)==0){
        
        n <- length( unique(df$membership) )
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        sub.col <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:length(unique(df$membership))]
        
    }

  
  
  
  p <- ggplot( data2plot_fig1, aes(x= my_organ_order, fill=order))+ geom_bar(position="fill", aes(fill = order))
  p <- p + scale_fill_manual(values = sub.col) + theme_Publication()+ labs(title="Proportions of samples in each subtype",x="Organs", y = "Proportions")
  
   
       # creating the bar plot
#    p <- ggplot(df, aes(x= organs, fill=membership))+ geom_bar(position="fill",aes(fill = membership))
#    p <- p + scale_fill_manual(values = sub.col)+ scale_colour_Publication()+ theme_Publication()+ labs(title="Proportions of samples in each subtype",x="Organs", y = "Proportions")
    
 
    ggsave(p, file=prop_pdf_file)
    
}
