create_subtype_col <-
function(){
    
    colours <- c( "#fdbb84",  "#e34a33", "#67000d", "#08519c","#08519c", "#addd8e","#addd8e","#addd8e", "#df65b0","#756bb1", "#addd8e", "#addd8e","#e34a33",  "#beaed4","grey" )
    subtype <- c("Basal","Classical", "Inflammatory", "Stem-like", "Stem.like", "Stem-PPAR", "SP","Stem.PPAR", "TA", "TA/SP", "Goblet-like", "Goblet.like","Enterocyte", "CSR","ND")
    
    sub.col <- data.frame(subtype, colours)
    sub.col <- sub.col[order(sub.col[,1]),]
    
    return(sub.col)
}
