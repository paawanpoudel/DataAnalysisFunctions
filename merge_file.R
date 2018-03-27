merge_file <-
function(file_list){
    
    
    final_df=data.frame()
    
     # Reading the individual files
    genes_1=read.delim2(file_list[1], sep="\t", row.names=NULL, header=TRUE)
    
     # checking if there are any genes with "NA"
    index=which(!is.na(genes_1[,1]))
    final_df=data.frame(genes_1[index,])
    
    colnames(final_df)[1]="Genes"
    
    if(length(file_list)>=2){
      
      for (f in 2:length(file_list))
        {
        
        data_1=read.delim2(file_list[f], sep="\t", row.names=NULL, header=TRUE)
        colnames(data_1)[1]="Genes"
        
        final_df=merge(final_df, data_1, by.x="Genes", by.y="Genes")
        
      }
      
      return(final_df)
      
    }else{

      return(final_df)
    
    }
}
