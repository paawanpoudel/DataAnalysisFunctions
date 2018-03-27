# this function creates the sample information file, which is used to color the organs and perform combat
create_sample_info <- function(file, output_dir){


        trainning_geo_tmp_1 <- dirname(file)
        trainning_geo_tmp <- tail(unlist( strsplit(trainning_geo_tmp_1, "/")), 1 )
        
        #locating of the sample info file for the trainning datasets. Sample info files are very important to create the plots
        sample_info_file=paste0(output_dir,"/",trainning_geo_tmp,"_sample_info.txt")
        
        if(file.exists(sample_info_file)){
        
           cat("Using the sample info file")
        
        }else{
		 
          cel.files <- read.table(file, nrows=1, sep="\t")
          
          org=rep(paste0(trainning_geo_tmp, length(cel.files)-1))
          
          s_info=data.frame(t(cel.files)[-1], org)
          
          colnames(s_info)=c("geo_accession", "label")
          
          write.table(s_info, sample_info_file, sep="\t", quote=FALSE, row.names=FALSE)

        }


}



