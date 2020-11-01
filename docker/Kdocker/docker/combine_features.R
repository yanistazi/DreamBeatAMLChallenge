combine_features <- function(feature_types,path="../data/features/") {
    
    combined_list<-list()
    for (i in feature_types){
        
        if (i=="RNA"){
            
            input <- read.table(paste(path,"prepared_data_rna_reduced.csv",sep=""),sep="\t",header=TRUE,row.names=1)
            
        }else if (i=="RNA_var"){
            
            input <- read.table(paste(path,"prepared_data_rna_vargenes.csv",sep=""),sep="\t",header=TRUE,row.names=1)
            
        }else if (i=="RNA_var_zscore"){
            
            input <- read.table(paste(path,"prepared_data_rna_vargenes_zscore.csv",sep=""),sep="\t",header=TRUE,row.names=1)

        }else if (i=="RNA_cor"){
            
            input <- read.table(paste(path,"prepared_data_rna_corgenes.csv",sep=""),sep="\t",header=TRUE,row.names=1)
            
        }else if (i=="RNA_cor_enriched"){
            
            input <- read.table(paste(path,"prepared_data_rna_corgenes_enriched.csv",sep=""),sep="\t",header=TRUE,row.names=1)
            
        }else if (i=="RNA_cor_enriched_zscore"){
            
            input <- read.table(paste(path,"prepared_data_rna_corgenes_enriched_zscore.csv",sep=""),sep="\t",header=TRUE,row.names=1)            

        }else if (i=="RNA_cor_zscore"){
            
            input <- read.table(paste(path,"prepared_data_rna_corgenes_zscore.csv",sep=""),sep="\t",header=TRUE,row.names=1)
            
        }else if (i=="RNA_var_cor"){
            
            input <- read.table(paste(path,"prepared_data_rna_varcorgenes.csv",sep=""),sep="\t",header=TRUE,row.names=1)

        }else if (i=="DNA"){
            
            input <- read.table(paste(path,"prepared_data_dna.csv",sep=""),sep="\t",header=TRUE,row.names=1)
  
        }else if (i=="clinical_categorical"){
            
            input <- read.table(paste(path,"prepared_data_clin_cat.csv",sep=""),sep="\t",header=TRUE,row.names=1)

        }else if (i=="clinical_categorical_reduced"){
            
            input <- read.table(paste(path,"prepared_data_clin_cat_reduced.csv",sep=""),sep="\t",header=TRUE,row.names=1)

        }else if (i=="clinical_numerical"){
            
            input <- read.table(paste(path,"prepared_data_clin_num.csv",sep=""),sep="\t",header=TRUE,row.names=1)

        }else if (i=="clinical_numerical_reduced"){
            
            input <- read.table(paste(path,"prepared_data_clin_num_reduced.csv",sep=""),sep="\t",header=TRUE,row.names=1)
    
        }else{
            
            print("The specified feature types are not supported\n")
            print('Allowed types: "RNA","DNA","clinical_numerical","clinical_numerical_reduced","clinical_categorical","clinical_categorical_reduced')
            return()
        }
        print(dim(input))
        combined_list[[i]]<-input
    }

    if  (length(feature_types)>1){
        
        #check that order of patients is the same across datasets
        ll <- combn(names(combined_list),2,simplify=FALSE)
        out <- lapply(ll, function(x) all.equal(rownames(x[[1]]), rownames(x[[2]])))
        if (all(unlist(out))){
        
            combined_df<-do.call(cbind,combined_list)
        
        }else{
        
            print("Patient order differs between feature types")
        }
    }else{
        
        combined_df<-combined_list[[i]]
    }
                  
    return(combined_df)
}


