combine_features <- function(feature_types,path="../data/features_SC2/") {
   # Combine different types of features in a single data frame
   #
   # INPUT : we also have auc now!
   # feature_types : a vector of strings (feasible lengths 1-4). The elements of the vector specify the feature types to be combined. Allowed strings: "RNA","DNA","clinical_numerical","clinical_categorical"  
   #                 E.g. feature_types <- c("RNA,DNA")
   #                 E.g. feature_types <- c("RNA,clinical_numerical")
   #                 E.g. feature_types <- c("clinical_numerical,clinical_categorical")
   #                 E.g. feature_types <- c("RNA") 
   # OUTPUT
   # a data.frame in which the individual features are merged 


   combined_list<-list()
   for (i in feature_types){

      if (i=="RNA"){ # RNA-PCA: 

	 input <- read.table(paste(path,"prepared_data_rna_reduced.csv",sep=""))

      }else if (i=="RNA_var"){ # RNA most variant genes

	 input <- read.table(paste(path,"prepared_data_rna_vargenes.csv",sep=""))

      }else if (i=="RNA_var10"){ # RNA most variant genes

	 input <- read.table(paste(path,"prepared_data_rna_vargenes10.csv",sep=""))

      }else if (i=="RNA_var20"){ # RNA most variant genes

	 input <- read.table(paste(path,"prepared_data_rna_vargenes20.csv",sep=""))

      }else if (i=="RNA_var30"){ # RNA most variant genes

	 input <- read.table(paste(path,"prepared_data_rna_vargenes30.csv",sep=""))

      }else if (i=="RNA_hr"){ # top significant genes in univariate Hazard Ratio analysis

	 input <- read.table(paste(path,"prepared_data_rna_hrgenes.csv",sep=""))

      }else if (i=="RNA_hr10"){ # top significant genes in univariate Hazard Ratio analysis

	 input <- read.table(paste(path,"prepared_data_rna_hrgenes10.csv",sep=""))

      }else if (i=="RNA_hr20"){ # top significant genes in univariate Hazard Ratio analysis

	 input <- read.table(paste(path,"prepared_data_rna_hrgenes20.csv",sep=""))

      }else if (i=="RNA_hr30"){ # top significant genes in univariate Hazard Ratio analysis

	 input <- read.table(paste(path,"prepared_data_rna_hrgenes30.csv",sep=""))

      }else if (i=="RNA_hr_var"){ # top significant genes in univariate Hazard Ratio analysis in RNA_var

	 input <- read.table(paste(path,"prepared_data_rna_hrvargenes.csv",sep=""))

      }else if (i=="mol"){ # clinical categorical and DNA mutation consolidated together as "molecular" (+sex + therapy-related)
	 input <- read.table(paste(path,"prepared_data_clin_cat_mut_clean.csv",sep=""))

      }else if (i=="clinical_numerical"){ # clinical numarical [Age Blast WBC] with imputation for Blast WBC

	 input <- read.table(paste(path,"prepared_data_clin_num_imputed.csv",sep=""))

      }else if (i=="AUC"){ # drug AUC with imputation

	 input <- read.table(paste(path,"prepared_data_auc_imputed.csv",sep=""))

      }else if (i=="combined_features"){ # features common with MSK data

	 input <- read.table(paste(path,"prepared_data_challenge_combined.csv",sep=""))
      
      }else if (i=="combined_features_reduced"){ # features common with MSK data

	 input <- read.table(paste(path,"prepared_data_challenge_combined_reduced.csv",sep=""))

      }else{

	 print("The specified feature types are not supported\n")
	 print('Allowed types: "RNA","DNA","clinical_numerical","clinical_numerical_reduced","clinical_categorical","clinical_categorical_reduced',"AUC","combined_features")
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


