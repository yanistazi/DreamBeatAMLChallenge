source("combine_features.R")

library(caret)
library(randomForest)
library(xgboost)
library(glmnet)
library(e1071)
library(kmr, lib.loc="fordocker/")

print("KMR")

predict_from_best_model_kmr_switch <- function(drug_df, pathkernel = "/kernels/" , namekernel="KrnarestrictMean.txt", drug_switch_file="drug_switch_kmr_rf.txt",metric){

   drug_switch = read.table(drug_switch_file, header=F, sep="\t", stringsAsFactors=F)$V1 

   df_output <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("lab_id", "inhibitor", "auc"))))

   list.output <- list()

   for (i in 1:nrow(drug_df)){

      print(paste0("Drug",i))
      print(drug_df[i,2])
      
      if (!drug_df[i,2]%in%drug_switch) {

	 print("KMR prediction")

	 our_model <- readRDS(paste("final_models","/loaded_",drug_df[i,2],".rds",sep=""))
	 newK = read.table(paste0(pathkernel,namekernel),sep="\t",header=T,stringsAsFactors=F)
	 newK = newK[,our_model$non_na_index]

	 # make predictions on "new data" using the final model
	 kmr.object <- our_model$object
	 kmr.lambda <- our_model$lambda
	 print(paste0("KMR Lambda is:",kmr.lambda))

	 final_predictions <- predict(object=kmr.object, newx=as.matrix(newK), lambda=kmr.lambda)[,i]

      } else {

	 print("not using kmr")
	 our_model <- readRDS(paste("final_models/loaded_",drug_df[i,2],".rds",sep=""))
     print(class(our_model))
	 features <- c("RNA_cor_enriched")
	 celllines <- data.matrix(combine_features(features , path="/features/" ))
	 final_predictions <- predict(our_model, celllines)

      }

      if (length(unique(final_predictions))==1) {
	 print("Constant Predictions!")
      }
      if (any(is.na(final_predictions))) {
	 print("You have NA!")
      }

      df.o = data.frame(lab_id=rownames(newK) , 
			inhibitor=rep(drug_df[i,1],length(final_predictions)) , 
			auc=as.vector(final_predictions))
      list.output[[i]] = df.o

      #df_output[(nrow(df_output)+1):(nrow(df_output)+length(final_predictions)),c("lab_id","inhibitor","auc")] <- c(unlist(as.list(rownames(celllines))),rep(rownames(data_matrix)[i],nrow(celllines)),as.vector(final_predictions))
   }

   df_output = do.call(rbind, list.output)
   print(dim(df_output))

   #matched_df <- read.table("/docker/drug_names_matched.tsv", sep="\t",header=T)
   #df_output$inhibitor <- as.character(sapply(as.character(df_output$inhibitor), function(y){matched_df[which(matched_df$modified==y),"given"]}))

   print(head(df_output,2))
   write.table(df_output,"/output/predictions.csv",sep=",",row.names=FALSE,quote=FALSE)  ### Output done on train data (just to see if it works)
   return()
}


