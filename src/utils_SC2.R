# Some helper functions
library(stringr)

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

      }else if (i=="RNA_hr"){ # top significant genes in univariate Hazard Ratio analysis

	 input <- read.table(paste(path,"prepared_data_rna_hrgenes.csv",sep=""))

      }else if (i=="RNA_hr10"){ # top significant genes in univariate Hazard Ratio analysis

	 input <- read.table(paste(path,"prepared_data_rna_hrgenes10.csv",sep=""))

      }else if (i=="RNA_hr20"){ # top significant genes in univariate Hazard Ratio analysis

	 input <- read.table(paste(path,"prepared_data_rna_hrgenes20.csv",sep=""))

      }else if (i=="RNA_hr30"){ # top significant genes in univariate Hazard Ratio analysis

	 input <- read.table(paste(path,"prepared_data_rna_hrgenes30.csv",sep=""))

      }else if (i=="mol"){ # clinical categorical and DNA mutation consolidated together as "molecular" (+sex + therapy-related)
	 input <- read.table(paste(path,"prepared_data_clin_cat_mut_clean.csv",sep=""))

      }else if (i=="clinical_numerical"){ # clinical numarical [Age Blast WBC] with imputation for Blast WBC

	 input <- read.table(paste(path,"prepared_data_clin_num_imputed.csv",sep=""))

      }else if (i=="AUC"){ # drug AUC with imputation

	 input <- read.table(paste(path,"prepared_data_auc_imputed.csv",sep=""))

      }else if (i=="combined_features"){ # features on data challenge common with MSK data

	 input <- read.table(paste(path,"prepared_data_challenge_combined.csv",sep=""))
          
      }else if (i=="combined_features_reduced"){ # features on data challenge  common with MSK data reduced

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

launch_prediction_SC2 <- function(features_model,predictors,str_predictors,nfolds=5 , nrepeats=10 , mc.cores=4,l_alpha=0,l_ntree=50,l_nodesize=20,
				  l_maxstepno=500 , l_K=10 , l_type="verweij" , l_penalty=100 , l_max.iter = 500,l_tol=0.01,
				  l_use_MSK=0, ...){

   # INPUT launch_prediction : (see full example in launch_prediction_example in folder prediction)

   # Mandatory Inputs:
   #  features_model: list of features to test
   #  drug_responses : list of drugs to test (default is "all") and test all drugs . User can specify either drug names or its index from the response dataframe
   #  predictors : list of predictors that the user wants to use 
   #  str_predictors : same as before but as character list (see example )
   # Optional Inputs
   # Any list of hyperparameters (l_alpha for elasticnet , l_ntree, l_nodesize for RF , l_max.depth,l_eta,l_nthread,l_nrounds for XGBoost and l_k for kNN)
   #If you launch a predictor without its hyperparameter, I set up some default parameters . You can specify a list for each of the hyperparameters and the algo will do a for loop         #overall of them
   # Other parameters as number of cores nfolds and nrepeats are also optional.
   # OUTPUT
   # One dataframe in folder predictions/results_SC1 by feature model with : nrows= number of predictors for each hyperparameters ( example : elasticnet with 10 values of alpha and RF with 5 values of trees will give 15 rows)
   #ncols= 7(numbers of metrics from eval CV) * number of drug_responses specified . If drug_responses = "all" : NCOLS = 7*122=854

   command_init <- "bsub -R 'rusage[mem=4]' -W 05:59 -n 4 -J 'model' -o "

   for (feature_model in features_model){

      combination <- paste(feature_model,collapse="with")

      if (!dir.exists(paste("./../prediction/results_SC2/",combination,sep=""))){
	 dir.create(paste("./../prediction/results_SC2/",combination,sep=""))
      }

      if (!dir.exists(paste("./../prediction/results_SC2/",combination,"/stdoutput",sep=""))){
	 dir.create(paste("./../prediction/results_SC2/",combination,"/stdoutput",sep=""))
      }

      outpath <- paste("./../prediction/results_SC2/",combination,sep="")
      str <- 1
      for (predictor in predictors){

	 if(identical(predictorGLM,predictor)){
	    for (alpha in l_alpha){
	       for (use_MSK in l_use_MSK){
		  filename <- paste(str_predictors[str],alpha,use_MSK,sep="_")
		  command <- paste(command_init, paste("'../prediction/results_SC2/",combination,"/stdoutput/",filename,".txt' ",sep=""),
				   "'Rscript evaluate_SC2_intermediate_parallel.R --predictor ",str_predictors[str]," --feature_model ",combination," --use_MSK ",use_MSK," --alpha ",alpha," --nrepeats ",nrepeats," --nfolds ",
				   nfolds, " --mc.cores ", mc.cores, " --outpath ",outpath, " --filename ",filename,"'",sep="") 
		  system(command)
		  print(command)
	       }
	    }
	 }

	 if(identical(predictorRF,predictor)){
	    for (ntree in l_ntree){
	       for (nodesize in l_nodesize){
		  for (use_MSK in l_use_MSK){
		     filename <- paste(str_predictors[str],ntree,nodesize,use_MSK,sep="_")
		     command <- paste(command_init, paste("'../prediction/results_SC2/",combination,"/stdoutput/",filename,".txt' ",sep=""),
				      "'Rscript evaluate_SC2_intermediate_parallel.R --predictor ",str_predictors[str]," --feature_model ", combination," --use_MSK ",use_MSK,
				      " --ntree ",ntree, " --nodesize ",nodesize," --nrepeats ",nrepeats," --nfolds ",
				      nfolds, " --mc.cores ", mc.cores, " --outpath ",outpath, " --filename ",filename,"'",sep="") 
		     system(command)
		     print(command)
		  }
	       }
	    }
	 }

	 if(identical(predictorBoost,predictor)){
	    for (maxstepno in l_maxstepno){
	       for (K in l_K){
		  for (type in l_type){
		     for (penalty in l_penalty){
			for (use_MSK in l_use_MSK){
			   filename <- paste(str_predictors[str],maxstepno,K,type,penalty,use_MSK,sep="_")
			   command <- paste(command_init, paste("'../prediction/results_SC2/",combination,"/stdoutput/",filename,".txt' ",sep=""),
					    "'Rscript evaluate_SC2_intermediate_parallel.R --predictor ",str_predictors[str]," --feature_model ", combination," --use_MSK ",use_MSK," --maxstepno ",maxstepno, " --K ",K," --type ",type, " --penalty ",penalty,
					    " --nrepeats ",nrepeats," --nfolds ", nfolds," --mc.cores ", mc.cores, " --outpath ",
					    outpath, " --filename ",filename,"'",sep="") 
			   system(command)
			   print(command)
			}
		     }
		  }
	       }
	    }
	 }

	 if(identical(predictorRFX,predictor)){
	    for (iter in l_max.iter){
	       for(tol in l_tol){
		  for (use_MSK in l_use_MSK){
		     filename <- paste(str_predictors[str],iter,tol,use_MSK,sep="_") 
		     command <- paste(command_init, paste("'../prediction/results_SC2/",combination,"/stdoutput/",filename,".txt' ",sep=""),
				      "'Rscript evaluate_SC2_intermediate_parallel.R --predictor ",str_predictors[str]," --feature_model ", combination," --use_MSK ",use_MSK,
				      " --iter ",iter," --tol ",tol, " --nrepeats ", nrepeats," --nfolds ", nfolds, " --mc.cores ", mc.cores,
				      " --outpath ",outpath, " --filename ",filename,"'",sep="")
		     system(command)
		     print(command)
		  }
	       }
	    }
	 }   

	 str <- str + 1
      }
   }
}




get_models <- function(path = "../prediction/results_SC2",dataset="combined_features"){

   combinations <- list.dirs(path = path, full.names = FALSE, recursive = FALSE)
   combinations <- combinations[which(combinations!=".ipynb_checkpoints")]
   combinations <- combinations[which(combinations %in% dataset)]
   list_files <- c()
   for (comb in combinations){
      comb_path <- paste(path,"/",comb,sep="")
      files_tosearch <- list.files(path = comb_path, pattern = "\\.csv$")
      list_files <- c(list_files,paste(comb_path,sep="/",files_tosearch))
   }
   #list_files <- list_files[str_detect(list_files,pattern=“.csv”)]
   df <- data.frame(mean=numeric(0),sd=numeric(0),min=numeric(0),max=numeric(0))
   list_files <- list_files[!str_detect(list_files,pattern="Boost")]   ### to remove when things work !!!!!!!
   for (l in list_files){
      to_add <- read.table(l)
      rownames(to_add) <- str_remove(l,path)
      df <- rbind(df,to_add)
   }
   return (df)
}              




save_trained_best_model_SC2 <- function(df_results,data_response,path.features="../data/features_SC2/"){

    save_rds <- function(final_model,name){
      saveRDS(final_model, paste(name,".rds",sep=""))
    }
    
    prepare_data_with_MSK <- function(features,response,use_MSK){
        if("NPM1" %in% colnames(features)){
          MSK_data <- read.table(paste0(path.features,"similarity_ordered_prepared_data_MSK.csv"))
            print("a")
        }else{
            print("reduced")
          MSK_data <- read.table(paste0(path.features,"similarity_ordered_prepared_data_MSK_reduced.csv"))
        }
        MSK_train <- MSK_data[,colnames(features)]

        MSK_response <- MSK_data[,names(MSK_data) %in%c("OS","OS_status")]
        MSK_response$time <- MSK_response$OS
        MSK_response$status <- MSK_response$OS_status
        MSK_response$OS <- NULL
        MSK_response$OS_status <- NULL
        num_rows <- round(nrow(MSK_data)*use_MSK/100)  # num of row to use based on use_MSK percentage 
        features <- rbind(features,MSK_train[1:num_rows,])
        response <- rbind(response,MSK_response[1:num_rows,c("time","status")])
        return(list(features=features,response=response))
    }
        

    

   best <- rownames(df_results[order(df_results$mean,decreasing=T),,drop=F])[1]
   best <- str_replace(best,".csv","")
   features <- as.list(strsplit(as.character(best), "/")[[1]])[2]
   features <- as.list(strsplit(as.character(features), "with")[[1]])
   features <- combine_features(features,path=path.features)
   response <- data_response
   response <- na.omit(response)
   response <- response[response$OS>0,]

   non_na <- intersect(rownames(na.omit(features)),rownames(na.omit(response)))
   features <- features[non_na,]
   response <- response[non_na,]
   response$time <- response$OS 
   response$status <- response$OS_status
   response$OS <- NULL
   response$OS_status <- NULL



   model_params <- as.list(strsplit(as.character(best), "/")[[1]])[3]
   model <- as.list(strsplit(as.character(model_params), "_")[[1]])[1]

   name <- substring(best,2) # to name the model
   name <- str_replace(name,"/predictor","withpredictor")
   if (model=="predictorRF"){
      nntree <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[2])
      nodesize <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[3])
      use_MSK <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[4])
      if(use_MSK>0){
          features <- prepare_data_with_MSK(features,response,use_MSK)$features
          response <- prepare_data_with_MSK(features,response,use_MSK)$response
      }
      final_model <- rfsrc(Surv(time, status) ~ .,
               data=data.frame(features,response), ntree=nntree, nodesize=nodesize, importance="none")
      save_rds(final_model,name)
   } else if(model=="predictorRFX"){
      max.iter <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[2])
      tol <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[3])
      use_MSK <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[4])
      if(use_MSK>0){
          features <- prepare_data_with_MSK(features,response,use_MSK)$features
          response <- prepare_data_with_MSK(features,response,use_MSK)$response
      }
      final_model <- CoxRFX(data.frame(features), Surv(time=response[,1],event =response[,2]) ,
            max.iter =max.iter,tol=tol,groups=rep(1, ncol(features)))
      save_rds(final_model,name)
   } else if(model=="predictorBoost"){
      maxstepno <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[2])
      K <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[3])
      type <- as.character(as.list(strsplit(as.character(model_params), "_")[[1]])[4])
      penalty <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[5])
      use_MSK <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[6])
      if(use_MSK>0){
          features <- prepare_data_with_MSK(features,response,use_MSK)$features
          response <- prepare_data_with_MSK(features,response,use_MSK)$response
      }       
      cv.res<-cv.CoxBoost(time=as.numeric(unlist(response[,1])),
              status=as.numeric(unlist(response[,2])),
              x=data.matrix(features),maxstepno=maxstepno,K=K, type=type , penalty=penalty )
      final_model <- CoxBoost(time=as.numeric(unlist(response[,1])),
                  status=as.numeric(unlist(response[,2])),
                  x=data.matrix(features),stepno=cv.res$optimal.step , penalty=penalty)
      save_rds(final_model,name)
   } else if(model=="predictorGLM"){
      alpha <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[2])
      use_MSK <- as.numeric(as.list(strsplit(as.character(model_params), "_")[[1]])[3])
      if(use_MSK>0){
          features <- prepare_data_with_MSK(features,response,use_MSK)$features
          response <- prepare_data_with_MSK(features,response,use_MSK)$response
      }
      final_model <- cv.glmnet(data.matrix(features), data.matrix(response),
                   family="cox", alpha=alpha, grouped=TRUE)

      save_rds(final_model,name)
   }

}



predict_best_model_SC2 <- function(rds_path,path="../data/features_SC2/"){

   model <- readRDS(rds_path)

   len <- length(as.list(strsplit(as.character(rds_path), "with")[[1]]))
   features <- as.list(strsplit(as.character(rds_path), "with")[[1]])[1:(len-1)]
   predictor <- str_replace(as.list(strsplit(as.character(rds_path), "with")[[1]])[len],".rds","")

   features <- na.omit(combine_features(features,path=path))    #### replace path with input/ submission
   if(grepl("predictorRF_",predictor)){
      return(data.frame(survival=predict(model,features)$predicted,row.names=rownames(features)))
   } else if (grepl("predictorRFX_",predictor)){
      print(head(t(-predict(model,features))))
      print(head(t(predict(model,features))))
      return(data.frame(survival= - predict(model,features)))
   } else if (grepl("predictorBoost",predictor)){
      return(data.frame(survival=t( - predict(model,features))))
   } else if (grepl("predictorGLM",predictor)){
      return(data.frame(survival= - predict(model,newx=data.matrix(features), s="lambda.min", type="response")[,1]))
   } 
   # I put sign minus for Boost RFX and GLM because the outputs correspond to beta*X 
}

