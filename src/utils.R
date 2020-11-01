# Some helper functions

compute_spearman <- function(responsePred , responseTest) {
	# Compute the Spearman Correlation between the prediction and the response
	#
	# INPUT
	# responsePred : qx1 matrix of predictions, q is the number of patients present in the test set  
	# responseTest : qx1 matrix of true values, q is the number of patients present in the test set  
	#
	# OUTPUT
	# a list
	# rho : correlation value
	# pvalue : significance of the correlation test
	
    
    responsePred <- as.matrix(responsePred)
    responseTest <- as.matrix(responseTest)
    
    #make sure that the order of patients in the two matrices is the same
    if (!isTRUE(all.equal(rownames(responsePred),rownames(responseTest)))){
        
        print("Order of patients is not the same")
        responsePred <- responsePred[rownames(responseTest),,drop=FALSE]
    }

    #corr <- cor.test(x=as.numeric(responsePred[,1]), y=as.numeric(responseTest[,1]), method = 'spearman')
    corr <- cor(x=as.numeric(responsePred[,1]), y=as.numeric(responseTest[,1]), method = 'spearman')

    return(list(rho=corr))
}

combine_features <- function(feature_types,path="../data/features/") {
	# Combine different types of features in a single data frame
	#
	# INPUT
	# feature_types : a vector of strings (feasible lengths 1-4). The elements of the vector specify the feature types to be combined. Allowed strings: "RNA","DNA","clinical_numerical","clinical_categorical"  
	#                 E.g. feature_types <- c("RNA,DNA")
	#                 E.g. feature_types <- c("RNA,clinical_numerical")
	#                 E.g. feature_types <- c("clinical_numerical,clinical_categorical")
	#                 E.g. feature_types <- c("RNA") 
	# OUTPUT
	# a data.frame in which the individual features are merged 
	

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

        }else if (i=="RNA_cor_zscore"){
            
            input <- read.table(paste(path,"prepared_data_rna_corgenes_zscore.csv",sep=""),sep="\t",header=TRUE,row.names=1)

        }else if (i=="RNA_cor_enriched_zscore"){
            
            input <- read.table(paste(path,"prepared_data_rna_corgenes_enriched_zscore.csv",sep=""),sep="\t",header=TRUE,row.names=1)            

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

launch_prediction <- function(features_model,drug_responses="all",predictors,str_predictors,nfolds=5 , nrepeats=10 , mc.cores=4,l_alpha=0,l_ntree=50,l_nodesize=20,
                              l_max.depth=2,l_eta=1,l_nthread=2,l_nrounds=20,l_k=5,hidden=c(50,20),algorithm="backprop",learningrate=1e-3, act.fct="tanh",l_kernel="radial"){

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
        
        if (!dir.exists(paste("./../prediction/results_SC1/",combination,sep=""))){
            dir.create(paste("./../prediction/results_SC1/",combination,sep=""))
        }
        
        if (!dir.exists(paste("./../prediction/results_SC1/",combination,"/stdoutput",sep=""))){
            dir.create(paste("./../prediction/results_SC1/",combination,"/stdoutput",sep=""))
        }
        
        outpath <- paste("./../prediction/results_SC1/",combination,sep="")
        str <- 1
        
        for (predictor in predictors){

            if(identical(predictorElasticNet,predictor)){
                for (alpha in l_alpha){
                    
                    filename <- paste(str_predictors[str],alpha,sep="_")
                    command <- paste(command_init, paste("'../prediction/results_SC1/",combination,"/stdoutput/",paste(filename,"_drugs_",paste(drug_responses,collapse="_"),sep=""),".txt' ",sep=""),
                                    "'Rscript evaluate_SC1_intermediate_parallel.R --predictor ",str_predictors[str]," --feature_model ", combination," --alpha ",alpha,
                                    " --nrepeats ",nrepeats," --nfolds ", nfolds, " --mc.cores ", mc.cores, " --outpath ",outpath, " --filename ",filename,"'",sep="") 
                    system(command)
                    print(command)
                }
            }

            if(identical(predictorRF,predictor)){
                for (ntree in l_ntree){
                    for (nodesize in l_nodesize){
                        
                        filename <- paste(str_predictors[str],ntree,nodesize,sep="_")
                        command <- paste(command_init, paste("'../prediction/results_SC1/",combination,"/stdoutput/",paste(filename,"_drugs_",paste(drug_responses,collapse="_"),sep=""),".txt' ",sep=""),
                                         "'Rscript evaluate_SC1_intermediate_parallel.R --predictor ",str_predictors[str]," --feature_model ", combination,
                                         " --ntree ",ntree, " --nodesize ",nodesize," --nrepeats ",nrepeats," --nfolds ", nfolds, " --mc.cores ", mc.cores, " --outpath ",outpath, " --filename ",filename,"'",sep="") 
                        system(command)
                        print(command)

                    }
                }
            }

            if(identical(predictorXGBoost,predictor)){
                for (max_depth in l_max.depth){
                    for (eta in l_eta){
                        for (nthread in l_nthread){
                            for (nrounds in l_nrounds){
                                
                                filename <- paste(str_predictors[str],max_depth,eta,nthread,nrounds,sep="_")
                                command <- paste(command_init, paste("'../prediction/results_SC1/",combination,"/stdoutput/",paste(filename,"_drugs_",paste(drug_responses,collapse="_"),sep=""),".txt' ",sep=""),
                                                 "'Rscript evaluate_SC1_intermediate_parallel.R --predictor ",str_predictors[str]," --feature_model ", 
                                                 combination," --max_depth ",max_depth, " --eta ",eta," --nthread ",nthread, " --nrounds ",nrounds, " --nrepeats ",1," --nfolds ", nfolds, " --mc.cores ", mc.cores, " --outpath ",outpath, 
                                                 " --filename ",filename,"'",sep="") 
                                system(command)
                                print(command)
                            }
                        }
                    }
                }
            }

            if(identical(predictorkNNReg,predictor)){
                for (k in l_k){

                    filename <- paste(str_predictors[str],k,sep="_") 
                    command <- paste(command_init, paste("'../prediction/results_SC1/",combination,"/stdoutput/",paste(filename,"_drugs_",paste(drug_responses,collapse="_"),sep=""),".txt' ",sep=""),
                                    "'Rscript evaluate_SC1_intermediate_parallel.R --predictor ",str_predictors[str]," --feature_model ", combination,
                                    " --k ",k, " --nrepeats ", nrepeats," --nfolds ", nfolds, " --mc.cores ", mc.cores, " --outpath ",outpath, " --filename ",filename,"'",sep="")
                    system(command)
                    print(command)
                }
            }  

            if(identical(predictorSVR,predictor)){
                for (kernel in l_kernel){
                    
                    filename <- paste(str_predictors[str],kernel,sep="_")
                    command <- paste(command_init, paste("'../prediction/results_SC1/",combination,"/stdoutput/",paste(filename,"_drugs_",paste(drug_responses,collapse="_"),sep=""),".txt' ",sep=""),
                                    "'Rscript evaluate_SC1_intermediate_parallel.R --predictor ",str_predictors[str]," --feature_model ", combination,
                                    " --kernel ",kernel, " --nrepeats ", nrepeats," --nfolds ", nfolds, " --mc.cores ", mc.cores, " --outpath ",outpath, " --filename ",filename,"'",sep="") 
                    system(command)
                    print(command)
                }
            } 
            
            str <- str + 1
        }
    }
}

                   
best_model_per_drug <- function(drug, choose = "all",algo = "all"){
    
    # INPUT drug : the name of a drug
    # OUTPUT : list declaring the maximum correlation and the corresponding combination of features and model

    path <- "./../prediction/results_SC1"
    
    # retrieve combinations
    combinations <- list.dirs(path = path, full.names = FALSE, recursive = FALSE)
    combinations <- combinations[which(combinations!='.ipynb_checkpoints')]
    
    if (choose!="all"){
        
        combinations <- combinations[which(combinations %in% choose)]
    }
    
    max_corr <- rep(-1,length(combinations))
    min_rmsd <- rep(10000,length(combinations))
    best_model <- rep("none",length(combinations))
    best_model_rmsd <- rep("none",length(combinations))
    
    names(max_corr) <- combinations
    names(min_rmsd) <- combinations
    names(best_model) <- combinations
    names(best_model_rmsd) <- combinations
    
    for (comb in combinations){
        
        comb_path <- paste(path,"/",comb,sep="")
        files_tosearch <- list.files(path = comb_path, pattern = "\\.csv$")
        
        if (algo!="all"){
        
            files_tosearch <- grep(paste(algo,collapse="|"),files_tosearch,value=TRUE)
        }
        
        files_tosearch <- c(grep("drugs_all",files_tosearch,value=TRUE), grep(drug,files_tosearch,value=TRUE))
        
        for (model in files_tosearch){
            
            model_path <- paste(comb_path,"/",model,sep="")
            res <- read.table(model_path)[drug,"meanSpearman"]
            res_rmsd <- read.table(model_path)[drug,"meanRMSD"]
            
            if (is.na(res)){ # hardcode NA results
                  res <- -1
            }
            if (is.na(res_rmsd)){ # hardcode NA results
                  res_rmsd <- 10000
            }
            
            if (res > max_corr[comb]){
                
                max_corr[comb] <- res
                best_model[comb] <- model
            }
            
            if (res_rmsd < min_rmsd[comb]){
                
                min_rmsd[comb] <- res_rmsd
                best_model_rmsd[comb] <- model
            }
        } 
    } 
    
    corr <- max(max_corr)
    rmsd <- min(min_rmsd)
    
    if (length(which(max_corr==corr))>1){ #in case 2 feature combinations have the same performance keep the first one
        
        model <- best_model[which(max_corr==corr)[1]]
        comb <- names(max_corr)[which(max_corr==corr)[1]]        
        
    }else{
    
        model <- best_model[which(max_corr==corr)]
        comb <- names(max_corr)[which(max_corr==corr)]
    }
    
    df <- data.frame("best_features"=comb,"best_model"=model,"max_corr"=corr)
    rownames(df) <- drug
    
    if (length(which(min_rmsd==rmsd))>1){ #in case 2 feature combinations have the same performance keep the first one
        
        model <- best_model_rmsd[which(min_rmsd==rmsd)[1]]
        comb <- names(min_rmsd)[which(min_rmsd==rmsd)[1]]        
        
    }else{
    
        model <- best_model_rmsd[which(min_rmsd==rmsd)]
        comb <- names(min_rmsd)[which(min_rmsd==rmsd)] 
    }
    
    df_rmsd <- data.frame("best_features"=comb,"best_model"=model,"min_rmsd"=rmsd)
    rownames(df_rmsd) <- drug
    
    toreturn <- list(df,df_rmsd)
    names(toreturn) <- c("spearman","rmsd")
    
    return(toreturn)
    
}

save_trained_best_models <- function (data_response,data_matrix,path="../docker/docker/final_models/loaded_"){

    ### Require : combine features
    ###
    ### Input : - data_response : all drug responses
    ###         - data_matrix : best models file containing best model and best feature combination for each drug
    ###
    
    ### Output : One RDS File by drug response

    
    ### Save RDS File function
    save_rds <- function(final_model,drug_name){
        saveRDS(final_model, paste(path,drug_name,".rds",sep=""))
    }
    ### Save RDS File function
    print(data_response)
    print(data_matrix)
    for (i in 1:nrow(data_matrix)){
        
        features <- as.list(strsplit(as.character(data_matrix$best_features[i]), "with")[[1]])
        celllines <- data.matrix(combine_features(features)) 
        response <- data.matrix(data_response[,rownames(data_matrix)[i],drop=F])

        ### Handle Nas
        non_na <- intersect(rownames(na.omit(celllines)),rownames(na.omit(response)))
        celllines <- celllines[non_na,]
        response <- response[non_na,,drop=F]
        ###

        hyper_params <- as.list(strsplit(as.character(data_matrix$best_model[i]), "_")[[1]])

        if (hyper_params[[1]]=="predictorRF"){

            nntree <- as.numeric(hyper_params[[2]])
            nodesize <- as.numeric(hyper_params[[3]])
            final_model <- randomForest(celllines, response, ntree = nntree, nodesize = nodesize)

            save_rds(final_model,rownames(data_matrix)[i])

        } else if (hyper_params[[1]]=="predictorElasticNet"){

            alpha <- as.numeric(hyper_params[[2]])

            rescv <- cv.glmnet(celllines, response , alpha = alpha)
            final_model <- glmnet(celllines, response , alpha = alpha, lambda = rescv$lambda.min)

            save_rds(final_model,rownames(data_matrix)[i])

        } else if (hyper_params[[1]]=="predictorkNNReg"){

            k <- as.numeric(hyper_params[[2]])

            final_model <- knnreg(celllines, response , k=k)

            save_rds(final_model,rownames(data_matrix)[i])

        } else if (hyper_params[[1]]=="predictorSVR"){

            kernel <- as.character(hyper_params[[2]])

            final_model <- tune(svm, type = "eps-regression", kernel =  kernel , epsilon = 0.1, train.x=celllines, train.y=response, 
                                  ranges = list(gamma = 10^(-5:1), cost = 10^(-2:4)),
                                  tunecontrol = tune.control(sampling = "cross", cross = 10))

            save_rds(final_model$best.model,rownames(data_matrix)[i])

        } else if (hyper_params[[1]]=="predictorXGBoost"){
            
            max.depth <- as.numeric(hyper_params[[2]])
            eta <- as.numeric(hyper_params[[3]]) 
            nthread <- as.numeric(hyper_params[[4]]) 
            nrounds <- as.numeric(hyper_params[[5]])

            final_model <- xgboost(data = celllines, label = response , max.depth = max.depth,
                                   eta = eta, nthread = nthread, nrounds = nrounds, objective = "reg:squarederror")
            save_rds(final_model,rownames(data_matrix)[i])
        }
        print(class(final_model))

    }
}
