evaluateCV_SC2_ensemble_transfer_learning <- function(mypredictor_challenge, mypredictor_challenge_MSK, 
						      celllines_challenge, celllines_challenge_MSK, 
						      response,weight_challenge=0.5,
						      MSK.path = "../data/features_SC2/",
						      nfolds=5, nrepeats=10, seed=17,
						      mc.cores=1,use_MSK=0,
						      alpha_challenge=0,ntree_challenge=500,nodesize_challenge=10,maxstepno_challenge=500,
						      K_challenge=10,type_challenge="verweij",penalty_challenge=100,max.iter_challenge = 500,
						      tol_challenge=0.01,
						      alpha_challenge_MSK=0,ntree_challenge_MSK=500,nodesize_challenge_MSK=10,maxstepno_challenge_MSK=500,
						      K_challenge_MSK=10,type_challenge_MSK="verweij",penalty_challenge_MSK=100,
						      max.iter_challenge_MSK = 500,tol_challenge_MSK=0.01,                       
						      ...) {

   response <- na.omit(response)
   response <- response[response$OS>0,]
   non_na <- intersect(intersect(rownames(na.omit(celllines_challenge)),rownames(na.omit(response))),rownames(na.omit(celllines_challenge_MSK)))
   celllines_challenge <- celllines_challenge[non_na,]
   celllines_challenge_MSK <- celllines_challenge_MSK[non_na,]
   response <- response[non_na,]
   response$time <- response$OS 
   response$status <- response$OS_status
   response$OS <- NULL
   response$OS_status <- NULL
   set.seed(seed)
   # Make folds
   n = nrow(celllines_challenge)
   folds <- list()
   for (i in seq(nrepeats)) {
      folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
   }
   nexp = length(folds) # the total number CV of experiments
   if("NPM1" %in% colnames(celllines_challenge_MSK)){
      MSK_data <- read.table(paste0(MSK.path,"similarity_ordered_prepared_data_MSK.csv"))
   }else{
      MSK_data <- read.table(paste0(MSK.path,"similarity_ordered_prepared_data_MSK_reduced.csv"))
   }
   MSK_train <- MSK_data[,colnames(celllines_challenge_MSK)]

   MSK_response <- MSK_data[,names(MSK_data) %in%c("OS","OS_status")]
   MSK_response$time <- MSK_response$OS
   MSK_response$status <- MSK_response$OS_status
   MSK_response$OS <- NULL
   MSK_response$OS_status <- NULL
   # Parallel CV
   print("start CV")
   rescv = mclapply(seq(nexp),
		    FUN=function(iexp) {
		       cat(".")
		       celllinesTrain_challenge = celllines_challenge[-folds[[iexp]],,drop=F]
		       celllinesTest_challenge = celllines_challenge[folds[[iexp]],,drop=F]

		       celllinesTrain_challenge_MSK = celllines_challenge_MSK[-folds[[iexp]],,drop=F]
		       celllinesTest_challenge_MSK = celllines_challenge_MSK[folds[[iexp]],,drop=F]

		       responseTrain = response[-folds[[iexp]],]
		       responseTest = response[folds[[iexp]],]                     
		       # Train and Predcit
		       if(use_MSK>0){
			  num_rows <- round(nrow(MSK_data)*use_MSK/100)  # num of row to use based on use_MSK percentage 
			  celllinesTrain_challenge_MSK <- rbind(celllinesTrain_challenge_MSK,MSK_train[1:num_rows,])
			  responseTrain_challenge_MSK <- rbind(responseTrain,MSK_response[1:num_rows,])
		       }


		       if(identical(predictorGLM,mypredictor_challenge)){
			  predict.test_challenge <- rank(mypredictor_challenge(celllinesTrain=celllinesTrain_challenge,
									       celllinesTest=celllinesTest_challenge, response=responseTrain, alpha=alpha_challenge))
		       } else if(identical(predictorRF,mypredictor_challenge)){
			  predict.test_challenge <- rank(mypredictor_challenge(celllinesTrain=celllinesTrain_challenge,
									       celllinesTest=celllinesTest_challenge, response=responseTrain, 
									       nntree=ntree_challenge,nodesize=nodesize_challenge))
		       } else if(identical(predictorRFX,mypredictor_challenge)){
			  predict.test_challenge <- rank(mypredictor_challenge(celllinesTrain=celllinesTrain_challenge,
									       celllinesTest=celllinesTest_challenge, response=responseTrain, 
									       max.iter=max.iter_challenge,tol=tol_challenge))
		       } else if(identical(predictorBoost,mypredictor_challenge)){
			  predict.test_challenge <- rank(mypredictor_challenge(celllinesTrain=celllinesTrain_challenge,
									       celllinesTest=celllinesTest_challenge, response=responseTrain,
									       maxstepno=maxstepno_challenge,K=K_challenge, type=type_challenge , penalty=penalty_challenge))  
		       }

		       # same for predictor challenge MSK    

		       if(identical(predictorGLM,mypredictor_challenge_MSK)){
			  predict.test_challenge_MSK <- rank(mypredictor_challenge_MSK(celllinesTrain=celllinesTrain_challenge_MSK,
										       celllinesTest=celllinesTest_challenge_MSK, response=responseTrain_challenge_MSK,
										       alpha=alpha_challenge_MSK))
		       } else if(identical(predictorRF,mypredictor_challenge_MSK)){
			  predict.test_challenge_MSK <- rank(mypredictor_challenge_MSK(celllinesTrain=celllinesTrain_challenge_MSK,
										       celllinesTest=celllinesTest_challenge_MSK, response=responseTrain_challenge_MSK,
										       nntree=ntree_challenge_MSK,nodesize=nodesize_challenge_MSK))
		       } else if(identical(predictorRFX,mypredictor_challenge_MSK)){
			  predict.test_challenge_MSK <- rank(mypredictor_challenge_MSK(celllinesTrain=celllinesTrain_challenge_MSK,
										       celllinesTest=celllinesTest_challenge_MSK, response=responseTrain_challenge_MSK, 
										       max.iter=max.iter_challenge_MSK,tol=tol_challenge_MSK))
		       } else if(identical(predictorBoost,mypredictor_challenge_MSK)){
			  predict.test_challenge_MSK <- rank(mypredictor_challenge_MSK(celllinesTrain=celllinesTrain_challenge_MSK,
										       celllinesTest=celllinesTest_challenge_MSK, response=responseTrain_challenge_MSK,
										       maxstepno=maxstepno_challenge_MSK,K=K_challenge_MSK, type=type_challenge_MSK , penalty=penalty_challenge_MSK))  
		       }    

		       predict.test <- weight_challenge * predict.test_challenge + (1-weight_challenge) * predict.test_challenge_MSK

		       # Evaluate CI on the test
		       ci.test = suppressWarnings(survConcordance(Surv(time,status) ~ predict.test, as.data.frame(responseTest)))
		       return(as.vector(ci.test$concordance))
		    },
		    mc.cores=mc.cores )


   return(unlist(rescv))

}
