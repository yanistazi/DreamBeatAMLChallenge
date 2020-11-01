source("combine_features.R")

library(glmnet)
#library(randomForestSRC)
#library(CoxBoost)
#library(CoxHD)


predict_from_best_model <- function(rds_path="final_model/", path="/features/", action.write=TRUE) {

   rdsfile <- list.files(rds_path)
   model <-  readRDS(paste(rds_path,rdsfile,sep="/"))

   len <- length(as.list(strsplit(as.character(rdsfile), "with")[[1]]))
   features <- as.list(strsplit(as.character(rdsfile), "with")[[1]])[1:(len-1)]
   predictor <- strsplit(as.character(rdsfile), "with")[[1]][len]
   features <- combine_features(features,path=path)

   print(dim(features))

   # NB : sign minus for Boost RFX and GLM because the outputs correspond to beta*
   if(grepl("predictorRF_",predictor)){

      mypred = -predict(model,features)$predicted 

   } else if (grepl("predictorRFX_",predictor)){

      mypred = - predict(model,features)

   } else if (grepl("predictorBoost",predictor)){

      mypred = t( - predict(model,features)) 

   } else if (grepl("predictorGLM",predictor)){

      mypred = - predict(model,newx=data.matrix(features), s="lambda.min", type="response")[,1]
   }

   dres = data.frame(lab_id = rownames(features) , survival = mypred)

   print(head(dres,10))
   if (action.write) {
      write.table(dres,"/output/predictions.csv",sep=",",row.names=FALSE,quote=FALSE)
   }

   return(dres)

}


predict_ensemble_learning <- function(rds_challenge_path="final_model_challenge/",
				      rds_MSKplusChallenge_path="final_model_MSK_challenge/",
				      path.predict="/features/",
				      weight_challenge=0.75) {

   df_challenge <- predict_from_best_model(rds_challenge_path,path=path.predict,action.write=FALSE)
   colnames(df_challenge) <- c("lab_id","survival_challenge")
   df_challenge$rank_challenge = rank(df_challenge$survival)

   df_MSK_challenge <- predict_from_best_model(rds_MSKplusChallenge_path,path=path.predict,action.write=FALSE)[,2,drop=F]
   colnames(df_MSK_challenge) <- "survival_MSKplusChallenge"
   df_MSK_challenge$rank_MSKplusChallenge = rank(df_MSK_challenge$survival)

   df <- cbind(df_challenge,df_MSK_challenge)
   df$survival <- weight_challenge * df$rank_challenge + (1-weight_challenge) * df$rank_MSKplusChallenge
   df$survival <- df$survival + rnorm(nrow(df),0,0.1)   ### handle equal values by adding random gaussian noise with mean 0 and var 0.1
   
   dres = df[,c("lab_id","survival")]

   write.table(dres,"/output/predictions.csv",sep=",",row.names=FALSE,quote=FALSE)

   return()

}
