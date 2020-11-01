library(glmnet)
library(parallel)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)


predictorGLM <- function(celllinesTrain , celllinesTest , response, alpha=0, ninternalfolds=5) {
    # alpha=1 --> l1 penalty
    # alpha=0 --> l2 penalty
    set.seed(17)
    # Train
    cvfit = cv.glmnet(data.matrix(celllinesTrain), data.matrix(response), family="cox", alpha=alpha, nfolds=ninternalfolds, grouped=TRUE)
    # Predict
    risk.predict = predict(cvfit, newx=data.matrix(celllinesTest), s="lambda.min", type="response")
    risk.predict = as.vector(risk.predict[,1])

    return(risk.predict)
}    

predictorBoost<-function(celllinesTrain , celllinesTest , response , maxstepno=500 , K=10 , type="verweij" , penalty=100 ){
    set.seed(17)
    

    cv.res<-cv.CoxBoost(time=as.numeric(unlist(response[,1])),
                  status=as.numeric(unlist(response[,2])),
                  x=data.matrix(celllinesTrain),maxstepno=maxstepno,K=K, type=type , penalty=penalty )
    
    cvfit <- CoxBoost(time=as.numeric(unlist(response[,1])),
                  status=as.numeric(unlist(response[,2])),
                  x=data.matrix(celllinesTrain),stepno=cv.res$optimal.step , penalty=penalty)

  
    risk.predict<-predict(cvfit,data.matrix(celllinesTest),newtime=as.numeric(unlist(response[,1])),newstatus=as.numeric(unlist(response[,2])),stepno=1000)
  
  return(as.vector(risk.predict))
}

predictorRF <- function(celllinesTrain , celllinesTest , response, nntree=10, nodesize=10) {
    set.seed(17 )
    # Train
    cvfit = rfsrc(Surv(time, status) ~ ., data=data.frame(celllinesTrain,response), ntree=nntree, nodesize=nodesize, importance="none")
    
    # Predict
    risk.predict = predict(cvfit, data.frame(celllinesTest), importance="none")$predicted
    
    return(risk.predict)
} 

predictorRFX <- function(celllinesTrain, celllinesTest, response, max.iter = 500,tol=0.01,groups = rep(1, ncol(celllinesTrain))) {
    set.seed(17)
    # Train
    cvfit = CoxRFX(data.frame(celllinesTrain), Surv(time=response[,1],event =response[,2]) , max.iter =max.iter,tol=tol,groups=groups)
    cvfit$Z <- NULL
    # Predict
    risk.predict<-predict(cvfit,data.frame(celllinesTest))
    
    return(risk.predict)
}
