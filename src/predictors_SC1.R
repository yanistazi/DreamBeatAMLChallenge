# In this file we have predictors functions, that takes at least 4 arguments : 
#		- celllinesTrain : n*d matrix of descriptors for n training cell lines
#		- celllinesTest : m*d matrix of descriptors for m test cell lines
#		- response : n*p matrix of response (i.e. drug AUC) values for the p chemicals on the n training cell lines
# and that outputs a m*p matrix of predicted response values for the p chemicals on the m test cell lines
#
# Please choose a function name that starts with "predictor'
# You can create additional files of the form predictors*.R to code more functions if you want


### 
library(caret)
#library(neuralnet)
library(xgboost)

## Zero predictor : always output zero
predictorZero <- function(celllinesTrain , celllinesTest , response) {
	ncelllinesTest <- dim(celllinesTest)[1]
	nchemicals <- dim(response)[2]
	pred <- matrix(data=0,nrow= ncelllinesTest,ncol= nchemicals)
	dimnames(pred) <- list(rownames(celllinesTest),rownames(chemicals))
	return(pred)
}

## Mean predictor : predict the mean response of each chemical
predictorMean <- function(celllinesTrain , celllinesTest , chemicals , response) {
	ncelllinesTest <- dim(celllinesTest)[1]
	#nchemicals <- dim(response)[2]
    nchemicals <- 1
	m <-apply(response,2,mean)
	pred <- t(matrix(data=m,nrow=nchemicals,ncol=ncelllinesTest))
	dimnames(pred) <- list(rownames(celllinesTest),rownames(chemicals))
	return(pred)
}

library(glmnet)
## ELASTIC-NET regression for each chemical (without using chemical features)
predictorElasticNet <- function(celllinesTrain , celllinesTest , response, alpha=0) {
	# alpha = 1 is LASSO
        # alpha = 0 is RIDGE

	# Train a lasso model
	cvob1=cv.glmnet(data.matrix(celllinesTrain), response , alpha = alpha)
	# Make predictions
	pred <- predict(cvob1,data.matrix(celllinesTest), s="lambda.min")


	pred1 <- data.frame(pred)
	#names(pred1) <- colnames(response)
	return(pred1)
}



## Random forest regression for each chemical (without using chemical features)
library(randomForest)
predictorRF <- function(celllinesTrain , celllinesTest , response, nntree,nodesize = 5) { # adding nodesize as well with default is 5

    # Train a random forest model
    model_rf=randomForest(celllinesTrain, response, ntree = nntree, nodesize = nodesize)
    # Make predictions
    pred <- predict(model_rf,celllinesTest)

	

	pred1 <- data.frame(pred)
	#names(pred1) <- colnames(response) 
	return(pred1)
}


## SV regression for each chemical (without using chemical features) based on a gaussian kernel with parameter gamma : 
## exp(- gamma |x-y|^2)
## The threshold epsilon is set to 0.1 but can be changed in the program
## Parameters gamma and cost are determined by cross validation
library(e1071)	

# add other kernels in svm  : polynomial and tanh
predictorSVR <- function(celllinesTrain , celllinesTest , response,kernel = "radial") {


		# Train svm by cross validation over parameters gamma and cost
	model_svm <- tune(svm, type = "eps-regression", kernel =  kernel , epsilon = 0.1, celllinesTrain, response, 
              				ranges = list(gamma = 10^(-5:1), cost = 10^(-2:4)),
              						tunecontrol = tune.control(sampling = "cross", cross = 10))
	# Make predictions
	pred <- predict(model_svm$best.model, celllinesTest)


	pred1 <- data.frame(pred)
	#names(pred1) <- colnames(response) 
	return(pred1)
}



predictorXGBoost <- function(celllinesTrain , celllinesTest , response, max.depth = 2, eta = 1, nthread = 2, nrounds = 20) {

	model_xgboost = xgboost(data = data.matrix(celllinesTrain), label = response , max.depth = max.depth, eta = eta, nthread = nthread, nrounds = nrounds, objective = "reg:squarederror")
	# Make predictions
	pred <- predict(model_xgboost, data.matrix(celllinesTest))


	pred1 <- data.frame(pred)
    row.names(pred1) <- rownames(celllinesTest)
	names(pred1) <- colnames(response) 
	return(pred1)
}



predictorkNNReg <- function(celllinesTrain , celllinesTest , response, k=5) {

	# Train kNN
	model_knn = knnreg(data.matrix(celllinesTrain), response , k=k)
	# Make predictions
	pred <- predict(model_knn, data.matrix(celllinesTest))


	pred1 <- data.frame(pred)
    row.names(pred1) <- rownames(celllinesTest)
	names(pred1) <- colnames(response) 
	return(pred1)
}

predictorNeuralNet <- function(celllinesTrain , celllinesTest , response, hidden =c(10,10),algorithm="backprop",learningrate =1e-3,act.fct = "tanh",
                linear.output = T) {

	# Train NeuralNet
    data <- cbind(celllinesTrain,responseTrain = response)
	model_NeuralNet = neuralnet(responseTrain~.,data=data, hidden =hidden,algorithm=algorithm,learningrate =learningrate ,act.fct = act.fct)
	# Make predictions
	pred <- compute(model_NeuralNet, celllinesTest)$net.result


	pred1 <- data.frame(pred)
	return(pred1)
}