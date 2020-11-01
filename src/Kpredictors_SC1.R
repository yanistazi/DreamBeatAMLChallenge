library(kmr)

KpredictorKMR <- function(celllinesTrain, celllinesTest, responseTrain,
			  task.number,
			  KernelTask,
			  kx_type = "linear", # linear gaussian precomputed
			  kx_option = list(sigma = 1), # for gaussian kernel
			  type.measure = "ci", # ci mse cor
			  nfolds.intern = 7,
			  nrepeats.intern = 4,
			  lambdas = 10^(-10:10)
			  ) {
   # Input:
   # celllinesTrain: ntrain patients x d descriptors
   # celllines Test: ntest patients x d descriptors
   # responseTrain: ntrain x ntask responseTrain (drug responseTrains) values
   # task.number: task for which the KMR prediction is returned
   # KernelTask ntask x ntask e.g. correlation matrix between the tasks
   # Output:
   # Vector of predicted values for ntest patients on the task #task.number

   # We assume that responseTrain[,task.number] does not contain NA
   if (any(is.na(responseTrain[,task.number]))) {
      print("Response for given task contains NA") }

   # Input NA if needed for other task(s) than task.number
   # ByRow ...  
   # TODO: improve?
   for (j in 1:nrow(responseTrain)) {
      i = which(is.na(responseTrain[j,]))
      if (length(i)>0) {
	 responseTrain[j,i] = mean(as.numeric(responseTrain[j,]),na.rm=T)
      }
   }

   #print("Before")
   # KMR cross-validation for regularization
   rescv = cv.kmr(x=celllinesTrain, y=responseTrain,
		  kx_type = kx_type,
		  kx_option = kx_option,
		  kt_type = "precomputed",
		  kt_option = list(kt=KernelTask),
		  type.measure = type.measure,
		  nfolds = nfolds.intern,
		  nrepeats = nrepeats.intern,
		  lambda = lambdas
   )
   #print("After")
   
   # Best Lambda for the task
   mylambda = rescv$bestlambda[task.number]
   #print(mylambda)

   # Train a KMR
   res = kmr(x=celllinesTrain, y=responseTrain,
	     kx_type = kx_type,
	     kx_option = kx_option,
	     kt_type = "precomputed",
	     kt_option = list(kt=KernelTask)
   )
   # Predict
   pred = predict(res, as.matrix(celllinesTest), lambda=mylambda)
   # Return
   return(pred[,task.number])
}



KtrainKMR <- function(celllinesTrain, responseTrain,
			  task.number,
			  KernelTask,
			  kx_type = "linear", # linear gaussian
			  kx_option = list(sigma = 1), # for gaussian kernel
			  type.measure = "ci", # ci mse cor
			  nfolds.intern = 7,
			  nrepeats.intern = 4,
			  lambdas = 10^(-10:10)
			  ) {
   # Input:
   # celllinesTrain: ntrain patients x d descriptors
   # responseTrain: ntrain x ntask responseTrain (drug responseTrains) values
   # task.number: task for which the KMR prediction is returned
   # KernelTask ntask x ntask e.g. correlation matrix between the tasks
   # Output:
   # Vector of predicted values for ntest patients on the task #task.number

   # We assume that responseTrain[,task.number] does not contain NA
   if (any(is.na(responseTrain[,task.number]))) {
      print("Response for given task contains NA") }

   # Input NA if needed for other task(s) than task.number
   # ByRow ...  
   # TODO: improve?
   for (j in 1:nrow(responseTrain)) {
      i = which(is.na(responseTrain[j,]))
      if (length(i)>0) {
	 responseTrain[j,i] = mean(as.numeric(responseTrain[j,]),na.rm=T)
      }
   }

   #print("Before")
   # KMR cross-validation for regularization
   rescv = cv.kmr(x=celllinesTrain, y=responseTrain,
		  kx_type = kx_type,
		  kx_option = kx_option,
		  kt_type = "precomputed",
		  kt_option = list(kt=KernelTask),
		  type.measure = type.measure,
		  nfolds = nfolds.intern,
		  nrepeats = nrepeats.intern,
		  lambda = lambdas
   )
   #print("After")
   
   # Best Lambda for the task
   mylambda = rescv$bestlambda[task.number]
   #print(mylambda)

   # Train a KMR
   res = kmr(x=celllinesTrain, y=responseTrain,
	     kx_type = kx_type,
	     kx_option = kx_option,
	     kt_type = "precomputed",
	     kt_option = list(kt=KernelTask)
   )
   
   # Return
   return(list(object=res, lambda=mylambda))
}

