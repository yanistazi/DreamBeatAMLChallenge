library(Hmisc)
library(parallel)

KevaluateCV <- function(myKpredictor , celllines , response , KernelTask , KpatientGiven=FALSE , KpatientType="linear" , nfolds=5 , nrepeats=10 , seed=9182456 , mc.cores=1, ...) {
   # Evaluate the performance of a predictor by cross-validation
   #
   # INPUT
   # mypredictor : a multitask Kernel predictor function that takes at least 5 arguments : 
   #		- celllinesTrain : matrix of descriptors for n chemicals (with n rows)
   #		- chemicals : matrix of descriptors for p chemicals (with p rows)
   #		- response : n*p matrix of response values for the p chemicals on the n training cell lines
   # 		and that outputs a m*p matrix of predicted toxicities for the p chemicals on the m test cell lines
   # celllines : a matrix of descriptors for celllines (1 row = 1 cell line)
   # response : a matrix of response values for all cell lines (in row) and all drugs/tasks (in column)
   # KernelTask : Kernel Task
   # nfolds : number of CV folds
   # nrepeats : number of times the cross-validation is repeated
   # seed : a seed number for the random number generator (use same seed to have the same CV splits)
   #
   # OUTPUT
   # a list with the performance of the method

   # Set random number generator seed
   set.seed(seed)


   # Number of chemicals
   nchemicals <- dim(response)[2]
   chemicals <- colnames(response)

   # Iterate for each drug 
   all_drugs_result <- list()
   for (dr in 1:nchemicals){

      # Remove the patients that have NA in the response for this specific drug
      drug <- chemicals[dr]

      print(drug)

      i_non_na <- which(!is.na(response[,drug]))
      celllines_dr <- celllines[i_non_na,]

      if (KpatientGiven) { # In that case celllines is a ncell x ncell similarity matrix
	 celllines_dr <- celllines_dr[,i_non_na]
      }

      # Keep all others drugs for the non-na drug patients to allow for information sharing
      drug_response <- response[i_non_na,,drop=F]
      ncelllines <- dim(celllines_dr)[1]

      # Make folds for the specific drug, based on the non_na 
      n <- ncelllines
      folds <- list()
      for (i in seq(nrepeats)) {
        folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
      }
      nexp <- length(folds)

      # Calculate Drug Correlation Kernel ...
      # WHAT WE NEED: recovery of symmetric positive semi-definite matrices [correlation matrix]
      # ie [correlation matrix completion problem]
      #kernelTask = cor(response, use="pairwise.complete.obs")
      #kernelTask = cor(response, use="pairwise.complete.obs")
      #zaza=CLRMC(kernelTask)
      #kernelTask = rcorr(as.matrix(response), type="spearman")

      ## Main CV loop (parallelized):
      resCV <- mclapply(seq(nexp), function(iexp){
			   cat('.')
			   mytrain = seq(n)[-folds[[iexp]]]
			   mytest = folds[[iexp]]

			   celllinesTrain <- celllines_dr[mytrain,]
			   celllinesTest <- celllines_dr[mytest,]

			   responseTrain <- drug_response[mytrain,]
			   responseTest <- drug_response[mytest,dr]

			   if (KpatientGiven) {
			      celllinesTrain <- celllinesTrain[,mytrain]
			      celllinesTest <- celllinesTest[,mytrain]
			      KpatientType <- "precomputed"
			   }
			   
			   # Calculate Drug Correlation Kernel on the Training Set
			   #kernelTask = cor(responseTrain, use="pairwise.complete.obs")
			   
			   responsePred <- myKpredictor(celllinesTrain = celllinesTrain ,
							celllinesTest = celllinesTest  ,
							responseTrain = responseTrain ,
							task.number = dr ,
							KernelTask = KernelTask, kx_type=KpatientType, ...)
			   #print(responsePred)

			   responseError <- responseTest - responsePred

			   #print(responseError)
			   # ~~~~HACK~~~~~ #
			   # if responsePred is a constant vector (which can happens for Lasso for examples...)
			   # then correlations returns NA ...
			   if (length(unique(responsePred))==1) {
			      print("Your shit is constant")
			      responsePred[1] = responsePred[1] + 0.1			   
			   }
			   responseSpearman <- cor(responsePred,responseTest,method = 'spearman')
			   # ~~~~~~~~~~~~~
			   res <- numeric(3)
			   res[1] <- sqrt(sum((responseError)^2)/(length(responseError)))
			   #res[2] <- mean(responseCI)
			   res[2] <- 0
			   res[3] <- responseSpearman
			   res
}, 
mc.cores=mc.cores)
      resCV <- unlist(resCV)
      err <- resCV[seq(1,nexp*3,by=3)]
      #print(err)
      corInd <- resCV[seq(2,nexp*3,by=3)]
      corrSpearman <- resCV[seq(3,nexp*3,by=3)]
      all_drugs_result[[drug]] <- c(meanRMSD=mean(err) , sdRMSD=sd(err) , meanCI=mean(corInd) , meanSpearman=mean(corrSpearman) , sdSpearman=sd(corrSpearman) , minSpearman=min(corrSpearman))
   }

   return(all_drugs_result)
}
