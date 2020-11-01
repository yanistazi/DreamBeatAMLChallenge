evaluateCV <- function(mypredictor , celllines , response , nfolds=5 , nrepeats=10 , seed=9182456 , mc.cores=1, ...) {
   # Evaluate the performance of a predictor by cross-validation
   #
   # INPUT
   # mypredictor : a predictor function that takes at least 4 arguments : 
   #		- celllines : matrix of descriptors for n chemicals (with n rows)
   #		- chemicals : matrix of descriptors for p chemicals (with p rows)
   #		- response : n*p matrix of response values for the p chemicals on the n training cell lines
   # 		and that outputs a m*p matrix of predicted toxicities for the p chemicals on the m test cell lines
   # celllines : a matrix of descriptors for cell lines (1 row = 1 cell line)
   # response : a matrix of response values for all cell lines (in row) and all chemicals (in column)
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
      non_na <- intersect(rownames(na.omit(celllines)),rownames(na.omit(response[drug])))
      celllines_dr <- celllines[non_na,]
      drug_response <- response[drug][non_na,,drop=F]
      ncelllines <- dim(celllines_dr)[1] # change place of ncelllines to match folds
      # Make folds for the specific drug, based on the non_na 
      n <- ncelllines
      folds <- list()
      for (i in seq(nrepeats)) {
        folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
      }
      nexp <- length(folds)

      ## Main CV loop (parallelized):
      resCV <- mclapply(seq(nexp), function(iexp){
			   cat('.')
			   celllinesTrain <- celllines_dr[-folds[[iexp]],]
			   celllinesTest <- celllines_dr[folds[[iexp]],]
			   responseTrain <- drug_response[-folds[[iexp]],]
			   responseTest <- drug_response[folds[[iexp]],,drop=F]
			   responsePred <- mypredictor(celllinesTrain , celllinesTest  , responseTrain, ...)
			   #print(responsePred)
			   responseError <- responseTest - responsePred
			   #responseCI <- apply(rbind(responsePred[,1], responseTest), 2, function(u){
			   # rcorr.cens(u[1:nrow(responsePred)], u[-(1:nrow(responsePred))], outx=FALSE)[[1]]})
			   # ~~~~HACK~~~~~
			   # if responsePred is a constant vector (which can happens for Lasso for examples...)
			   # then correlations returns NA ...
			   if (length(unique(responsePred[,1]))==1) {
			      print("zaza")
			      responsePred[1,1] = responsePred[1,1] + 0.1			   
			   }
			   responseSpearman <- compute_spearman(responsePred,responseTest)
			   # ~~~~~~~~~~~~~
			   res <- numeric(3)
			   res[1] <- sqrt(sum((responseError)^2)/(dim(responseError)[1] * dim(responseError)[2]))
			   #res[2] <- mean(responseCI)
			   res[2] <- 0
			   res[3] <- responseSpearman[["rho"]]
			   res
}, 
mc.cores=mc.cores)
      resCV <- unlist(resCV)
      err <- resCV[seq(1,nexp*3,by=3)]
      #print(err)
      corInd <- resCV[seq(2,nexp*3,by=3)]
      corrSpearman <- resCV[seq(3,nexp*3,by=3)]
      all_drugs_result[[drug]] <- c(meanRMSD=mean(err), sdRMSD=sd(err), meanCI=mean(corInd), meanSpearman=mean(corrSpearman), sdSpearman=sd(corrSpearman))
   }

   return(all_drugs_result)
}
