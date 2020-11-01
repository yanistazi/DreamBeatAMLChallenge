evaluateCV_SC2 <- function(mypredictor, celllines, response, nfolds=5, nrepeats=10, seed=17,
                           mc.cores=1,use_MSK=0, ...) {

    
    # use_MSK : percentage of MSK data to use (the dataframe is ordered based on  the similarity with challenge data)
    # similar as SC1 for concordance index
    #
    # output a list of size the number of CV experiments (eg 50) (= nfolds x nrepeats)
    
    # Remove the patients that have NA in survival or that have time =0 in survival ! (from 213 to 193)
    response <- na.omit(response)
    response <- response[response$OS>0,]
    non_na <- intersect(rownames(na.omit(celllines)),rownames(na.omit(response)))
    celllines <- celllines[non_na,]
    response <- response[non_na,]
    response$time <- response$OS 
    response$status <- response$OS_status
    response$OS <- NULL
    response$OS_status <- NULL
    set.seed(seed)
    # Make folds
    n = nrow(celllines)
    folds <- list()
    for (i in seq(nrepeats)) {
        folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
    }
    nexp = length(folds) # the total number CV of experiments

    if (use_MSK>0) {
       if("NPM1" %in% colnames(celllines)){
	  MSK_data <- read.table("../data/features_SC2/similarity_ordered_prepared_data_MSK.csv")
       }else{
	  MSK_data <- read.table("../data/features_SC2/similarity_ordered_prepared_data_MSK_reduced.csv")
       }
       MSK_train <- MSK_data[,colnames(celllines)]
       MSK_response <- MSK_data[,names(MSK_data) %in%c("OS","OS_status")]
       MSK_response$time <- MSK_response$OS
       MSK_response$status <- MSK_response$OS_status
       MSK_response$OS <- NULL
       MSK_response$OS_status <- NULL
    }
    # Parallel CV
    print("start CV")
    rescv = mclapply(seq(nexp),
                   FUN=function(iexp) {
                       cat(".")
                       celllinesTrain = celllines[-folds[[iexp]],,drop=F]
                       celllinesTest = celllines[folds[[iexp]],,drop=F]
                       responseTrain = response[-folds[[iexp]],]
                       responseTest = response[folds[[iexp]],]                     
                       # Train and Predcit
                       if(use_MSK>0){
                           num_rows <- round(nrow(MSK_data)*use_MSK/100)  # num of row to use based on use_MSK percentage 
                           celllinesTrain <- rbind(celllinesTrain,MSK_train[1:num_rows,])
                           responseTrain <- rbind(responseTrain,MSK_response[1:num_rows,])
                       }

                       predict.test = mypredictor(celllinesTrain=celllinesTrain, celllinesTest=celllinesTest, response=responseTrain, ...)
                       
                       # Evaluate CI on the test
                       ci.test = suppressWarnings(survConcordance(Surv(time,status) ~ predict.test, as.data.frame(responseTest)))
                       return(as.vector(ci.test$concordance))
                   },
                   mc.cores=mc.cores )
                   
    return(unlist(rescv))

}




