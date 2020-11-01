# Script to evaluate the performance of a prediction method.

# OUTPUT
# a list with the performance of the method

setwd("./")
source("myutils.R")
source("myevaluate_SC1.R")
source("predictors_SC1.R")

library(parallel) 
library(optparse)

option_list <- list(
                    make_option("--predictor", action="store_true", default=FALSE, type="character", help="Predictor to use"),
                    make_option("--drug_responses", action="store_true", default="all", type="character", help="Drugs to valuate [default %default]"),
                    make_option("--index_train", action="store_true", default=FALSE, type="character", help="Patient indexes for the training set"),
                    make_option("--index_test", action="store_true", default=FALSE, type="character", help="Patient indexes for the test set"),
                    make_option("--times", action="store_true", default=FALSE, type="integer", help="Index of folds in cross validation"),
                    make_option("--nfolds", action="store_true", default=5, type="integer", help="Number of folds in CV, [default %default]"),
                    make_option("--nrepeats", type="integer", default=10, help="Number of repeats in CV [default %default]"),
                    make_option("--seed", type="integer", default=9182456, help="Seed for reproducibility [default %default]"),
                    make_option("--mc.cores", type="integer", default=4, help="Number of cores for parallelization [default %default]"),
                    make_option("--feature_model", type="character", default=FALSE, help="Vector of features"),
                    make_option("--alpha", type="numeric", default=0, help="alpha for elastic net [default %default]"),
                    make_option("--ntree", type="integer", default=50, help="Number of trees in RF [default %default]"),
                    make_option("--nodesize", type="integer", default=20, help="Node size in RF [default %default]"),
                    make_option("--max_depth", type="integer", default=2, help="Maximum depth in XGBoost [default %default]"),
                    make_option("--eta", type="integer", default=1, help="Eta in XGBoost [default %default]"),
                    make_option("--nthread", type="integer", default=2, help="Number of threads in XGBoost [default %default]"),
                    make_option("--nrounds", type="integer", default=20, help="Number of rounds in XGBoost [default %default]"),
                    make_option("--k", type="integer", default=5, help="Number of neighbors in kNN [default %default]"),
                    make_option("--hidden", type="character", default="c(50,20)", help="Number of nodes per hidden layer in NN [default %default]"),
                    make_option("--algorithm", type="character", default="backprop", help="Learning Algorithm in NN [default %default]"),
                    make_option("--learning_rate", type="numeric", default=0.001, help="Learning rate in NN [default %default]"),
                    make_option("--activation_function", type="character", default="tanh", help="Activation function in NN [default %default]"),
                    make_option("--kernel", type="character", default="radial", help="Kernel in SVR [default %default]"),
                    make_option("--filename", type="character", default="FALSE", help="Filename to save the output"),
                    make_option("--outpath", type="character", default="FALSE", help="Path to write the output")
                )

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser)

mypredictor <- eval(parse(text=args$predictor))
drug_responses <- args$drug_responses
index_train <- args$index_train
index_test <- args$index_test
time <- args$times
nfolds <- args$nfolds
nrepeats <- args$nrepeats
seed <- args$seed
mc.cores <- args$mc.cores
combination <- args$feature_model
alpha <- args$alpha
ntree <- args$ntree
nodesize <- args$nodesize
max_depth <- args$max_depth
eta <- args$eta
nthread <- args$nthread
nrounds <- args$nrounds
k <- args$k
hidden <- eval(parse(text=args$hidden))
algorithm <- args$algorithm
learning_rate <- args$learning_rate
activation_function <- args$activation_function
kernel <- args$kernel
filename <- args$filename
outpath <- args$outpath

i_train <- as.numeric(strsplit(index_train,"_")[[1]])
i_test <- as.numeric(strsplit(index_test,"_")[[1]])

feature_model <- strsplit(combination,"with")[[1]]
celllines <- combine_features(feature_model)

if(drug_responses=="all"){
    
    response <- read.table("../data/responses/prepared_data_auc.csv")

}else{
        
    response <- read.table("../data/responses/prepared_data_auc.csv")[,drug_responses]
}


if(identical(predictorElasticNet,mypredictor)){
    
    result <- evaluateCV(mypredictor, celllines, response, i_train, alpha=alpha, nrepeats=nrepeats, nfolds=nfolds, mc.cores=mc.cores)

}else if (identical(predictorRF,mypredictor)){
    
    result <- evaluateCV(mypredictor, celllines, response, i_train, nntree=ntree, nodesize=nodesize, nrepeats=nrepeats, nfolds=nfolds, mc.cores=mc.cores)

}else if (identical(predictorXGBoost,mypredictor)){
    
    result <- evaluateCV(mypredictor, celllines, response, i_train, max.depth=max_depth, eta=eta, nthread=nthread, nrounds=nrounds, nrepeats=nrepeats, nfolds=nfolds, mc.cores=mc.cores)
    
}else if (identical(predictorkNNReg,mypredictor)){
    
    result <- evaluateCV(mypredictor, celllines, response, i_train, k=k, nrepeats=nrepeats, nfolds=nfolds, mc.cores=mc.cores)
    
}else if (identical(predictorSVR,mypredictor)){
    
    result <- evaluateCV(mypredictor, celllines, response, i_train, nrepeats=nrepeats, nfolds=nfolds, mc.cores=mc.cores, kernel=kernel)
}

df_tosave <- do.call(rbind,result)
chemicals <- colnames(response)
rownames(df_tosave)<-chemicals

filename <- paste(filename,"_drugs_",paste(drug_responses,collapse="_"),".csv",sep="")
print(paste(outpath,"/",filename,sep=""))
write.table(df_tosave,paste(outpath,"/",filename,sep=""),row.names=T,col.names=T,sep="\t")


