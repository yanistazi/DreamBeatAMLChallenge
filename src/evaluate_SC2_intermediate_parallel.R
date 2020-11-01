# Script to evaluate the performance of a prediction method.

# OUTPUT
# a list with the performance of the method

setwd("./")
source("utils_SC2.R")
source("evaluate_SC2.R")
source("predictors_SC2.R")

library(parallel) 
library(optparse)

option_list <- list(
                    make_option("--predictor", action="store_true", default=FALSE, type="character", help="Predictor to use"),
                    make_option("--drug_responses", action="store_true", default="all", type="character", help="Drugs to valuate [default %default]"),
                    make_option("--nfolds", action="store_true", default=5, type="integer", help="Number of folds in CV, [default %default]"),
                    make_option("--use_MSK", type="integer", help="Percentage of MSK Data to use, [default %default]"),
                    make_option("--nrepeats", type="integer", default=10, help="Number of repeats in CV [default %default]"),
                    make_option("--seed", type="integer", default=17, help="Seed for reproducibility [default %default]"),
                    make_option("--mc.cores", type="integer", default=4, help="Number of cores for parallelization [default %default]"),
                    make_option("--feature_model", type="character", default=FALSE, help="Vector of features"),
                    make_option("--alpha", type="numeric", default=0, help="alpha for elastic net [default %default]"),
                    make_option("--ntree", type="integer", default=50, help="Number of trees in RF [default %default]"),
                    make_option("--nodesize", type="integer", default=20, help="Node size in RF [default %default]"),
                    make_option("--maxstepno", type="integer", default=500, help="Maximum stepth in CoxBoost CV [default %default]"),
                    make_option("--type", type="character", default="verweij", help="Type in Cox Boost [default %default"),
                    make_option("--penalty", type="integer", default=100, help="Penalty in Cox Boost [default %default]"),
                    make_option("--K", type="integer", default=10, help="K in Cox Boost [default %default]"),
                    make_option("--iter", type="integer", default=500, help="Max number of iterations in Random Effects [default %default]"),
                    make_option("--tol", type="numeric", default=0.01, help="tolerance in Random Effects [default %default]"),
                    make_option("--nthread", type="integer", default=2, help="Number of threads in XGBoost [default %default]"),
                    make_option("--nrounds", type="integer", default=20, help="Number of rounds in XGBoost [default %default]"),
                    make_option("--filename", type="character", default="FALSE", help="Filename to save the output"),
                    make_option("--outpath", type="character", default="FALSE", help="Path to write the output")
                )

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser)

mypredictor <- eval(parse(text=args$predictor))
nfolds <- args$nfolds
nrepeats <- args$nrepeats
seed <- args$seed
mc.cores <- args$mc.cores
combination <- args$feature_model
alpha <- args$alpha
ntree <- args$ntree
nodesize <- args$nodesize
maxstepno <- args$maxstepno
K <- args$K
type <- args$type
penalty <- args$penalty
iter <- args$iter
tol <- args$tol
nthread <- args$nthread
nrounds <- args$nrounds
filename <- args$filename
outpath <- args$outpath
use_MSK <- args$use_MSK

feature_model <- strsplit(combination,"with")[[1]]
celllines <- combine_features(feature_model)
print(dim(celllines))

    
response <- read.table("../data/responses/prepared_data_response.csv")

if(identical(predictorGLM,mypredictor)){
    
    result <- evaluateCV_SC2(mypredictor, celllines, response, alpha=alpha, nrepeats=nrepeats, nfolds=nfolds, mc.cores=mc.cores,use_MSK=use_MSK)

}else if (identical(predictorRF,mypredictor)){
    
    result <- evaluateCV_SC2(mypredictor, celllines, response, nntree=ntree, nodesize=nodesize, nrepeats=nrepeats, nfolds=nfolds, mc.cores=mc.cores,use_MSK=use_MSK)

}else if (identical(predictorBoost,mypredictor)){
    
    result <- evaluateCV_SC2(mypredictor, celllines, response, maxstepno=maxstepno , K=K , type=type , penalty=penalty , nrepeats=nrepeats, nfolds=nfolds, mc.cores=mc.cores,use_MSK=use_MSK)
    
}else if (identical(predictorRFX,mypredictor)){
    
    result <- evaluateCV_SC2(mypredictor, celllines, response, max.iter = iter ,tol=tol, nrepeats=nrepeats, nfolds=nfolds, mc.cores=mc.cores,use_MSK=use_MSK)
    
}

df_tosave <- data.frame(mean=mean(result),sd=sd(result),min=min(result),max=max(result))

filename <- paste(filename,".csv",sep="")
print(paste(outpath,"/",filename,sep=""))
write.table(df_tosave,paste(outpath,"/",filename,sep=""),row.names=T,col.names=T,sep="\t")


