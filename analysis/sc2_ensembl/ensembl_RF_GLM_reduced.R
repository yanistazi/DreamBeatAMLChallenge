# To launch from dream-beat/src

#.libPaths("/ifs/res/papaemme/users/eb2/dream-beat/repos")

source("../../src/evaluate_SC2.R")
source("../../src/predictors_SC2.R")
source("../../src/utils_SC2.R")
source("../../src/evaluateCV_SC2_ensemble_transfer_learning.R")

mc.cores = 10

celllines_challenge <- combine_features(c("clinical_numerical","mol"),path="../../data/features_SC2/")
celllines_challenge_MSK <- combine_features("combined_features_reduced",path="../../data/features_SC2/")
response <- read.table("../../data/responses/prepared_data_response.csv")


#2 ) RF -GLM
mypredictor_challenge <- predictorRF
mypredictor_challenge_MSK <- predictorGLM

ntree_challenge <- c(100,300,500)
nodesize_challenge <- c(5,10,25)
alpha_challenge_MSK <- c(0,0.25,0.75)
weight_challenge <- c(0.5,0.75)
use_MSK <- c(100)

df1 <- data.frame()

for (a in ntree_challenge){
   for (n in nodesize_challenge){
      for (b in alpha_challenge_MSK){
	 for (c in weight_challenge) {
	    for (d in use_MSK) {
	       res <- evaluateCV_SC2_ensemble_transfer_learning(mypredictor_challenge,mypredictor_challenge_MSK,
								    celllines_challenge, celllines_challenge_MSK, response,
								    MSK.path = "../../data/features_SC2/",
								    nfolds=5, nrepeats=5, seed=17,,weight_challenge=c,
								    mc.cores=mc.cores,use_MSK=d,ntree_challenge=a,
								    nodesize_challenge=n,alpha_challenge_MSK=b)

	       df1[paste("RF_GLM",a,n,b,c,d,sep="_"),c("mean_ci","sd_ci","min_ci")] <- c(mean(res),sd(res),min(res))

	    }
	 }
      }
   }
}
write.table(df1,"RF_GLM_reduced.tsv",sep="\t", quote=F)
