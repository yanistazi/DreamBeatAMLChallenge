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


#1 ) GLM -GLM
mypredictor_challenge <- predictorGLM
mypredictor_challenge_MSK <- predictorGLM

#alpha_challenge <- c(0,0.25,0.5,0.75,1)
#alpha_challenge_MSK <- c(0,0.25,0.5,0.75,1)
#weight_challenge <- c(0,0.25,0.5,0.75,1)
#use_MSK <- c(25,50,75,100)

alpha_challenge <- c(0,0.25,0.5,0.75)
alpha_challenge_MSK <- c(0,0.25,0.5,0.75)
weight_challenge <- c(0.25,0.5,0.75)
use_MSK <- c(100)

df <- data.frame()

for (a in alpha_challenge){
   for (b in alpha_challenge_MSK){
      for (c in weight_challenge) {
	 for (d in use_MSK) {
	    res <- evaluateCV_SC2_ensemble_transfer_learning(mypredictor_challenge,mypredictor_challenge_MSK,
								 celllines_challenge, celllines_challenge_MSK, response,
								 MSK.path = "../../data/features_SC2/",
								 nfolds=5, nrepeats=5, seed=17,,weight_challenge=c,
								 mc.cores=mc.cores,use_MSK=d,alpha_challenge=a,alpha_challenge_MSK=b)

	    df[paste("GLM_GLM",a,b,c,d,sep="_"),c("mean_ci","sd_ci","min_ci")] <- c(mean(res),sd(res),min(res))

	 }
      }
   }
}
write.table(df,"GLM_GLM_reduced.tsv", sep="\t", quote=F)
