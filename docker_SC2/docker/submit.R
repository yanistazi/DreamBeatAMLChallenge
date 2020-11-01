# Script to execute inside docker image

setwd("/docker")

source("/docker/docker_data_prep.R")
source("/docker/predict.R")

create_features()
list.files("/input")

#predict_from_best_model(rds_path="final_model/", path="/features/")

#predict_from_best_model(rds_path="final_model_challenge/", path="/features/")
#predict_from_best_model(rds_path="final_model_MSK_challenge/", path="/features/")

weight_challenge=0.50
predict_ensemble_learning(rds_challenge_path="final_model_challenge/",
			  rds_MSKplusChallenge_path="final_model_MSK_challenge/",
			  path.predict="/features/",
			  weight_challenge=weight_challenge)

