# Script to execute inside docker image

setwd("/docker")

source("/docker/docker_data_prep.R")
source("/docker/docker_kernels_prep.R")
source("/docker/predict.R")

create_features()
create_kernels()
list.files("/input")
list.files("/kernels")

#metric <- "" # if the spearman metric needs to be used, please use metric <- ""
#data_matrix <- read.table(paste("best_models","",".csv",sep=""),sep="\t")
#predict_from_best_model(data_matrix, metric = metric)

drug_file = "/docker/fordocker/drug_names_matched.tsv"
drug_df = read.table(drug_file, sep="\t", stringsAsFactors=F, header=T)
drug_switch_file="drug_switch_kmr_rf.txt"
predict_from_best_model_kmr_switch(drug_df, pathkernel = "/kernels/" , namekernel="KrnarestrictMean.txt",drug_switch_file=drug_switch_file)

