source("Kpredictors_SC1.R")
library(randomForest)
source("./../docker/Kdocker/docker/combine_features.R")
source("utils.R")

response <- read.table("./../data/responses/prepared_data_auc.csv")
data_matrix <- read.table("./../docker/Kdocker/docker/best_models_rmsd.csv")

# Drug Kernel
Kdrug = read.table("./../data/kernels/Kempirical.txt", sep="\t",header=T,stringsAsFactors=F)

# not used as using pre-computed kernels
# Celllines
#X = read.table("../../../data/features/prepared_data_rna_vargenes.csv", sep="\t",header=T)

# Patient Kernel
Kpatient = read.table("./../data/kernels/KrnarestrictMean.txt",header=T,stringsAsFactors=F)

### Save RDS File function
save_rds <- function(final_model,drug_name){
   
   saveRDS(final_model, paste("./../docker/Kdocker/docker/final_models/loaded_",sep=drug_name,".rds"))
}

drug_switch = read.table("./../docker/Kdocker/docker/drug_switch_kmr_rf.txt", header=F, sep="\t", stringsAsFactors=F)$V1 

system("mkdir -p ./../docker/Kdocker/docker/final_models")

for (i in 1:ncol(response)){

   print(paste0("Drug",i))

   if (!colnames(response)[i]%in%drug_switch) {

      print("KMR prediction")

      ### Handle Nas
      non_na <- which(!is.na(response[,i]))   
      #celllines_dr <- X[non_na,]
      response_dr <- response[non_na,,drop=F]
      Kpatient_dr = Kpatient[non_na,non_na]

      final_model = KtrainKMR(celllinesTrain = Kpatient_dr ,
			      responseTrain = response_dr ,
			      task.number = i ,
			      KernelTask = as.matrix(Kdrug) ,
			      kx_type = "precomputed", # linear gaussian
			      #kx_option = list(sigma = 1), # for gaussian kernel
			      type.measure = "ci", # ci mse cor
			      nfolds.intern = 3 , #7,
			      nrepeats.intern = 3 , #4,
			      lambdas = 10^(-10:10)#10^(-10:10)
      )

      final_model$non_na_index = non_na
      save_rds(final_model,colnames(response)[i])
       
   }else{
       
       print("no KMR prediction")
       save_trained_best_models(data.matrix(response[,rownames(data_matrix)[i],drop=F]),data_matrix[i,],path="../docker/Kdocker/docker/final_models/loaded_")
   }
}
