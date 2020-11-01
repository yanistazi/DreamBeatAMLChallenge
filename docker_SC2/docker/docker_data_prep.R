# Data preparation R file for real time processing of the data.

# Inside the Docker image there will be:
# 1) an input folder where raw data is going to be [as per challenge organizers]
# 2) a features folder where the processed features will be written by our script
# 3) a fordocker folder where we will put|copy all the needed files/info from the training data
# 4) an output folder where the predictions.csv will be written by our script

library(reshape2)

create_features <- function(trainingPcaPath="/docker/fordocker/trainingPCA.RData",
			    trainingMedianPath="/docker/fordocker/training_median.tsv",
			    trainingGenePath="/docker/fordocker/training_gene_mutation.csv",
# 			    vargenesPath="/docker/fordocker/vargenes.csv",
# 			    vargenes10Path="/docker/fordocker/vargenes10.csv",
# 			    vargenes20Path="/docker/fordocker/vargenes20.csv",
# 			    vargenes30Path="/docker/fordocker/vargenes30.csv",
# 			    hrgenesPath="/docker/fordocker/hrgenes.csv",
# 			    hrgenes10Path="/docker/fordocker/hrgenes10.csv",
# 			    hrgenes20Path="/docker/fordocker/hrgenes20.csv",
# 			    hrgenes30Path="/docker/fordocker/hrgenes30.csv",
# 			    hrvargenesPath="/docker/fordocker/hrvargenes.csv",
			    mskfeaturesPath="/docker/fordocker/names_combined_features.tsv",
			    mskfeaturesreducedPath="/docker/fordocker/names_combined_features_reduced.tsv",
			    inputpath="/input",
			    featurepath="/features"
			    ) {

   
   #----------------- Read Input -------------------------#
   trainingMedian = read.table(trainingMedianPath, sep="\t", stringsAsFactors=F, header=T)
   trainingGene = read.table(trainingGenePath, sep="\t", stringsAsFactors=F, header=T)

   #----------------- clinical_numerical -----------------#
   df_clin_num <- read.csv(paste0(inputpath,"/clinical_numerical.csv"))
   rownames(df_clin_num) <- df_clin_num$lab_id
   df_clin_num$lab_id <- NULL
   df_clin_num_imputed <- df_clin_num
   df_clin_num_imputed <- df_clin_num_imputed[,-which(colnames(df_clin_num_imputed)=="ageAtDiagnosis")]
   # Fill NA if needed:
   for (cc in colnames(df_clin_num_imputed)){
      j = which(is.na(df_clin_num_imputed[,cc]))
      if (length(j)>0) {
	 df_clin_num_imputed[j,cc] = trainingMedian[trainingMedian$name==cc,"value"] 
      }
   }
   write.table(df_clin_num_imputed,paste0(featurepath,"/prepared_data_clin_num_imputed.csv"),row.names=T,col.names=T,sep="\t")


   #----------------- AUC -----------------#
   df_auc <- read.csv(paste0(inputpath,"/aucs.csv"))
   data_auc <- (dcast(df_auc, lab_id ~ inhibitor,value.var="auc"))
   rownames(data_auc) <- data_auc$lab_id
   jm = match(rownames(df_clin_num),data_auc$lab_id)
   data_auc = data_auc[jm,]
   data_auc$lab_id <- NULL
   data_auc_imputed <- apply(data_auc,2,as.numeric)
   rownames(data_auc_imputed) <- rownames(data_auc)
   # Fill NA if needed:
   for (cc in colnames(data_auc_imputed)){
      j = which(is.na(data_auc_imputed[,cc]))
      if (length(j)>0) {
	 data_auc_imputed[j,cc] = trainingMedian[trainingMedian$name==cc,"value"] 
      }
   }
   write.table(data_auc_imputed,paste0(featurepath,"/prepared_data_auc_imputed.csv"),row.names=T,col.names=T,sep="\t")


   # MOL = INTEGRATION clinical_categorical & DNA mutation
   #----------------- clinical_categorical -----------------#
   df_clin_cat <- read.csv(paste0(inputpath,"/clinical_categorical.csv"))
   rownames(df_clin_cat) <- df_clin_cat$lab_id
   jm = match(rownames(df_clin_num),df_clin_cat$lab_id)
   df_clin_cat = df_clin_cat[jm,]
   df_clin_cat$lab_id <- NULL
   df_clin_cat$priorMalignancyMyeloid = as.numeric(df_clin_cat$priorMDS | 
						   df_clin_cat$priorMPN |
						   df_clin_cat$priorMDSMPN)
   keep = c(
	    "priorMalignancyNonMyeloid",
	    "priorMalignancyMyeloid",
	    "specificDxAtAcquisition",
	    "consensus_sex",
	    "FLT3.ITD"
   )
   df_clin_cat = df_clin_cat[,keep]
   df_clin_cat$inv_16 = 0 
   df_clin_cat$inv_16[which(df_clin_cat$specificDxAtAcquisition==0)] = 1
   df_clin_cat$inv_3 = 0 
   df_clin_cat$inv_3[which(df_clin_cat$specificDxAtAcquisition==1)] = 1
   df_clin_cat$CEBPAbi = 0 
   df_clin_cat$CEBPAbi[which(df_clin_cat$specificDxAtAcquisition==3)] = 1
   df_clin_cat$NPM1 = 0 
   df_clin_cat$NPM1[which(df_clin_cat$specificDxAtAcquisition==4)] = 1
   df_clin_cat$t_8_21 = 0  # consistent with the fusion field 
   df_clin_cat$t_8_21[which(df_clin_cat$specificDxAtAcquisition==6)] = 1
   df_clin_cat$t_9_11 = 0  # consistent with the fusion field 
   df_clin_cat$t_9_11[which(df_clin_cat$specificDxAtAcquisition==7)] = 1
   df_clin_cat$t_15_17 = 0  # consistent with the fusion field 
   df_clin_cat$t_15_17[which(df_clin_cat$specificDxAtAcquisition==11)] = 1
   save_clin_cat = df_clin_cat
   df_clin_cat$AML_therapy = 0
   df_clin_cat$AML_therapy[which(df_clin_cat$specificDxAtAcquisition==12)] = 1
   df_clin_cat$AML_therapy_secondary = 0
   df_clin_cat$AML_therapy_secondary[which(df_clin_cat$priorMalignancyMyeloid==1 | df_clin_cat$AML_therapy==1)] = 1
   df_clin_cat$AML_therapy = NULL
   df_clin_cat$inv_3 = NULL
   df_clin_cat$priorMalignancyNonMyeloid = NULL
   df_clin_cat$priorMalignancyMyeloid = NULL
   df_clin_cat$specificDxAtAcquisition = NULL
   df_clin_cat$specificDxAtAcquisition = NULL
   df_clin_cat$ITD = df_clin_cat$FLT3.ITD
   df_clin_cat$FLT3.ITD = NULL
   #----------------- DNAseq -----------------#
   # Here we need to fill a dataframe npatient_test X ngene_train
   df_dna <- read.csv(paste0(inputpath,"/dnaseq.csv"))
   data_dna = matrix(0, nrow=nrow(df_clin_num), ncol=length(trainingGene$gene))
   data_dna = as.data.frame(data_dna)
   colnames(data_dna) = trainingGene$gene
   rownames(data_dna) = rownames(df_clin_num)
   for (gg in trainingGene$gene) {
      j = which(df_dna$Hugo_Symbol==gg)
      if (length(j)>0) {
	 gopatient = as.vector(unique(df_dna$lab_id[j]))
	 data_dna[rownames(data_dna)%in%gopatient,gg] = 1
      }
   }
   # Merge dataframe:
   df_clin_mut = df_clin_cat
   df_clin_mut$NPM1[which(data_dna$NPM1==1)] = 1
   data_dna$NPM1 = NULL # NPM1 redundant
   df_clin_mut = cbind(df_clin_mut, data_dna)
   # Fill NA if needed:
   for (cc in colnames(df_clin_mut)){
      j = which(is.na(df_clin_mut[,cc]))
      if (length(j)>0) {
	 df_clin_mut[j,cc] = trainingMedian[trainingMedian$name==cc,"value"] 
      }
   }
   write.table(df_clin_mut,paste0(featurepath,"/prepared_data_clin_cat_mut_clean.csv"),row.names=T,col.names=T,sep="\t")


   #----------------- MSK DATA -----------------#
   msk.features = read.table(mskfeaturesPath, sep="\t", stringsAsFactors=F, header=F)$V1
   msk.genes = msk.features[1:21]
   tmp_dna = matrix(0, nrow=nrow(df_clin_num), ncol=length(msk.genes))
   tmp_dna = as.data.frame(tmp_dna)
   colnames(tmp_dna) = msk.genes 
   rownames(tmp_dna) = rownames(df_clin_num)
   for (gg in msk.genes) {
      j = which(df_dna$Hugo_Symbol==gg)
      if (length(j)>0) {
	 gopatient = as.vector(unique(df_dna$lab_id[j]))
	 tmp_dna[rownames(tmp_dna)%in%gopatient,gg] = 1
      }
   }
   tmp_dna$NPM1[which(df_clin_mut$NPM1==1)]=1
   tmp2 = df_clin_cat[,msk.features[22:23]]
   tmp2$inv_3 = 0
   tmp2$inv_3[which(save_clin_cat$inv_3==1)] = 1
   tmp2$CEBPA_bi = df_clin_mut$CEBPAbi
   tmp3 = df_clin_mut[msk.features[26:28]]
   tmp3$gender = df_clin_cat$consensus_sex
   tmp3$bm_blasts = df_clin_num_imputed$X..Blasts.in.PB
   tmp3$wbc = df_clin_num_imputed$WBC.Count
   tmp3$age = df_clin_num_imputed$ageAtSpecimenAcquisition
   df_for_combined = cbind(tmp_dna,tmp2,tmp3)
   write.table(df_for_combined,paste0(featurepath,"/prepared_data_challenge_combined.csv"),row.names=T,col.names=T,sep="\t")

   print("NA?")
   print(apply(df_for_combined,2,function(x) sum(is.na(x))))

   msk.features.reduced = read.table(mskfeaturesreducedPath, sep="\t", stringsAsFactors=F, header=F)$V1
   df_for_combined_reduced = df_for_combined[,msk.features.reduced]
   write.table(df_for_combined_reduced,paste0(featurepath,"/prepared_data_challenge_combined_reduced.csv"),row.names=T,col.names=T,sep="\t")
   

#    #----------------- RNAseq -----------------#
#    # Here we project the test gene expression data into the training set PC
#    df_rna <- read.csv(paste0(inputpath,"/rnaseq.csv"))
#    df_rna <- t(df_rna)
#    # Use genes as column names
#    colnames(df_rna) <- df_rna[2,]
#    data_rna <- data.frame(df_rna[3:nrow(df_rna),])
#    # Change row name to match with patient names in other datasets
#    rownames(data_rna) <- gsub("X","",rownames(data_rna))
#    rownames(data_rna) <- gsub("\\.","-",rownames(data_rna))
#    # Make data.frame values numeric 
#    indx <- sapply(data_rna, is.factor)
#    data_rna[indx] <- lapply(data_rna[indx], function(x) as.numeric(as.character(x)))
#    save_data_rna <- data_rna
#    # Load training PCA object
#    load(file=trainingPcaPath) #---> x.pca object
#    gene.pca = rownames(x.pca$rotation)
#    data_rna = data_rna[,gene.pca]
#    reduced_rna <- predict(x.pca, newdata=data_rna)
#    n_components <- 71 
#    reduced_rna <- reduced_rna[,1:n_components]
#    # identical(rownames(df_clin_num),rownames(reduced_rna)) # Checked True
#    write.table(reduced_rna,paste0(featurepath,"/prepared_data_rna_reduced.csv"),row.names=T,col.names=T,sep="\t")

#    #----------------- RNAseq -----------------#
#    # Here we select the genes with the highest variance in gene expression based on the training set 
#    vargenes <- read.csv(vargenesPath,header=F)                            
#    data_rna_vargenes <- save_data_rna[,as.character(vargenes$V1)]
#    write.table(data_rna_vargenes,paste0(featurepath,"/prepared_data_rna_vargenes.csv"),row.names=T,col.names=T,sep="\t")
#    # 10 20 30
#    vargenes10 <- read.csv(vargenes10Path,header=F)                            
#    data_rna_vargenes10 <- save_data_rna[,as.character(vargenes10$V1)]
#    write.table(data_rna_vargenes10,paste0(featurepath,"/prepared_data_rna_vargenes10.csv"),row.names=T,col.names=T,sep="\t")

#    vargenes20 <- read.csv(vargenes20Path,header=F)                            
#    data_rna_vargenes20 <- save_data_rna[,as.character(vargenes20$V1)]
#    write.table(data_rna_vargenes20,paste0(featurepath,"/prepared_data_rna_vargenes20.csv"),row.names=T,col.names=T,sep="\t")

#    vargenes30 <- read.csv(vargenes30Path,header=F)                            
#    data_rna_vargenes30 <- save_data_rna[,as.character(vargenes30$V1)]
#    write.table(data_rna_vargenes30,paste0(featurepath,"/prepared_data_rna_vargenes30.csv"),row.names=T,col.names=T,sep="\t")

#    #----------------- RNAseq -----------------#
#    # Here we select the genes with most significant hazard ratio
#    hrgenes <- read.csv(hrgenesPath,header=F)                            
#    data_rna_hrgenes <- save_data_rna[,as.character(hrgenes$V1)]
#    write.table(data_rna_hrgenes,paste0(featurepath,"/prepared_data_rna_hrgenes.csv"),row.names=T,col.names=T,sep="\t")
#    # 10 20 30
#    hrgenes10 <- read.csv(hrgenes10Path,header=F)                            
#    data_rna_hrgenes10 <- save_data_rna[,as.character(hrgenes10$V1)]
#    write.table(data_rna_hrgenes10,paste0(featurepath,"/prepared_data_rna_hrgenes10.csv"),row.names=T,col.names=T,sep="\t")

#    hrgenes20 <- read.csv(hrgenes20Path,header=F)                            
#    data_rna_hrgenes20 <- save_data_rna[,as.character(hrgenes20$V1)]
#    write.table(data_rna_hrgenes20,paste0(featurepath,"/prepared_data_rna_hrgenes20.csv"),row.names=T,col.names=T,sep="\t")

#    hrgenes30 <- read.csv(hrgenes30Path,header=F)                            
#    data_rna_hrgenes30 <- save_data_rna[,as.character(hrgenes30$V1)]
#    write.table(data_rna_hrgenes30,paste0(featurepath,"/prepared_data_rna_hrgenes30.csv"),row.names=T,col.names=T,sep="\t")

#    # HR genes that are also VAR genes (n=41...)
#    hrvargenes <- read.csv(hrvargenesPath,header=F)                            
#    data_rna_hrvargenes <- save_data_rna[,as.character(hrvargenes$V1)]
#    write.table(data_rna_hrvargenes,paste0(featurepath,"/prepared_data_rna_hrvargenes.csv"),row.names=T,col.names=T,sep="\t")

   return()
}
