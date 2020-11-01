# Data preparation R file for real time processing of the data.

# Inside the Docker image there will be:
# 1) an input folder where raw data is going to be [as per challenge organizers]
# 2) a features folder where the processed features will be written by our script
# 3) a fordocker folder where we will put|copy all the needed files/info from the training data
# 4) an output folder where the predictions.csv will be written by our script

library(reshape2)

create_features <- function(
			    trainingMedianPath="/docker/fordocker/training_median_clinical.tsv",
			    trainingGenePath="/docker/fordocker/training_gene_mutation.tsv",
			    corgenesPath="/docker/fordocker/corgenes.csv", 
			    corgenesPath_enriched="/docker/fordocker/corgenes_enriched.csv",
			    inputpath="/input",
			    featurepath="/features"
			    ) {
   # Inputs:
   #  - trainingPcaPath: path for the trainingPCA RData object
   #  - trainingMedianPath: path for the files with median training set values for clinical variables
   #  - trainingGenePath: path for the gene used in training for mutation features
   # Output:
   # - write the feature dataframes in features/

   #system("mkdir -p features") # does not hurt

   #----------------- Read Input -------------------------#
   trainingMedian = read.table(trainingMedianPath, sep="\t", stringsAsFactors=F, header=T)
   trainingGene = read.table(trainingGenePath, sep="\t", stringsAsFactors=F, header=T)

   #----------------- clinical_numerical -----------------#
   df_clin_num <- read.csv(paste0(inputpath,"/clinical_numerical.csv"))
   rownames(df_clin_num) <- df_clin_num$lab_id
   df_clin_num$lab_id <- NULL
   # Fill NA if needed:
   for (cc in colnames(df_clin_num)){
      j = which(is.na(df_clin_num[,cc]))
      if (length(j)>0) {
	 df_clin_num[j,cc] = trainingMedian[trainingMedian$name==cc,"value"] 
      }
   }
   write.table(df_clin_num,paste0(featurepath,"/prepared_data_clin_num.csv"),row.names=T,col.names=T,sep="\t")


   #----------------- clinical_numerical_reduced -----------------#
   write.table(df_clin_num[,c(2,5)],paste0(featurepath,"/prepared_data_clin_num_reduced.csv"),row.names=T,col.names=T,sep="\t")


   #----------------- clinical_categorical -----------------#
   df_clin_cat <- read.csv(paste0(inputpath,"/clinical_categorical.csv"))
   rownames(df_clin_cat) <- df_clin_cat$lab_id
   jm = match(rownames(df_clin_num),df_clin_cat$lab_id)
   df_clin_cat = df_clin_cat[jm,]
   # identical(rownames(df_clin_num),rownames(df_clin_cat)) # Checked True
   df_clin_cat$lab_id <- NULL
   # Fill NA if needed:
   for (cc in colnames(df_clin_cat)){
      j = which(is.na(df_clin_cat[,cc]))
      if (length(j)>0) {
	 df_clin_cat[j,cc] = trainingMedian[trainingMedian$name==cc,"value"] 
      }
   }
   write.table(df_clin_cat,paste0(featurepath,"/prepared_data_clin_cat.csv"),row.names=T,col.names=T,sep="\t")

   #----------------- clinical_categorical_reduced -----------------#
   # We can simplify the clinical categorical data a little bit
   # As is it overall known that "prior MDS | MPN | MDSMPN" are very agressive
   df_clin_cat$priorMalignancyMyeloid = as.numeric(df_clin_cat$priorMDS | 
						   df_clin_cat$priorMPN |
						   df_clin_cat$priorMDSMPN)
   # Note that "dxAtInclusion" and "dxAtSpecimenAcquisition" 
   # is redundant with "priorMDS|MPN|MDSMPN"                                     
   keep = c(
	    "priorMalignancyNonMyeloid",
	    "priorMalignancyType",
	    "priorMalignancyMyeloid",
	    "specificDxAtInclusion",
	    "specificDxAtAcquisition",
	    "specimenType",
	    "consensus_sex",
	    "FAB.Blast.Morphology",
	    "Karyotype",
	    "FLT3.ITD",
	    "finalFusion"
   )
   df_clin_cat_reduced = df_clin_cat[,keep]
   write.table(df_clin_cat_reduced,paste0(featurepath,"/prepared_data_clin_cat_reduced.csv"),row.names=T,col.names=T,sep="\t")


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
   write.table(data_dna,paste0(featurepath,"/prepared_data_dna.csv"),row.names=T,col.names=T,sep="\t")

   #----------------- RNAseq -----------------#
   # Here we project the test gene expression data into the training set PC
   df_rna <- read.csv(paste0(inputpath,"/rnaseq.csv"))
   df_rna <- t(df_rna)
   # Use genes as column names
   colnames(df_rna) <- df_rna[2,]
   data_rna <- data.frame(df_rna[3:nrow(df_rna),])
   # Change row name to match with patient names in other datasets
   rownames(data_rna) <- gsub("X","",rownames(data_rna))
   rownames(data_rna) <- gsub("\\.","-",rownames(data_rna))
   # Make data.frame values numeric 
   indx <- sapply(data_rna, is.factor)
   data_rna[indx] <- lapply(data_rna[indx], function(x) as.numeric(as.character(x)))
   write.table(data_rna,paste0(featurepath,"/prepared_data_rna_full.csv"),row.names=T,col.names=T,sep="\t")
   
   data_rna_corgenes <- data_rna
   data_rna_corgenes_enriched <- data_rna
   save_data_rna <- data_rna

   # Load training PCA object
   #load(file=trainingPcaPath) #---> x.pca object
   #gene.pca = rownames(x.pca$rotation)
   #data_rna = data_rna[,gene.pca]
   #reduced_rna <- predict(x.pca, newdata=data_rna)
   #n_components <- 71 
   #reduced_rna <- reduced_rna[,1:n_components]
   # identical(rownames(df_clin_num),rownames(reduced_rna)) # Checked True
   #save the reduced RNA dataset
   #write.table(reduced_rna,paste0(featurepath,"/prepared_data_rna_reduced.csv"),row.names=T,col.names=T,sep="\t")

   #----------------- RNAseq -----------------#
   # Here we select the genes with the highest correlation between gene expression and auc 
   corgenes <- read.csv(corgenesPath,header=F)                            
   data_rna_corgenes <- data_rna_corgenes[,as.character(corgenes$V1)]
   write.table(data_rna_corgenes,paste0(featurepath,"/prepared_data_rna_corgenes.csv"),row.names=T,col.names=T,sep="\t")

#    m <- apply(data_rna_corgenes,2,mean)
#    s <- apply(data_rna_corgenes, 2, sd)
#    data_rna_corgenes <- sweep(data_rna_corgenes,2,m)
#    data_rna_corgenes <- sweep(data_rna_corgenes, 2, s, "/")

#    write.table(data_rna_corgenes,paste0(featurepath,"/prepared_data_rna_corgenes_zscore.csv"),row.names=T,col.names=T,sep="\t")

   #----------------- RNAseq -----------------#
   # Here we select the genes with the highest correlation between gene expression and auc & enriching them with other genes too 
   corgenes_enriched <- read.csv(corgenesPath_enriched,header=F)                            
   data_rna_corgenes_enriched <- data_rna_corgenes_enriched[,as.character(corgenes_enriched$V1)]
   write.table(data_rna_corgenes_enriched,paste0(featurepath,"/prepared_data_rna_corgenes_enriched.csv"),row.names=T,col.names=T,sep="\t")

#    m <- apply(data_rna_corgenes_enriched,2,mean)
#    s <- apply(data_rna_corgenes_enriched, 2, sd)
#    data_rna_corgenes_enriched <- sweep(data_rna_corgenes_enriched,2,m)
#    data_rna_corgenes_enriched <- sweep(data_rna_corgenes_enriched, 2, s, "/")

#    write.table(data_rna_corgenes_enriched,paste0(featurepath,"/prepared_data_rna_corgenes_enriched_zscore.csv"),row.names=T,col.names=T,sep="\t")


   #----------------- RNAseq -----------------#
   # Here we select the genes with the highest correlation between gene expression and auc & genes with highest variance

#    data_var <- read.table(paste0(featurepath,"/prepared_data_rna_vargenes.csv"),header=T,sep="\t")
#    data_cor <- read.table(paste0(featurepath,"/prepared_data_rna_corgenes.csv"),header=T,sep="\t")

#    data_rna_var_cor <- save_data_rna[,union(colnames(data_var),colnames(data_cor))]
#    write.table(data_rna_var_cor,paste0(featurepath,"/prepared_data_rna_varcorgenes.csv"),row.names=T,col.names=T,sep="\t")

   return()
}
