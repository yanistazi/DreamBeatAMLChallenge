# Data preparation R file for real time processing of the data.

# Inside the Docker image there will be:
# 1) an input folder where raw data is going to be [as per challenge organizers]
# 2) a features folder where the processed features will be written by our script
# 3) a fordocker folder where we will put|copy all the needed files/info from the training data
# 4) an output folder where the predictions.csv will be written by our script

library(reshape2)

create_kernels <- function(
			   featurepath="/features" ,
			   kernelpath="/kernels" ,
			   training_X_gene_restrict_path="/docker/forkernels/training_X_gene_restrict.txt" ,
			   training_Kscale_gene_restrict_path="/docker/forkernels/training_Kscale_gene_restrict.txt"
			   ) {


   system(paste0("mkdir -p ",kernelpath)) # does not hurt


   #----------------- RNA-seq RBF Kernel -------------------------#
   # ---------------- using a "restrict" ~5000 genes -------------#
   xtrain = read.table(training_X_gene_restrict_path, sep="\t", header=T, stringsAsFactors=F)
   xtest = read.table(paste0(featurepath,"/prepared_data_rna_full.csv"), sep="\t", header=T, stringsAsFactors=F)
   xtest = xtest[,colnames(xtrain)]

   ddscale = read.table(training_Kscale_gene_restrict_path, header=F, sep="\t")$V1

   xx = rbind(xtrain, xtest)
   mytrain = 1:nrow(xtrain)
   mytest = (1+nrow(xtrain)):(nrow(xtest)+nrow(xtrain))
   d <- dist(xx)
   dm <- as.matrix(d)
   dd <- dm^2
   dd <- dd/ddscale

   ddtest <- dd[mytest, mytrain]

   siglist <- c(1/20 , 1/10 , 1/5 , 2/5 , 3/5 , 4/5 , 1, 2, 4, 8)

   list.rbf <- list()
   for (i in seq(length(siglist))) {
      
      k <- exp(-ddtest/siglist[i])
      list.rbf[[i]] = k

      paste0(kernelpath,'/KrnarestrictRbf',i,'.test')

      write.table(k,paste0(kernelpath,'/KrnarestrictRbf',i,'.txt'),quote=FALSE,sep='\t')
   }

   # Mean rnaseq kernel
   KrnaseqMean <- list.rbf[[1]]
   for(i in 2:length(list.rbf)) {
      KrnaseqMean <- KrnaseqMean + list.rbf[[i]]
   }
   KrnaseqMean <- KrnaseqMean/length(list.rbf)
   write.table(KrnaseqMean , paste0(kernelpath,'/KrnarestrictMean','.txt') ,quote=FALSE,sep='\t')

   KrnaseqMean510 <- list.rbf[[5]]
   for(i in 6:length(list.rbf)) {
      KrnaseqMean510 <- KrnaseqMean510 + list.rbf[[i]]
   }
   KrnaseqMean510 <- KrnaseqMean510/length(list.rbf[5:length(list.rbf)])
   write.table(KrnaseqMean510 , paste0(kernelpath,'/KrnarestrictMean510','.txt') ,quote=FALSE,sep='\t')

   return()
}
