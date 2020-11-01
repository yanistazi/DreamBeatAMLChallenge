# Function to identify drugs where the RF predictions
# are consistently better than the KMR predictions
# Consistent: across ranges of kernels | parameters | ie not outlier better

output_drugs_kmr_switch <- function(
				    path.kmr.cv, 
				    path.rf.cv,
				    keep.kernels = c(paste0("KrnarestrictRbf",6:9),"KrnarestrictMean"),
				    cutoff = 0.05,
				    drug_switch_file = "drug_switch_kmr_rf.txt",
				    drug_switch_path = "./../docker/Kdocker/docker"
				    ) {

   # KMR CV RESULTS
   myfiles = list.files(path.kmr.cv)
   lres = lapply(myfiles, function(x) read.table(paste(path.kmr.cv,x,sep="/"),sep="\t",header=T,stringsAsFactors=F))
   dres = do.call(rbind, lres)
   dgood = dres[dres$Kpatient%in%keep.kernels,]
   # Mean perf per drug across all decent kernel trials
   kkl = lapply(unique(as.vector(dgood$drug)), function(dd) { 
		   x=dgood[dgood$drug==dd,]
		   y = data.frame(meanRMSD=mean(x$meanRMSD),sdRMSD=sd(x$sdRMSD),
				  meanSpearman=mean(x$meanSpearman),sdmeanSpearman=sd(x$meanSpearman),drug=dd) } )
   dd.kmr = do.call(rbind, kkl)


   # RF CV RESULTS
   rfiles = list.files(path.rf.cv)
   r = rfiles[grepl(paste(c("predictorRF","predictorkNNReg"),collapse="|"),rfiles)]
   rl = lapply(r, function(x) {
		  y=read.table(paste(path.rf.cv,x,sep="/"),sep="\t",header=T,stringsAsFactors=F) 
		  y$drug = rownames(y)
		  return(y)
		   })
   drl = do.call(rbind, rl)
   # Mean perf per drug across all parameters trials
   rrl = lapply(unique(drl$drug), function(dd) {
		   x=drl[drl$drug==dd,]
		   y = data.frame(meanRMSD=mean(x$meanRMSD),sdRMSD=sd(x$sdRMSD),
				  meanSpearman=mean(x$meanSpearman),sdmeanSpearman=sd(x$meanSpearman),drug=dd) } )
   dd.rf = do.call(rbind, rrl)

   # IDENTIFY "BAD" KMR DRUGS
   jbad = which((dd.kmr$meanSpearman+dd.kmr$sdmeanSpearman) < (dd.rf $meanSpearman-dd.rf $sdmeanSpearman-cutoff))
   baddrug = as.vector(dd.kmr$drug[jbad])

   write.table(data.frame(baddrug) , paste(drug_switch_path,drug_switch_file,sep="/") , col.names=F, row.names=F, sep="\t",quote=F)

   return()

}
