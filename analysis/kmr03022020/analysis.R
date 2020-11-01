library(ggplot2)
library(gplots)
library(ggpubr)
source("ggstyles.R")

myfiles = list.files("CVresults")
lres = lapply(myfiles, function(x) read.table(paste0("CVresults/",x),sep="\t",header=T,stringsAsFactors=F))
dres = do.call(rbind, lres)


# Overall Results Averaged Across Drugs
#----------------------------------------
lav = lapply(lres, function(x) data.frame(meanRMSD=mean(x$meanRMSD),sdRMSD=sd(x$sdRMSD),
				    meanSpearman=mean(x$meanSpearman),sdmeanSpearman=sd(x$meanSpearman),
				    minSpearman=mean(x$minSpearman),sdminSpearman=sd(x$minSpearman),
				    Kpatient=x$Kpatient[1], Kdrug=x$Kdrug[1]
				    ))
dav = do.call(rbind, lav)

a = dav ; a$Kpatient = factor(a$Kpatient, levels=dav$Kpatient[sort(dav$meanSpearman,index.r=T)$ix])
g1 = ggplot(a,aes(x=Kpatient,y=meanSpearman)) + geom_bar(stat="identity") + 
   geom_errorbar(aes(ymin=meanSpearman-sdmeanSpearman,ymax=meanSpearman+sdmeanSpearman),size=.3,width=.2)+
   coord_flip() + 
   theme0 + ggtitle("Mean Spearman Across Drugs \n of the mean drug CV Spearman")

a = dav ; a$Kpatient = factor(a$Kpatient, levels=dav$Kpatient[sort(dav$minSpearman,index.r=T)$ix])
g2 = ggplot(a,aes(x=Kpatient,y=minSpearman)) + geom_bar(stat="identity") + 
   geom_errorbar(aes(ymin=minSpearman-sdminSpearman,ymax=minSpearman+sdminSpearman),size=.3,width=.2)+
   coord_flip() + 
   theme0 + ggtitle("Mean Spearman Across Drugs \n of the min drug CV Spearman")

a = dav ; a$Kpatient = factor(a$Kpatient, levels=dav$Kpatient[sort(dav$meanRMSD,index.r=T,decreasing=T)$ix])
g3 = ggplot(a,aes(x=Kpatient,y=meanRMSD),) + geom_bar(stat="identity") + 
   geom_errorbar(aes(ymin=meanRMSD-sdRMSD,ymax=meanRMSD+sdRMSD),size=.3,width=.2)+
   coord_flip() + 
   theme0 + ggtitle("Mean RMSD Across Drugs \n of the mean drug CV RMSD")

ga = ggarrange(g1,g2,g3, ncol=3)
system("mkdir -p plots")
ggsave(ga, file="plots/AveragePerf.pdf", width=18, height=6)


# Look at good and bad drugs
#----------------------------------------

# Keep the good kernels:
dgood = dres[as.vector(dres$Kpatient)%in%c("KrnarestrictRbf6","KrnarestrictRbf7","KrnarestrictRbf8","KrnarestrictRbf9","KrnarestrictMean"),]
tmp = dgood[dgood$Kpatient%in%"KrnarestrictRbf9",]

mylevels = tmp$drug[sort(tmp$meanSpearman,index.r=T)$ix]
dgood$drug = factor(dgood$drug, levels=mylevels)

ggd = ggplot(dgood,aes(x=drug,y=meanSpearman)) + geom_bar(stat="identity") +
   geom_errorbar(aes(ymin=meanSpearman-sdSpearman,ymax=meanSpearman+sdSpearman),size=.3,width=.2) +
   coord_flip() + facet_grid(.~Kpatient) +  
   theme0 #+ ggtitle("Mean Spearman Across Drugs \n of the min drug CV Spearman")

ggd2 = ggplot(dgood,aes(x=drug,y=minSpearman)) + geom_bar(stat="identity") +
   geom_errorbar(aes(ymin=minSpearman-sdSpearman,ymax=minSpearman+sdSpearman),size=.3,width=.2) +
   coord_flip() + facet_grid(.~Kpatient) +  
   theme0 #+ ggtitle("Mean Spearman Across Drugs \n of the min drug CV Spearman")

ggsave(ggd, file="plots/DrugPerf.pdf", width=28, height=12)

#ggarrange(ggd,ggd2,ncol=1)

dgood$lowlimit = dgood$meanSpearman - dgood$sdSpearman
dgood$is_good = "decent"
dgood$is_good[dgood$lowlimit < 0.1] = "poor"

ggdg = ggplot(dgood,aes(x=drug,y=meanSpearman)) + geom_bar(stat="identity") +
   geom_errorbar(aes(ymin=meanSpearman-sdSpearman,ymax=meanSpearman+sdSpearman),size=.3,width=.2) +
   coord_flip() + facet_grid(is_good~Kpatient, scales="free_y") +
   theme0 #+ ggtitle("Mean Spearman Across Drugs \n of the min drug CV Spearman")

ggsave(ggdg, file="plots/DrugPerfgood.pdf", width=28, height=12)

rrdg = ggplot(dgood,aes(x=drug,y=meanRMSD)) + geom_bar(stat="identity") +
   geom_errorbar(aes(ymin=meanRMSD-sdRMSD,ymax=meanRMSD+sdRMSD),size=.3,width=.2) +
   coord_flip() + facet_grid(is_good~Kpatient, scales="free_y") +
   theme0 #+ ggtitle("Mean Spearman Across Drugs \n of the min drug CV Spearman")

ggsave(rrdg, file="plots/DrugPerfgoodRMSD.pdf", width=28, height=12)


# Look at good and bad drugs | INTEGRATION??
#----------------------------------------
r = list.files("../../prediction/results_SC1/RNA_cor_enriched/")
r = r[grepl("predictorRF",r)]
rl = lapply(r, function(x) { 
	       y=read.table(paste0("../../prediction/results_SC1/RNA_cor//",x),sep="\t",header=T,stringsAsFactors=F) 
	       y$drug = rownames(y)
	       return(y)
})
drl = do.call(rbind, rl)
rrl = lapply(unique(drl$drug), function(dd) { x=drl[drl$drug==dd,] ; y = data.frame(meanRMSD=mean(x$meanRMSD),sdRMSD=sd(x$sdRMSD),
				    meanSpearman=mean(x$meanSpearman),sdmeanSpearman=sd(x$meanSpearman),drug=dd)
				    } )

# MEAN PERF PER DRUG OF RF ACROSS MANY PARAMETER
dr = do.call(rbind, rrl)
dr$algo = "RF"

# MEAN PERF PER DRUG OF KMR ACROSS DIFFERENT NON-STUPID KERNELS
kkl = lapply(unique(as.vector(dgood$drug)), function(dd) { x=dgood[dgood$drug==dd,] ; y = data.frame(meanRMSD=mean(x$meanRMSD),sdRMSD=sd(x$sdRMSD),
				    meanSpearman=mean(x$meanSpearman),sdmeanSpearman=sd(x$meanSpearman),drug=dd)
				    } )
kr = do.call(rbind, kkl)
kr$algo="KMR"
identical(dr$drug, kr$drug)

dkr = rbind(dr,kr)

dkr$drug = factor(dkr$drug, levels=levels(dgood$drug))
gcomp1 = ggplot(dkr,aes(x=drug,y=meanSpearman,color=algo)) + 
   geom_pointrange(aes(ymin=meanSpearman-sdmeanSpearman,ymax=meanSpearman+sdmeanSpearman),size=.3) +
   coord_flip() + theme0

ggsave(gcomp1, file="plots/compareRain.pdf",height=12)

jbad = which((kr$meanSpearman+kr$sdmeanSpearman) < (dr$meanSpearman-dr$sdmeanSpearman-0.05))
baddrug = as.vector(kr$drug[jbad])

dkr$to_keep = "keep KMR"
dkr$to_keep[dkr$drug%in%baddrug] = "switch RF"

gcomp2 = ggplot(dkr,aes(x=drug,y=meanSpearman,color=algo)) + 
   geom_pointrange(aes(ymin=meanSpearman-sdmeanSpearman,ymax=meanSpearman+sdmeanSpearman),size=.3) +
   coord_flip() + theme0 + facet_grid(.~to_keep)

ggsave(gcomp2, file="plots/compareRain2.pdf",height=12, width=9)

gcomp3 = ggplot(dkr,aes(x=drug,y=meanRMSD,color=algo)) + 
   geom_pointrange(aes(ymin=meanRMSD-sdRMSD,ymax=meanRMSD+sdRMSD),size=.3) +
   coord_flip() + theme0 + facet_grid(.~to_keep)

write.table(data.frame(baddrug) , "drug_switch_kmr_rf.txt" , col.names=F, row.names=F, sep="\t",quote=F)
