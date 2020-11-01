##--------------------------------------
## Kernel Ridge Regression
## Try all the kernels combinaisons
##--------------------------------------
args <- commandArgs(TRUE)
patientKernel <- args[1]
drugKernel <- args[2]
mycores <- args[3]

source("../../src/Kpredictors_SC1.R")
source("../../src/Kevaluate_SC1.R")

#patientKernel="KrnarestrictLinear.txt"
#drugKernel="Kempirical.txt"

path="../../data/kernels/"

Kpatient = read.table(paste0(path,patientKernel), sep="\t", header=T, stringsAsFactor=F)
dim(Kpatient)

Kdrug = read.table(paste0(path,drugKernel), sep="\t", header=T, stringsAsFactor=F)
dim(Kdrug)

response = read.table("../../data/responses/prepared_data_auc.csv", sep="\t", header=T, stringsAsFactors=F)

print("Start CV")
myname = paste(gsub(".txt","",patientKernel), gsub(".txt","",drugKernel), sep="_")
print(myname)


# ! ...NOTE THAT THIS IS 5 folds CV with only 2 repreats... !  TODO 10
# ! ...NOTE THAT THIS IS INTERNALE 4 folds CV with only 21 repreats... !  TODO 2


rescv = KevaluateCV(myKpredictor = KpredictorKMR ,
		    celllines = Kpatient ,
		    response = response ,
		    KernelTask = as.matrix(Kdrug) ,
		    mc.cores=mycores ,
		    nfolds=5 ,
		    nrepeats=2 ,
		    KpatientGiven=TRUE ,
		    KpatientType="precomputed" ,
		    nfolds.intern = 4 ,
		    nrepeats.intern = 1 ,
		    lambdas = 10^(-10:10)
)

if (1==0) {
rescv = KevaluateCV(myKpredictor = KpredictorKMR ,
		    celllines = Kpatient ,
		    response = response[,1:50] ,
		    KernelTask = as.matrix(Kdrug[1:50,1:50]) ,
		    mc.cores=1 ,
		    nfolds=5 ,
		    nrepeats=1 ,
		    KpatientGiven=TRUE ,
		    KpatientType="precomputed" ,
		    nfolds.intern = 4 ,
		    nrepeats.intern = 1 ,
		    lambdas = 10^(-10:10)
)
}


resdf = as.data.frame(do.call(rbind,rescv))
resdf$drug = rownames(resdf)
resdf$Kpatient = gsub(".txt","",patientKernel)
resdf$Kdrug = gsub(".txt","",drugKernel)

system("mkdir -p CVresults")
write.table(resdf, file=paste0(myname,".csv"), quote=F, sep="\t")
write.table(resdf, file=paste0("CVresults/",myname,".txt"), quote=F, sep="\t")
