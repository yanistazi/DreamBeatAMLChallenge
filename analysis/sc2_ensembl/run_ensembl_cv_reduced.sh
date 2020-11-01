#!/bin/sh

echo "GLM-GLM"
bsub -J glmglm -n 10 -oo oo_glmglm -eo eo_glmglm /work/isabl/bin/Rscript_3.6.1 ensembl_GLM_GLM_reduced.R


echo "GLM-RF"
bsub -J glmrf -n 10 -oo oo_glmrf -eo eo_glmrf /work/isabl/bin/Rscript_3.6.1 ensembl_GLM_RF_reduced.R


echo "RF-GLM"
bsub -J rfglm -n 10 -oo oo_rfglm -eo eo_rfglm /work/isabl/bin/Rscript_3.6.1 ensembl_RF_GLM_reduced.R


echo "RF-RF"
bsub -J rfrf -n 10 -oo oo_rfrf -eo eo_rfrf /work/isabl/bin/Rscript_3.6.1 ensembl_RF_RF_reduced.R
