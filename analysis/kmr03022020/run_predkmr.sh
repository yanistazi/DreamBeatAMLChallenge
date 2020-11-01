#!/bin/sh

#allKpatient=$(ls ../../data/kernels/ | grep KrnarestrictRb)
#allKpatient="KrnarestrictLinear.txt Kgenetic.txt Kdirac.txt Kuniform.txt"

allKrbf=$(ls ../../data/kernels/ | grep KrnarestrictRb)

allKpatient="$allKrbf KrnarestrictMean.txt KrnarestrictMean510.txt"

allKdrug="Kempirical.txt"

mcores=10

for Kpatient in $allKpatient
do
  for Kdrug in $allKdrug
  do
    name=${Kpatient}_${Kdrug}
    echo ${name}
    bsub -J ${name} -n ${mcores} -eo eo_${name} -oo oo_${name} "Rscript PredKMR_rscript.R ${Kpatient} ${Kdrug} ${mcores}"
  done
done
