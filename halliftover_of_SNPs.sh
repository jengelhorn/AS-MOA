#!/bin/bash
#run by ./halliftover_of_SNPs.sh genotype ouput-directory hal-file SNP-file
EXPECTED_ARGS=4
E_BADARGS=4
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./halliftover_of_SNPs.sh genotype ouput-directory hal-file SNP-file "
  exit $E_BADARGS
fi


g=$1
dir=$2
hal=$3
SNPs=$4

export PATH="/netscratch/dep_psl/grp_frommer/Thomas/bin/bedops/bin:$PATH"

 cd ${dir}

##tsv file contains 0-based positions

echo "halLiftover"
 
halLiftover  --hdf5InMemory ${hal} B73 ${SNPs} ${g} ${g}.B73.hallifted.all.bed






