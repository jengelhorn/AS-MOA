#!/bin/bash
#run by ./deduplicate_and_translate_hallifted_SNPs.sh genotype inDIR genome_size_file outDIR
EXPECTED_ARGS=4
E_BADARGS=4
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./deduplicate_and_translate_hallifted_SNPs.sh genotype inDIR genome_size_file outDIR "
  exit $E_BADARGS
fi


g=$1
Idir=$2
gen=$3
Odir=$4

export PATH="/netscratch/dep_psl/grp_frommer/Thomas/bin/bedops/bin:$PATH"


 cd ${Idir}

#${g}.B73.hallifted.all.bed


gawk -v OFS='\t' '{print $1$2}' ${g}.B73.hallifted.all.bed | sort | uniq -d > ${g}.B73.hallifted.all.dups_${g}coord.txt

gawk -v OFS='\t' '{print $4}' ${g}.B73.hallifted.all.bed | sort | uniq -d > ${g}.B73.hallifted.all.dups_B73coord.txt



gawk -v OFS='\t' 'NR==FNR{a[$1]=$1; next} !($1$2 in a){print $0}' ${g}.B73.hallifted.all.dups_${g}coord.txt  ${g}.B73.hallifted.all.bed > ${g}.B73.hallifted.all.clean1.bed

gawk -v OFS='\t' 'NR==FNR{a[$1]=$1; next} !($4 in a){print $0}'  ${g}.B73.hallifted.all.dups_B73coord.txt ${g}.B73.hallifted.all.clean1.bed > ${g}.B73.hallifted.all.clean2.bed



mkdir $Odir

echo "renaming chromosomes"

gawk -v OFS='\t' -v ALTn=$g '{print ALTn"-"$1,$2,$3,$4}' ${g}.B73.hallifted.all.clean2.bed | sortBed -g ${gen} > ${Odir}/${g}.B73.ID.hallifted.all.bed

gawk -v OFS='\t'  '{{split($4,var,"."); print "B73-"var[1],var[2]-1,var[2],$1"."$3"."var[4]"."var[3]}}' ${g}.B73.hallifted.all.clean2.bed  | sortBed -g ${gen}  > ${Odir}/B73.${g}.ID.hallifted.all.bed
