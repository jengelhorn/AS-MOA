#!/bin/bash
#run by ./hallifted_bed_to_count.sh genotype dir 1:1map_SNPs genome_size_file bedgraph_condition1 bedgraph_condition2 cond1 cond2
EXPECTED_ARGS=8
E_BADARGS=8
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./hallifted_bed_to_count.sh genotype dir 1:1map_SNPs genome_size_file bedgraph_condition1 bedgraph_condition2 cond1 cond2 "
  exit $E_BADARGS
fi


g=$1
dir=$2
SNPs=$3
gen=$4
bed1=$5
bed2=$6
C1=$7
C2=$8

 cd ${dir}
 
 
echo "retain only deduplicated, biallelic SNPs that occur in two lines from hallifted file"

gawk -v OFS='\t' 'NR==FNR{a[$1$2]=$1$2; next} ($1$2 in a){print $0}' $SNPs B73.${g}.ID.hallifted.all.bed > B73.${g}.ID.hallifted.all.final.bed 
 
gawk -v OFS='\t' 'NR==FNR{a[$1$3]=$1$3; next} {split($4,pos,".")};("B73-"pos[1]pos[2] in a){print $0}' $SNPs ${g}.B73.ID.hallifted.all.bed > ${g}.B73.ID.hallifted.all.final.bed

echo "bedtool mapping, counting Signal at each SNP position in both alleles"

zcat ${bed1} | bedtools map -a B73.${g}.ID.hallifted.all.final.bed -b - -g ${gen} -c 4 > B73.${g}.${C1}.counts.bed

zcat ${bed1} | bedtools map -a ${g}.B73.ID.hallifted.all.final.bed -b - -g ${gen} -c 4 > ${g}.B73.${C1}.counts.bed

zcat ${bed2} | bedtools map -a B73.${g}.ID.hallifted.all.final.bed -b - -g ${gen} -c 4 > B73.${g}.${C2}.counts.bed

zcat ${bed2} | bedtools map -a ${g}.B73.ID.hallifted.all.final.bed -b - -g ${gen} -c 4 > ${g}.B73.${C2}.counts.bed

