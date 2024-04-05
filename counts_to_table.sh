#!/bin/bash
#run by ./counts_to_table.sh genotype dir cond1 cond2 Peaks_cond1 Peaks_cond2 bedgraph_condition1 bedgraph_condition2 SNPs.tsv
EXPECTED_ARGS=9
E_BADARGS=9
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./counts_to_table.sh genotype dir cond1 cond2 Peaks_cond1 Peaks_cond2 bedgraph_condition1 bedgraph_condition2 SNPs.tsv "
  exit $E_BADARGS
fi


g=$1
dir=$2
C1=$3
C2=$4
PeakC1=$5
PeakC2=$6
bed1=$7
bed2=$8
tsv=$9

cd ${dir}

#Start with : Chr start stop ID Value

#To add peak information

echo "sort and name peaks"

#Peaks are taken from Narrowpeak files of MACS3

for tr in ${C1};  do sort -k1,1 -k2,2n ${PeakC1} | gawk -v OFS='\t' '{if($2<$3 && $1!~"scaf-alt"){print $1,$2,$3,$4"_"$5}}' > ${g}.${tr}.peaks.names.sorted.bed ; done

for tr in ${C2};  do sort -k1,1 -k2,2n ${PeakC2} | gawk -v OFS='\t' '{if($2<$3 && $1!~"scaf-alt"){print $1,$2,$3,$4"_"$5}}' > ${g}.${tr}.peaks.names.sorted.bed ; done

echo "sort count files"

for tr in ${C1} ${C2}; do gawk -v OFS='\t' -v g=$g '{if($4!~"scaf-alt"){print $0} else{gsub("scaf-alt","scaf",$4);print $1,$2,$3,$4,"het"}}' ${dir}/B73.${g}.${tr}.counts.bed | sort -k1,1 -k2,2n  > B73.${g}.${tr}.counts.sorted.bed; done

for tr in ${C1} ${C2};  do gawk -v OFS='\t' -v g=$g '{if($1!~"scaf-alt"){print $0} else{gsub("scaf-alt","scaf",$1);print $1,$2,$3,$4,"het"}}' ${dir}/${g}.B73.${tr}.counts.bed | sort -k1,1 -k2,2n  > ${g}.B73.${tr}.counts.sorted.bed; done



#What we want: Chr(B73)  STop=POS(B73)   ID(NAMchr_NAMPOS)   REFal(B73)   ALTal(NAM)   Count(Either B73 or NAM)    and we want this B73 coordinate sorted
echo "intersect peaks and sort"

for tr in ${C1} ${C2}; do intersectBed -a B73.${g}.${tr}.counts.sorted.bed -b ${g}.${tr}.peaks.names.sorted.bed -wa -wb -sorted -loj| gawk -v OFS='\t' '{split($4,var,"\."); print $1,$3,var[1]":"var[2],var[4], var[3], $5, $9}' >  B73.${g}.${tr}.ID.count.sort.csv; done

for tr in ${C1} ${C2}; do intersectBed -a ${g}.B73.${tr}.counts.sorted.bed -b ${g}.${tr}.peaks.names.sorted.bed -wa -wb -sorted -loj|gawk -v OFS='\t' '{split($4,var,"\."); print var[1],var[2],$1":"$3,var[3], var[4], $5, $9}' | sort -k1,1 -k2,2n  > ${g}.B73.${tr}.ID.count.sort.csv; done



#Now make new file with both ref and alt, one WW, one DS:

# Chr(B73)  STop=POS(B73)   ID(NAMchr_NAMPOS)   REFal(B73)   ALTal(NAM)   Count(B73)   Count(NAM) Peak(B73) Peak(NAM)  PF(PF, if both 0 print n.r.)

echo "preparing postfrequency files"

for tr in ${C1} ${C2}; do paste B73.${g}.${tr}.ID.count.sort.csv ${g}.B73.${tr}.ID.count.sort.csv | gawk -v OFS='\t' '{if($2==$9){if($6+$13==0){print "B73-"substr($1,5,length($1)-4),$2,$10,$4,$5,$6,$13,$7,$14,"n.r."} else {print "B73-"substr($1,5,length($1)-4),$2,$10,$4,$5,$6,$13,$7,$14,$6/($6+$13)}}}' > B73.${g}.${tr}.PF.csv; done

echo "adding genotypes and depth "

#Add column with genotype, i.e. if that SNP is present in the List

#Chr(B73)  STop=POS(B73) REFal(B73)   ALTal(NAM) ID(NAMchr_NAMPOS) Genotype (0/0 if B73, 1/1 if NAM) Count(B73)   Count(NAM) PF (n.r. if both 0)

#Add depth information, i.e. calculate number of reads B73 and NAM and change PF to n.p. post frequency if none of the alleles has a peak at this position

#Chr(B73)  STop=POS(B73) REFal(B73)  ALTal(NAM) ID(NAMchr_NAMPOS) Genotype (0/0 if B73, 1/1 if NAM) norm_Count(B73)   norm_Count(NAM) PEAK(B73) Peak(NAM)  PF (n.p. if no peak in either B73 or NAM) read_Count(B73) read_Count(NAM)

#Determine value per read first

CW=`zcat ${bed1} | awk -v OFS="\t" 'BEGIN{n=10000} {if($4>0 && $4<n){n=$4}} END { print n }' `;
CD=`zcat ${bed2} | awk -v OFS="\t" 'BEGIN{n=10000} {if($4>0 && $4<n){n=$4}} END { print n }' `;

for tr in ${C1} ; do gawk -v OFS='\t' -v g=$g -v CW=${CW} '{if(NR==FNR) {a[$1"ü"$2]=$1"ü"$2; text[$1"ü"$2]=$0; next} if(substr($1,5,length$1-4)"ü"$2-1 in a){if($8~"peak" || $9~"peak"){print $1,$2,$4,$5,$3,"1/1",$6,$7,$8,$9,$10,int(($6/CW)+0.5),int(($7/CW)+0.5)} else {print $1,$2,$4,$5,$3,"1/1",$6,$7,$8,$9,"n.p."$10,int(($6/CW)+0.5),int(($7/CW)+0.5)}} else {if($8~"peak" || $9~"peak"){print $1,$2,$4,$5,$3,"0/0",$6,$7,$8,$9,$10,int(($6/CW)+0.5),int(($7/CW)+0.5)} else{print $1,$2,$4,$5,$3,"0/0",$6,$7,$8,$9,"n.p."$10,int(($6/CW)+0.5),int(($7/CW)+0.5)}}}' ${tsv} B73.${g}.${tr}.PF.csv  > B73.${g}.${tr}.PF.GT.RN.csv; done

for tr in ${C2} ; do gawk -v OFS='\t' -v g=$g -v CD=${CD} '{if(NR==FNR) {a[$1"ü"$2]=$1"ü"$2; text[$1"ü"$2]=$0; next} if(substr($1,5,length$1-4)"ü"$2-1 in a){if($8~"peak" || $9~"peak"){print $1,$2,$4,$5,$3,"1/1",$6,$7,$8,$9,$10,int(($6/CD)+0.5),int(($7/CD)+0.5)} else {print $1,$2,$4,$5,$3,"1/1",$6,$7,$8,$9,"n.p."$10,int(($6/CD)+0.5),int(($7/CD)+0.5)}} else {if($8~"peak" || $9~"peak"){print $1,$2,$4,$5,$3,"0/0",$6,$7,$8,$9,$10,int(($6/CD)+0.5),int(($7/CD)+0.5)} else{print $1,$2,$4,$5,$3,"0/0",$6,$7,$8,$9,"n.d."$10,int(($6/CD)+0.5),int(($7/CD)+0.5)}}}' ${tsv} B73.${g}.${tr}.PF.csv  > B73.${g}.${tr}.PF.GT.RN.csv; done

