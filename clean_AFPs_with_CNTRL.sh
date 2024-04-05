#!/bin/bash
#run by ./clean_AFPs_with_CNTRL.sh genotype dir Control_counts_B73 Control_counts_NAM Binom_test_output Cond1 Cond2 Control_counts_B73_cond2 Control_counts_NAM_cond2 Binom_test_output_cond2 outdir

EXPECTED_ARGS=11
E_BADARGS=11
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./clean_AFPs_with_CNTRL.sh genotype dir Control_counts_B73 Control_counts_NAM Binom_test_output Cond1 Cond2 Control_counts_B73_cond2 Control_counts_NAM_cond2 Binom_test_output_cond2 outdir "
  exit $E_BADARGS
fi


g=$1
dir=$2
C1=$3
#bedfile containing window which was counted in, e.g. 65 bp and NAM coordinate plus mean count in control data for B73
C2=$4
#bedfile containing window which was counted in, e.g. 65 bp and NAM coordinate plus mean count in control data for NAM
Bino=$5
#Output file of Binomial_fdr_SNPs_1222.R for cond1
cond1=$6
cond2=$7
C3=$8
#bedfile containing window which was counted in, e.g. 65 bp and NAM coordinate plus mean count in control data for B73, in the case here the condition 2 non-significant SNPs are not used to determine the threshold so it is sufficient to just give AFPs with counts here
C4=$9
#same as C3 but for NAM coordinate
Bino2=${10}
#Output file of Binomial_fdr_SNPs_1222.R for cond2
out=${11}


cd ${dir}

echo "Paste CNTRL files for "$g

gawk -v OFS='\t' 'NR==FNR {a[$4]=$4; text[$4]=$0; next} $4 in a {print text[$4],$0}' $C1 $C2 |  gawk -v OFS='\t' '{if($5+$10>0) print $0, ($5/($5+$10))}' > ${g}.CTRL_both.bed

echo "Add AFP info back"

gawk -v OFS='\t' 'NR==FNR {a[$4]=$4; text[$4]=$11; next} substr($7,2,length($7)-2) in a {print $0, text[substr($7,2,length($7)-2)]}' ${g}.CTRL_both.bed $Bino | sort -k18,18n > ${g}.CTRL_both.AFP_info.csv


echo "Count nonAFP FPs and find upper and lower percentile and cut-off values"


FP=`gawk -v OFS='\t' -v g=$g 'BEGIN{n=0}{if($17>=0.01){n++}}END{print n}' ${g}.CTRL_both.AFP_info.csv `;
UP=`gawk -v OFS='\t' -v FP=$FP 'BEGIN{print int(FP*0.05+0.5)}' `;
DOWN=`gawk -v OFS='\t' -v FP=$FP 'BEGIN{print int(FP*0.95+0.5)}' `;

echo "5 percent top:" $UP
echo "5 percent bottom:" $DOWN



UPT=`gawk -v OFS='\t' -v g=$g -v UP=$UP 'BEGIN{n=0}{if($17>=0.01){n++;if(n==UP){print $18}}}' ${g}.CTRL_both.AFP_info.csv `;
echo "Upper 5 percentile limit:" $UPT
DOWNT=`gawk -v OFS='\t' -v g=$g -v DOWN=$DOWN 'BEGIN{n=0}{if($17>=0.01){n++;if(n==DOWN){print $18}}}' ${g}.CTRL_both.AFP_info.csv `;

echo "Lower 5 percentile limit:" $DOWNT

echo "Paste CNTRL files cond2"

gawk -v OFS='\t' 'NR==FNR {a[$4]=$4; text[$4]=$0; next} $4 in a {print text[$4],$0}' $C3 $C4 |  gawk -v OFS='\t' '{if($5+$10>0) print $0, ($5/($5+$10))}' > ${g}.CTRL_both.${cond2}.bed

echo "Add AFP info back cond2"

gawk -v OFS='\t' 'NR==FNR {a[$4]=$4; text[$4]=$11; next} substr($7,2,length($7)-2) in a {print $0, text[substr($7,2,length($7)-2)]}' ${g}.CTRL_both.${cond2}.bed $Bino2 | sort -k18,18n > ${g}.CTRL_both.AFP_info.${cond2}.csv


echo "Preparing cleaned files"

 gawk -v OFS='\t' -v DOWNT=$DOWNT -v UPT=$UPT '{if($17<0.01){if($18<DOWNT && $18>UPT)print $0}}' ${g}.CTRL_both.AFP_info.csv > ${out}/${g}.AFPs.cleaned.${cond1}.csv
 
  gawk -v OFS='\t' -v DOWNT=$DOWNT -v UPT=$UPT '{if($17<0.01){if($18<DOWNT && $18>UPT)print $0}}' ${g}.CTRL_both.AFP_info.${cond2}.csv > ${out}/${g}.AFPs.cleaned.${cond2}.csv
