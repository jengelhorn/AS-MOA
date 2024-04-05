# AS-MOA
Scripts for allele-specific analysis of transcription factor (TF) binding
## Description 
These scripts were developed to count TF binding data (e.g. generated with MOA-seq) on corresponding positions in F1 hybrids. They assume omni hybrids with B73 as a common mother, so for the moment B73 as one allele is hard coded but this will be flexibilised in the future. We further expect TF binding to be analysed in two conditions, e.g. well-watered and drought. If only one condition was analysed, all parts with cond2 can be commented out and cond2 in the commands can be specified as 0.


-------------------
## halliftover_of_SNPs.sh 

Example commant: 
```bash
./halliftover_of_SNPs.sh genotype output-directory hal-file SNP-file
```


This program generates corresponding coordinates of B73 Positions to be analysed in the genome of the other parent, called "genotype"

- expects a hal file where one genome is called B73 and the other has the same name as specified in "genotype", e.g. Mo17
- expects an SNP file in the format of XXX
- generates an output [genotype].B73.hallifted.all.bed in bed formate XXX

To do: Verify from server, check if PATH is needed, improve annotation

---------------------
## deduplicate_and_translate_hallifted_SNPs.sh

Usage: ./deduplicate_and_translate_hallifted_SNPs.sh genotype inDIR genome_size_file outDIR 

This program removes duplicated positions from the [genotype].B73.hallifted.all.bed file, which indicate a position in one genome is ambiguous in the other. It also adds identifiers to the chromosome names of the parental genomes and generates IDs of the paternal and maternal positions for translation.

- expects [genotype].B73.hallifted.all.bed file in inDIR
- expects a genome file size for a diploid genome in the form CHROM LENGTH for each Chromosome, where the two parental chromosomes are named B73-[chromosomename] and [genotype]-[chromosomename], chromosomename being the original names used in the hal and bed files in the previous steps
- generates an output B73.[genotype].ID.hallifted.all.bed in bed formate XXX

To do: Verify from server, check if PATH is needed, improve annotation

----------------------
The output of deduplicate_and_translate_hallifted_SNPs.sh can be used to exclude positions that are not 1:1 mappable in at least two lines from the further analysis. For this, a list of all 1:1 mappable in two lines is generated:

for g in [genotype1] [genotype2] ...; do cat B73.${g}.ID.hallifted.all.bed >> all.ID.hallifted.all.bed; done

gawk -v OFS='\t' '{print $1,$2,$3}' all.ID.hallifted.all.bed |sort| uniq -c > all.ID.hallifted.all.count.csv

gawk -v OFS='\t' '{if($1>1){print $2,$3,$4}}'  all.ID.hallifted.all.count.csv > All_SNPs_1_to_1_mapp_min_2_lines.bed

This last file will be employed in the next program

-------------------------
## hallifted_bed_to_count.sh 

Usage: ./hallifted_bed_to_count.sh genotype dir 1:1map_SNPs genome_size_file bedgraph_condition1 bedgraph_condition2 cond1 cond2

First, this script removes all SNPs not found in at least two lines from the analysis. If this is not desired, the original SNP file can be specified as 1:1map_SNPs. The script then counts the signal on each specified position from the SNP file in both parental genomes and pastes the information together.

- expects a working directory containing B73.[genotype].ID.hallifted.all.bed files from deduplicate_and_translate_hallifted_SNPs.sh
- expects a SNP file in the format XXX 
- expects a genome-size file in the same format as mentioned above
- expects gzipped bedgraph files containing the counts of an experiment mapped to the diploid genome
- expects names for two conditions, e.g. WW and DS which should stay the same for all following scripts
- generates four files: B73.[genotype].[cond1].counts.bed  [genotype].B73.[cond1].counts.bed B73.[genotype].[cond2].counts.bed  [genotype].B73.[cond2].counts.bed in the formate XXX

To do: improve annotation, maybe accommodate non-zipped files

--------------------------
## counts_to_table.sh

Usage: ./counts_to_table.sh genotype dir cond1 cond2 Peaks_cond1 Peaks_cond2 bedgraph_condition1 bedgraph_condition2 SNPs.tsv

This script compiles counting information from the two genomes obtained by hallifted_bed_to_count.sh, calculates binding frequencies (BF: B73 counts/ (B73 counts + genotype counts)), read counts per allele, adds information about the position being located in a peak and the presence of an SNP at each position in the respective hybrid (note that values are counted for all SNPs/Pos that are fed in in the beginning (e.g. all biallecic, 1:1 mappable SNPs in the population analysed, not only SNPs in the specific hybrid).  

- expects peak files of the formate XXX
- expects to find the four output files of hallifted_bed_to_count.sh in the working directory [dir]
- expects gzipped bedgraph files containing the counts of an experiment mapped to the diploid genome (they are used to convert the normalised values back into read count values)
- expects a list of positions that carry an SNP between the two alleles of the specific hybrid in tsv formate from cactus (0-based)
- generates B73.[genotype].[condition].PF.GT.RN.csv files in the formate Chr(B73)  STop=POS(B73) REFal(B73)  ALTal ID([genotype]chr_[genotype]POS) Genotype (0/0 if [genotype] allele is B73, 1/1 if [genotype] allele is ALTal) normalised_Count(B73)   normalised_Count([genotype]) PEAK(B73) Peak([genotype])  BF (n.p. if no peak in either B73 or [genotype]) read_Count(B73) read_Count([genotype])

To do: check annotation, improve naming of peaks so other formats can be included, correct n.d. to n.p., include next step here with adding a read/normalised count cutoff and only get 1/1 into the command

--------------------------

for g in Ki3	Ki11 A619	B97   CML277	CML322	CML333	CML69	HP301	IL14H Ky21	M162W	Mo17	Mo18W	NC358	Oh43	;do for tr in WW DS; do gawk -v OFS='\t' -v g=$g 'BEGIN{print "NR","Chr","Pos","REFal","ALTal","PosNAM","Genotype","EGcount_B73","EG_counts_NAM","Peak_B73","Peak_NAM","PF","Counts_B73","Counts_NAM"}{if($6=="1/1" && $0~"peak"){if($11>0 && $11<1){if($7>7 || $8>7){print NR,$0}}}}' /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/custom/counts_STAR_80_EG/${g}/new_25_lines_22_12/B73.${g}.${tr}.q255.PF.GT.RN.csv > ${g}.${tr}.PF.q255.GT.SNPs.CPM7.txt; done; done

----------------------------------------------
To perform binomial testing of allele-specific binding:

Rscript --vanilla Binomial_fdr_SNPs_1222.R ${g}.${tr}.PF.q255.GT.SNPs.CPM7.txt ${g}.${tr}.q255.bino.fdr.CPM7.txt; done; done

Adjust names accordingly

----------------------------------------------
To account for mapping bias, a control can be used to eliminate those significant sites that are also biased in control data, e.g. short-read sequencing genomic data. Reads of this data should be shortened to match the read length of the original experiment.

Counting of control reads is performed for the same positions as the experimental data, thus halliftover and translation of SNPs does not have to be repeated. The working directory should be the same as for the original analysis because the outputs of the previous steps are required.

--------------------------------
No script yet for the counting, only commands:

These commands take the list of eligible sites for allele-specific binding (all that went into the binomial test, in our case all inside a peak, with a normalised read count larger than 7 and at least one read in each allele), creates a window in which the control data should be checked (e.g. 65 bp for MOA) and counts the mean normalised read count in this window for each of the two alleles. This is done for one condition only in our case (WW) to have the same threshold for both conditions. For cond2 (DS) we only count the control data on the allele-specific sites to use this data later to exclude sites with above threshold deviations from 1:1 in the control data.


for g in  [genotype]	; do gawk -v OFS='\t' -v g=$g '{if(NR>1){print substr($3,2,length($3)-2),$4-33,$4+32,substr($7,2,length($7)-2)}}' ../${g}.WW.q255.bino.fdr.CPM7.txt |  gawk -v OFS='\t' '{if($2<0){print $1,0,$3,$4}else{print $0}}' | sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt > B73.${g}.65bp.bed; done


for g in [genotype] ; do gawk -v OFS='\t' -v g=$g '{if(NR>1){split($7,pos,":"); print substr(pos[1],2,length(pos[1])-1),substr(pos[2],1,length(pos[1])-1)-33,substr(pos[2],1,length(pos[2])-1)+32,substr($7,2,length($7)-2)}}' ../${g}.WW.q255.bino.fdr.CPM7.txt  |  gawk -v OFS='\t' '{if($2<0){print $1,0,$3,$4}else{print $0}}'| sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt > ${g}.65bp.bed; done


for g in  [genotype]; do gawk -v OFS='\t' -v g=$g '{if(NR>1){print substr($3,2,length($3)-2),$4-33,$4+32,substr($7,2,length($7)-2)}}' ../${g}.DS.q255.bino.fdr.01.CPM7.txt |  gawk -v OFS='\t' '{if($2<0){print $1,0,$3,$4}else{print $0}}' | sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt > B73.${g}.DS.01.65bp.bed; done

for g in [genotype]	; do gawk -v OFS='\t' -v g=$g '{if(NR>1){split($7,pos,":"); print substr(pos[1],2,length(pos[1])-1),substr(pos[2],1,length(pos[1])-1)-33,substr(pos[2],1,length(pos[2])-1)+32,substr($7,2,length($7)-2)}}' ../${g}.DS.q255.bino.fdr.01.CPM7.txt  |  gawk -v OFS='\t' '{if($2<0){print $1,0,$3,$4}else{print $0}}'| sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt > ${g}.DS.01.65bp.bed; done


for g in IL14H; do time gawk -v OFS='\t' '{if($1!~"B73"){a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3-i),($3-i+1), $4}} else{ print $0 }}}' /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bamCoverage/CTRL/EG_norm/B73.${g}.CTRL.RPGCq3.b1.bedgraph| bedtools map -a ./new/${g}.DS.01.65bp.bed -b - -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt -c 4 -o mean > ./new/${g}.DS.01.CTRL.65bp.bed; done


for g in 		IL14H ; do time gawk -v OFS='\t' '{if($1!~"B73"){a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3-i),($3-i+1), $4}} else{ print $0 }}}' /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bamCoverage/CTRL/EG_norm/B73.${g}.CTRL.RPGCq3.b1.bedgraph| bedtools map -a ${g}.65bp.bed -b - -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt -c 4 -o mean > ./new/${g}.CTRL.65bp.bed; done



for g in	IL14H	; do time gawk -v OFS='\t' '{if($1~"B73"){a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3-i),($3-i+1), $4}} else{ print $0 }}}' /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bamCoverage/CTRL/EG_norm/B73.${g}.CTRL.RPGCq3.b1.bedgraph | bedtools map -a B73.${g}.DS.01.65bp.bed -b - -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt -c 4 -o mean > B73.${g}.DS.01.CTRL.65bp.bed; done


for g in	IL14H; do time gawk -v OFS='\t' '{if($1~"B73"){a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3-i),($3-i+1), $4}} else{ print $0 }}}' /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bamCoverage/CTRL/EG_norm/B73.${g}.CTRL.RPGCq3.b1.bedgraph | bedtools map -a B73.${g}.65bp.bed -b - -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt -c 4 -o mean > B73.${g}.CTRL.65bp.bed; done


To do: write script for this part 

---------------------------------
## clean_AFPs_with_CNTRL.sh

Usage: ./clean_AFPs_with_CNTRL.sh genotype dir Control_counts_B73_cond1 Control_counts_[genotype]_cond1 Binom_test_output_cond1 Cond1 Cond2 Control_counts_B73_cond2 Control_counts_[genotype]_cond2 Binom_test_output_cond2 outdir 

This script takes the control values from both alleles and calculates a binding frequency for each window (around an MP position for WW and AMP for DS). It then filters for non-AMP MPs and calculates the BF cut-off to exclude the upper and lower 5 percentile of biased binding. Thus, an upper and lower control BF cut-off is calculated. All AMPs with control BF above the upper and below the lower threshold are then filtered out.

- expects four bed files containing control counts for MP windows for cond1 for B73 and [genotype] and AMP windows for cond2.
- expects the binomial test output for both conditions
- expects names of the conditions (e.g. WW and DS)

To do: improve comments, add cond1 to usage statement
