# AS-MOA
Scripts for allele-specific analysis of transcription factor (TF) binding
## Description 
These scripts were developed to count TF binding data (e.g. generated with MOA-seq) on corresponding positions in F1 hybrids. They were developed in Maize assuming omni hybrids with B73 as a common mother, so for the moment B73 as one allele is hard coded but this will be flexibilised in the future. The variable parent is specified as "genotype" in the following commands. We further expect TF binding to be analysed in two conditions, e.g. well-watered and drought. If only one condition was analysed, all parts with cond2 can be commented out and cond2 in the commands can be specified as 0.
## Install
Most functions are standard on a Linux system. The scripts are in bash and use mainly basic functions. External functions called include 
- halliftover (https://github.com/ComparativeGenomicsToolkit/hal)
- bedtools (https://bedtools.readthedocs.io/en/latest/index.html)
- R (libraries dplyr and broom) 

-------------------
## halliftover_of_SNPs.sh 

Example command: 
```bash
./halliftover_of_SNPs.sh genotype output-directory hal-file SNP-file
```


This program generates corresponding coordinates of B73 Positions to be analysed in the genome of the other parent, called "genotype"

- expects a hal file where one genome is called B73 and the other has the same name as specified in "genotype", e.g. Mo17
- expects an SNP file in the format of a bed file with B73 (reference) coordinates with Chr.Position.RefrenceAllele.AlternativeAllele as a 4th column:
 ``` 
chr1    15      16      chr1.16.T.A
chr1    16      17      chr1.17.A.C
chr1    41      42      chr1.42.T.G
chr1    51      52      chr1.52.T.C
```
- generates an output [genotype].B73.hallifted.all.bed in bed formate with the positions of SNPs in the paternal genome with ChrB73.PositionB73.RefrenceAllele.AlternativeAllele as 4th column:
 ```
chr9    68048   68049   chr1.4721.C.T
chr9    68033   68034   chr1.4736.A.G
chr9    68031   68032   chr1.4760.A.C
chr9    68006   68007   chr1.4785.T.G
 ``` 

---------------------
## deduplicate_and_translate_hallifted_SNPs.sh

Example command: 
```bash
./deduplicate_and_translate_hallifted_SNPs.sh genotype inDIR genome_size_file outDIR 
```
This program removes duplicated positions from the [genotype].B73.hallifted.all.bed file, which indicate a position in one genome is ambiguous in the other. It also adds identifiers to the chromosome names of the parental genomes and generates IDs of the paternal and maternal positions for translation.

- expects [genotype].B73.hallifted.all.bed file in inDIR
- expects a genome file size for a diploid genome in the form CHROM LENGTH for each Chromosome, where the two parental chromosomes are named B73-[chromosomename] and [genotype]-[chromosomename], chromosomename being the original names used in the hal and bed files in the previous steps
- outDIR will be created
- generates an output B73.[genotype].ID.hallifted.all.bed and [genotype].B73.ID.hallifted.all.bed in bed formate one with B73 coordinates and paternal position as ChrPat.PositionPat.RefrenceAllele.AlternativeAllele and one with the paternal position and ChrB73.PositionB73.RefrenceAllele.AlternativeAllele

----------------------
## Generate a file with 1:1 mappable SNPs
This step is only necessary when working with several hybrids to generate a pan-cistrome. When working with one hybrid, the SNP file 

The output of deduplicate_and_translate_hallifted_SNPs.sh can be used to exclude positions that are not 1:1 mappable in at least two lines from the further analysis. For this, a list of all 1:1 mappable positions in two lines is generated:

for g in [genotype1] [genotype2] ...; do cat B73.${g}.ID.hallifted.all.bed >> all.ID.hallifted.all.bed; done

gawk -v OFS='\t' '{print $1,$2,$3}' all.ID.hallifted.all.bed |sort| uniq -c > all.ID.hallifted.all.count.csv

gawk -v OFS='\t' '{if($1>1){print $2,$3,$4}}'  all.ID.hallifted.all.count.csv > All_SNPs_1_to_1_mapp_min_2_lines.bed

This last file will be employed in the next program

-------------------------
## hallifted_bed_to_count.sh 

Example command: 

```bash
./hallifted_bed_to_count.sh genotype dir 1:1map_SNPs genome_size_file bedgraph_condition1 bedgraph_condition2 cond1 cond2
```

First, this script removes all SNPs not found in at least two lines from the analysis. If this is not desired, the original SNP file can be specified as 1:1map_SNPs. The script then counts the signal on each specified position from the SNP file in both parental genomes and pastes the information together.

- expects a working directory containing B73.[genotype].ID.hallifted.all.bed files from deduplicate_and_translate_hallifted_SNPs.sh
- expects a SNP file in bed file format with the reference Prefix in the chromosome, e.g. B73-chr10       10000000        10000001
- expects a genome-size file in the same format as mentioned above
- expects gzipped bedgraph files containing the counts of an experiment mapped to the diploid genome, we assume normalised data here mapped to a concatenated genome. For determination of allele-specific binding sites we recommend using only reads that map exactly once to this genome. We recommend normalising as reads per genome coverage e.g. with bamcoverage (https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)
- expects names for two conditions, e.g. WW and DS which should stay the same for all following scripts
- generates four files: B73.[genotype].[cond1].counts.bed  [genotype].B73.[cond1].counts.bed B73.[genotype].[cond2].counts.bed  [genotype].B73.[cond2].counts.bed in bed formate with the count as a 4th column


--------------------------
## counts_to_table.sh

Example command: 



```bash
./counts_to_table.sh genotype dir cond1 cond2 Peaks_cond1 Peaks_cond2 bedgraph_condition1 bedgraph_condition2 SNPs.tsv
```

This script compiles counting information from the two genomes obtained by hallifted_bed_to_count.sh, calculates binding frequencies (BF: B73 counts/ (B73 counts + genotype counts)), read counts per allele, adds information about the position being located in a peak and the presence of an SNP at each position in the respective hybrid (note that values are counted for all SNPs/Pos that are fed in in the beginning (e.g. all biallecic, 1:1 mappable SNPs in the population analysed, not only SNPs in the specific hybrid).  

- expects peak files in bed formate
- expects to find the four output files of hallifted_bed_to_count.sh in the working directory [dir]
- expects gzipped bedgraph files containing the counts of an experiment mapped to the diploid genome (they are used to convert the normalised values back into read count values)
- expects a list of positions that carry an SNP between the two alleles of the specific hybrid in tsv formate from cactus (0-based)
- generates B73.[genotype].[condition].PF.GT.RN.csv files in the formate Chr(B73)  STop=POS(B73) REFal(B73)  ALTal ID([genotype]chr_[genotype]POS) Genotype (0/0 if [genotype] allele is B73, 1/1 if [genotype] allele is ALTal) normalised_Count(B73)   normalised_Count([genotype]) PEAK(B73) Peak([genotype])  BF (n.p. if no peak in either B73 or [genotype]) read_Count(B73) read_Count([genotype])



--------------------------
## Binomial testing of allele-specific binding:

For each hybrid, prepare a table that can be read into R with only sites that carry a SNP in that line, have at least one read on each allele and a count higher than 7 (adjust to a reasonable number, e.g. corresponding to 25 reads in one line, a cut off with read numbers can lead to lower numbers of allele-specific binding sites in lines with lower coverage).

```
for g in [genotype];do for tr in cond1 cond2; do gawk -v OFS='\t' -v g=$g 'BEGIN{print "NR","Chr","Pos","REFal","ALTal","PosPat","Genotype","EGcount_B73","EG_counts_Pat","Peak_B73","Peak_Pat","PF","Counts_B73","Counts_Pat"}{if($6=="1/1" && $0~"peak"){if($11>0 && $11<1){if($7>7 || $8>7){print NR,$0}}}}' dir/B73.${g}.${tr}.PF.GT.RN.csv > ${g}.${tr}.PF.GT.SNPs.CPM7.txt; done; done
```

```bash
for g in [genotype];do for tr in cond1 cond2; do Rscript --vanilla Binomial_fdr_SNPs_git.R ${g}.${tr}.PF.GT.SNPs.CPM7.txt ${g}.${tr}.bino.fdr.CPM7.txt; done; done
```
To get significant ones at FRD corrected p-value < 0.01:

```bash
for g in [genotype];do do for tr in cond1 cond2; do gawk -v OFS='\t' -v g=$g '{if(NR==1){print $0};if(NR>1 && $17<0.01){print $0}}' ${g}.${tr}.bino.fdr.CPM7.txt > ${g}.${tr}.bino.fdr.01.CPM7.txt; done;done
```
----------------------------------------------
## Controling for eventual mapping bias
To account for mapping bias, a control can be used to eliminate those significant sites that are also biased in control data, e.g. short-read sequencing genomic data. Reads of this data should be shortened to match the read length of the original experiment.

Counting of control reads is performed for the same positions as the experimental data, thus halliftover and translation of SNPs does not have to be repeated. The working directory should be the same as for the original analysis because the outputs of the previous steps are required.

--------------------------------
#get_control_counts.sh

These commands take the list of eligible sites for allele-specific binding (all that went into the binomial test, in our case all inside a peak, with a normalised read count larger than 7 and at least one read in each allele), creates a window in which the control data should be checked (e.g. 65 bp for MOA) and counts the mean normalised read count in this window for each of the two alleles. This is done for one condition only in our case (WW) to have the same threshold for both conditions. For cond2 (DS) we only count the control data on the allele-specific sites to use this data later to exclude sites with above threshold deviations from 1:1 in the control data.

Example code:

```bash
./get_control_counts.sh [genotype] [dir] [window] [cond1] [cond2] [SNPs_tested_bino_cond1] [genome_size_file] [sig_SNPs_cond2] [control_bedgraph]
```
- SNPs_tested_bino_cond1 sig_SNPs_cond2 control_bedgraph are expected to be in the format delivered by Binomial_fdr_SNPs_git.R



---------------------------------
## clean_AFPs_with_CNTRL.sh

Example command: 

```bash
./clean_AFPs_with_CNTRL.sh genotype dir Control_counts_B73_cond1 Control_counts_[genotype]_cond1 Binom_test_output_cond1 Cond1 Cond2 Control_counts_B73_cond2 Control_counts_[genotype]_cond2 Binom_test_output_cond2 outdir 
```

This script takes the control values from both alleles and calculates a binding frequency for each window (around an MP position for WW and AMP for DS). It then filters for non-AMP MPs and calculates the BF cut-off to exclude the upper and lower 5 percentile of biased binding. Thus, an upper and lower control BF cut-off is calculated. All AMPs with control BF above the upper and below the lower threshold are then filtered out.

- expects four bed files containing control counts for MP windows for cond1 for B73 and [genotype] and AMP windows for cond2.
- expects the binomial test output for both conditions
- expects names of the conditions (e.g. WW and DS)



