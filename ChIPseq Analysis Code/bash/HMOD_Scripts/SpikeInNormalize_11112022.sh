#!/bin/bash
# see http://core.sam.pitt.edu/node/5678
#SBATCH -N 1
#SBATCH -t 0-02:00:00 #2hrs
#SBATCH -c 8 #8cores
#SBATCH --mem=64g
#SBATCH --job-name=Normalize

##USE THIS TO OPEN AN ENVIRONMENT FOR TESTING
# srun -n1 -t02:00:00 -c 8 --mem=16g --pty bash

#LOAD MODULES
module load parallel/2017-05-22
module load deeptools/3.3.0

parallel --citation will cite


#SET DIRECTORIES FOR CURRENT ANALYSIS
#SET DIRECTORIES FOR CURRENT ANALYSIS
dir=/bgfs/karndt/Spt6_AID_Histone_ChIP_analysis
cd ${dir}

#MAKE DIRECTORY FOR BIGWIGS
mkdir ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs

#SET PATH VARIABLES
spike_in_bams=${dir}/Bowtie2_Alignments/Spomb/bam
spike_in_stats=${dir}/Bowtie2_Alignments/Spomb/stats

scer_bams=${dir}/Bowtie2_Alignments/Scer/bam
scer_stats=${dir}/Bowtie2_Alignments/Scer/stats


#GET READ COUNTS FOR SPIKEIN SAMPLES
# rm ${spike_in_stats}/names.txt ${spike_in_stats}/reads.txt

echo "SampleID" > ${spike_in_stats}/names.txt
echo "SpikeInReads" > ${spike_in_stats}/reads.txt

for i in `ls ${spike_in_stats}/*_markdups_stats.txt`
do

j=`echo $i | awk -F '[/]' '{print $8}' | awk -F '[_]' '{print $1 "_" $2}'`

echo ${j} >> ${spike_in_stats}/names.txt

head -n1 ${spike_in_stats}/${j}_markdups_stats.txt | cut -d " " -f 4 >> ${spike_in_stats}/reads.txt

done

#Combine 
paste ${spike_in_stats}/names.txt ${spike_in_stats}/reads.txt > ${spike_in_stats}/names_and_reads.txt



#GET READ COUNTS FOR SCER SAMPLES
rm ${scer_stats}/names.txt ${spike_in_stats}/reads.txt

echo "SampleID" > ${scer_stats}/names.txt
echo "Reads" > ${scer_stats}/reads.txt

for i in `ls ${scer_stats}/*_markdups_stats.txt`
do

j=`echo $i | awk -F '[/]' '{print $8}' | awk -F '[_]' '{print $1 "_" $2}'`

echo ${j} >> ${scer_stats}/names.txt

head -n1 ${scer_stats}/${j}_markdups_stats.txt | cut -d " " -f 4 >> ${scer_stats}/reads.txt

done


#Combine 
paste ${scer_stats}/names.txt ${scer_stats}/reads.txt > ${scer_stats}/names_and_reads.txt

#COMBINE ALL INFO
join -1 1 -2 1 ${dir}/Scripts/hmod_sample_info_11112022.txt ${scer_stats}/names_and_reads.txt > ${dir}/Scripts/updated_hmod_sample_info_11112022.tmp
join -1 1 -2 1 ${dir}/Scripts/updated_hmod_sample_info_11112022.tmp ${spike_in_stats}/names_and_reads.txt > ${dir}/Scripts/updated_hmod_sample_info_11112022.txt


#CALCULATE NORMALIZATION FACTOR
awk ' $4=="IP" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "_" $3  "_" $5}' ${dir}/Scripts/updated_hmod_sample_info_11112022.txt | sort -k 8 > ${dir}/Scripts/IP.tmp
awk ' $4=="input" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "_" $3  "_" $5}' ${dir}/Scripts/updated_hmod_sample_info_11112022.txt > ${dir}/Scripts/Input.tmp

#Modifications to reuse IP for H3K4me2
awk ' $4=="input" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "_" $3  "_" $5}' ${dir}/Scripts/updated_hmod_sample_info_11112022.txt \
| awk ' $3=="H3" {print}' | awk ' $3="H3K4me2" {print}' \
| awk ' $8="H3K4me2" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "_" $3  "_" $5}' > ${dir}/Scripts/Input.tmp2

# head ${dir}/Scripts/Input.tmp2 | awk '{print print  "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' 

#Modifications to reuse IP for H3K4me3
awk ' $4=="input" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "_" $3  "_" $5}' ${dir}/Scripts/updated_hmod_sample_info_11112022.txt \
| awk ' $3=="H3" {print}' | awk ' $3="H3K4me3" {print}' \
| awk ' $8="H3K4me3" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "_" $3  "_" $5}' > ${dir}/Scripts/Input.tmp3

#Modifications to reuse IP for Myc
awk ' $4=="input" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "_" $3  "_" $5}' ${dir}/Scripts/updated_hmod_sample_info_11112022.txt \
| awk ' $3=="H3" {print}' | awk ' $3="Myc" {print}' \
| awk ' $8="Myc" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "_" $3  "_" $5}' > ${dir}/Scripts/Input.tmp4

#combnine
cat ${dir}/Scripts/Input.tmp ${dir}/Scripts/Input.tmp2 ${dir}/Scripts/Input.tmp3 ${dir}/Scripts/Input.tmp4 | sort -k 8 > ${dir}/Scripts/Input.tmp5


join -1 8 -2 8 ${dir}/Scripts/IP.tmp ${dir}/Scripts/Input.tmp5 > ${dir}/Scripts/NormParalleleInfo.tmp
awk '{print $1 "\t" $2 "\t" $9 "\t" $7 "\t" $8 "\t" $14 "\t" $15 "\t" $3 "\t" $4 "\t" $5 "\t" $12 "\t" $6}' ${dir}/Scripts/NormParalleleInfo.tmp > ${dir}/Scripts/NormParalleleInfo.tmp2

#INFO ON CALCULATING NORMALIZATION FACTORS

# Jeronimo et al., 2019, Cell Reports 28, 1206–1218
# July 30, 2019 ª 2019 The Author(s).
# https://doi.org/10.1016/j.celrep.2019.06.097

# Reads for each biological replicate were independently aligned to the S. cerevisiae (UCSC sacCer3) and S. pombe (Downloaded from
# https://www.pombase.org/ on March 10th 2018) reference genomes using the short read aligner Bowtie 2 (version 2.2.4) (Langmead
# and Salzberg, 2012). Coverage for each base pair of the S. cerevisiae genome was computed using genomeCoverageBed from
# BEDTools (version 2.19.1) (Quinlan and Hall, 2010) and normalized using S. pombe reads as follow. The read density at each position
# of the S. cerevisiae genome was multiplied by a normalization factor N defined as:
# 
# 						N = 1,000,000 * P / R * C
# 
# where: 1,000,000 is an arbitrary chosen number used for convenience;
# 		R is the total number of reads in the IP sample that mapped to the S. pombe genome;
# 		C is the total number of reads in the Input sample that mapped to the S. cerevisiae genome;
# 		P is the total number of reads in the Input sample that mapped to the S. pombe genome.

#CALCULATE NORMALIZATION FACTORS
# Here
# R = $5 (IPSpReads)
# P = $7 (InputSpReads)
# C = $6 (InputScReads)
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" ((1000000 * $7) / ($5 * $6)) }' ${dir}/Scripts/NormParalleleInfo.tmp2 > ${dir}/Scripts/NormParalleleInfo.tmp3

#build header and add to file
echo "ChIPDescription_IPSampleID_InputSampleID_IPScReads_IPSpReads_InputScReads_InputSpReads_Treatment_ChIP_TypeIP_TypeInput_Rep_NormFactor" > header.tmp
awk -F "_" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13}' header.tmp > header.tmp2
cat header.tmp2 ${dir}/Scripts/NormParalleleInfo.tmp3 > ${dir}/Scripts/NormParalleleInfo.txt

#REMOVE TMP FILES
rm ${spike_in_stats}/names.txt ${spike_in_stats}/reads.txt ${scer_stats}/names.txt ${scer_stats}/reads.txt
rm ${dir}/Scripts/updated_hmod_sample_info_11112022.tmp
rm ${dir}/Scripts/IP.tmp ${dir}/Scripts/Input.tmp ${dir}/Scripts/Input.tmp2 ${dir}/Scripts/Input.tmp3 ${dir}/Scripts/Input.tmp4 ${dir}/Scripts/Input.tmp5
rm ${dir}/Scripts/NormParalleleInfo.tmp ${dir}/Scripts/NormParalleleInfo.tmp2 ${dir}/Scripts/NormParalleleInfo.tmp3
rm header.tmp header.tmp2


#SETUP FOR PARALLEL

#LOAD FILE
sample_index=${dir}/Scripts/NormParalleleInfo.txt

#print the file name of the master file to for the user
echo "Processing all samples in file: ${sample_index}"

#get lists for parallel
ChIPDescription=`awk '{print $1}' ${sample_index}`
echo $ChIPDescription
IPSampleID=`awk '{print $2}' ${sample_index}`
echo $InputSampleID
InputSampleID=`awk '{print $3}' ${sample_index}`
echo $IPSampleID
IPScReads=`awk '{print $4}' ${sample_index}`
echo $IPScReads
IPSpReads=`awk '{print $5}' ${sample_index}`
echo $IPSpReads
InputScReads=`awk '{print $6}' ${sample_index}`
echo $InputScReads
InputSpReads=`awk '{print $7}' ${sample_index}`
echo $InputSpReads
Treatment=`awk '{print $8}' ${sample_index}`
echo $Treatment
ChIP=`awk '{print $9}' ${sample_index}`
echo $ChIP
TypeIP=`awk '{print $10}' ${sample_index}`
echo $TypeIP
TypeInput=`awk '{print $11}' ${sample_index}`
echo $TypeInput
Rep=`awk '{print $12}' ${sample_index}`
echo $Rep
NormFactor=`awk '{print $13}' ${sample_index}`
echo $NormFactor



# GENERATE NORMALIZED BIGWIGS
parallel --header : bamCoverage \
--bam ${dir}/Bowtie2_Alignments/Scer/renamed_bams/{IPSampleID}_{Treatment}_{ChIP}_{TypeIP}_{Rep}.bam \
--scaleFactor {NormFactor} \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/{Treatment}_{ChIP}_{TypeIP}_{Rep}_normalized.bw \
--outFileFormat bigwig \
--binSize 1 \
--extendReads \
--numberOfProcessors "max" \
::: ${ChIPDescription} :::+ ${InputSampleID} :::+ ${IPSampleID} :::+ ${IPScReads} :::+ ${IPSpReads} :::+ ${InputScReads} :::+ ${InputSpReads} :::+ ${Treatment} :::+ ${ChIP} :::+ ${TypeIP} :::+ ${TypeInput} :::+ ${Rep} :::+ ${NormFactor}
 



