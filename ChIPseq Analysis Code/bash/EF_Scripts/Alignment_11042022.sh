#!/bin/bash
# see http://core.sam.pitt.edu/node/5678
#SBATCH -N 1
#SBATCH -t 2-00:00:00 #2days
#SBATCH -c 8 #8cores
#SBATCH --mem=128g
#SBATCH --job-name=Align

#SET DIRECTORIES FOR CURRENT ANALYSIS
fastq_dir=/bgfs/karndt/shared/karndt_rfe4_wam30/Arndt_Spt6_AID_10_4_2022
dir=/bgfs/karndt/Spt6_AID_ChIPSeq_analysis

##USE THIS TO OPEN AN ENVIRONMENT FOR TESTING
# srun -n1 -t06:00:00 -c 8 --mem=128g --pty bash

##LOAD MODULES
module load bowtie2/2.4.5
module load gcc/8.2.0
module load samtools/1.9

#MAKE SCRIPT STOP ON ANY ERROR
# set -ue

#SET DIRECTORY
cd ${dir}

#MAKE OUTPUT DIRECTORYS FOR ALIGNMENTS
mkdir ${dir}/Bowtie2_Alignments

mkdir ${dir}/Bowtie2_Alignments/Spomb
mkdir ${dir}/Bowtie2_Alignments/Spomb/bam
mkdir ${dir}/Bowtie2_Alignments/Spomb/stats

mkdir ${dir}/Bowtie2_Alignments/Scer
mkdir ${dir}/Bowtie2_Alignments/Scer/FASTQ
mkdir ${dir}/Bowtie2_Alignments/Scer/bam
mkdir ${dir}/Bowtie2_Alignments/Scer/stats


#SET PATH VARIABLES
spike_in_index=/bgfs/karndt/Annotations/Schizosaccharomyces_pombe_Ensembl_EF2/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/Bowtie2Index/genome
sac_cer_index=/bgfs/karndt/Annotations/Saccharomyces_cerevisiae_Ensembl_R64-1-1/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/Bowtie2Index/genome

spike_in_bams=${dir}/Bowtie2_Alignments/Spomb/bam
spike_in_stats=${dir}/Bowtie2_Alignments/Spomb/stats

scer_fastq=${dir}/Bowtie2_Alignments/Scer/FASTQ
scer_bams=${dir}/Bowtie2_Alignments/Scer/bam
scer_stats=${dir}/Bowtie2_Alignments/Scer/stats

#START LOOP
for i in `ls ${fastq_dir} | awk -F "_" '{print $1 "_" $2}'`
do


#ALIGN TO SPIKEIN
bowtie2 \
	-p 8 \
	-q \
	-k 2 \
	--no-unal \
	--no-mixed \
	--no-discordant \
	--no-overlap \
	-x ${spike_in_index} \
	-1 ${fastq_dir}/${i}_R1_001.fastq.gz \
	-2 ${fastq_dir}/${i}_R2_001.fastq.gz \
	--un-conc-gz ${scer_fastq}/${i}_no_aln_to_spkin_R%.fastq.gz \
	--met-file ${spike_in_stats}/${i}_all_metrics_spkin.txt \
	2> ${spike_in_stats}/${i}_all_stats_spkin.txt > ${spike_in_bams}/${i}_spkin.sam

#CALCULATE INITIAL BAM STATISTICS
samtools flagstat -@8 ${spike_in_bams}/${i}_spkin.sam > ${spike_in_stats}/${i}_initial_sam_stats_spkin.txt

#SORT FOR FIXMATE
samtools sort -n -o ${spike_in_bams}/${i}_spkin_sorted.bam -O BAM -T ${i}.temp ${spike_in_bams}/${i}_spkin.sam

#REMOVE ORIGINAL ALIGNMENT FILE
rm ${spike_in_bams}/${i}_spkin.sam

#MARK PCR DUPLICATES
samtools fixmate -r -m -O BAM -@8 ${spike_in_bams}/${i}_spkin_sorted.bam ${spike_in_bams}/${i}_spkin_fixmate.bam

#RESORT BAM FOR MARKDUP
samtools sort -o ${spike_in_bams}/${i}_spkin_sorted2.bam -O BAM -T ${i}.temp ${spike_in_bams}/${i}_spkin_fixmate.bam

#REMOVE PCR DUPLICATES
samtools markdup -r -s -T nodups -O BAM -@8 ${spike_in_bams}/${i}_spkin_sorted2.bam ${spike_in_bams}/${i}_spkin_nodups.bam 2> ${spike_in_stats}/${i}_markdups_stats.txt

#REMOVE EXTRA BAM FILES
rm ${spike_in_bams}/${i}_spkin_sorted.bam ${spike_in_bams}/${i}_spkin_fixmate.bam ${spike_in_bams}/${i}_spkin_sorted2.bam

#CALCULATE FINAL BAM STATISTICS
samtools flagstat -@8 ${spike_in_bams}/${i}_spkin_nodups.bam > ${spike_in_stats}/${i}_final_bam_stats_spkin.txt



	#\
# 	| samtools sort -n -O BAM -@8 -T ${i}.temp - \
# 	| samtools fixmate - -r -O BAM -@8 - \
# 	| samtools sort -c -O BAM -@8 -T ${i}.temp - \ 
	
	
# 	\
# 	| samtools markdup -r -s -T nodups 

# ERROR The input file must be coordinate sorted and must have gone through fixmates with the mate scoring option o
# n.

#ALIGN TO SCER
bowtie2 \
	-p 8 \
	-q \
	-k 2 \
	--no-unal \
	--no-mixed \
	--no-discordant \
	--no-overlap \
	-x ${sac_cer_index} \
	-1 ${scer_fastq}/${i}_no_aln_to_spkin_R1.fastq.gz \
	-2 ${scer_fastq}/${i}_no_aln_to_spkin_R2.fastq.gz \
	--met-file ${scer_stats}/${i}_all_metrics.txt \
	2> ${scer_stats}/${i}_all_stats.txt > ${scer_bams}/${i}.sam
	
#CALCULATE INITIAL BAM STATISTICS
samtools flagstat -@8 ${scer_bams}/${i}.sam > ${scer_stats}/${i}_initial_sam_stats.txt

#SORT FOR FIXMATE
samtools sort -n -o ${scer_bams}/${i}_sorted.bam -O BAM -T ${i}.temp ${scer_bams}/${i}.sam

#REMOVE ORIGINAL ALIGNMENT FILE
rm ${scer_bams}/${i}.sam

#MARK PCR DUPLICATES
samtools fixmate -r -m -O BAM -@8 ${scer_bams}/${i}_sorted.bam ${scer_bams}/${i}_fixmate.bam

#RESORT BAM FOR MARKDUP
samtools sort -o ${scer_bams}/${i}_sorted2.bam -O BAM -T ${i}.temp ${scer_bams}/${i}_fixmate.bam

#REMOVE PCR DUPLICATES
samtools markdup -r -s -T nodups -O BAM -@8 ${scer_bams}/${i}_sorted2.bam ${scer_bams}/${i}_nodups.bam 2> ${scer_stats}/${i}_markdups_stats.txt

#REMOVE EXTRA BAM FILES
rm ${scer_bams}/${i}_sorted.bam ${scer_bams}/${i}_fixmate.bam ${scer_bams}/${i}_sorted2.bam

#CALCULATE FINAL BAM STATISTICS
samtools flagstat -@8 ${scer_bams}/${i}_nodups.bam > ${scer_stats}/${i}_final_bam_stats.txt

#INDEX BAM
samtools index ${scer_bams}/${i}_nodups.bam

done


exit
