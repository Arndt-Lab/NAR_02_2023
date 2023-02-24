#!/bin/bash
# see http://core.sam.pitt.edu/node/5678
#SBATCH -N 1
#SBATCH -t 0-00:20:00 #2days
#SBATCH -c 8 #8cores
#SBATCH --mem=16g
#SBATCH --job-name=Rename

##USE THIS TO OPEN AN ENVIRONMENT FOR TESTING
# srun -n1 -t02:00:00 -c 8 --mem=16g --pty bash

#LOAD MODULES
module load parallel/2017-05-22
parallel --citation will cite

#SET DIRECTORIES FOR CURRENT ANALYSIS
dir=/bgfs/karndt/Spt6_AID_Histone_ChIP_analysis
cd ${dir}

#SETUP FOR PARALLEL
sample_index=${dir}/Scripts/hmod_sample_info_11112022.txt

#print the file name of the master file to for the user
echo "Processing all samples in file: ${sample_index}"

#get lists for parallel
#get file prefix list
SampleID=`awk '{print $1}' ${sample_index}`
#get stain list
Treatment=`awk '{print $2}' ${sample_index}`
#get protein list
ChIP=`awk '{print $3}' ${sample_index}`
#get rep listfile_prefix = awk '{print $1}' $1
Type=`awk '{print $4}' ${sample_index}`
#get list of qubit readings
Rep=`awk '{print $5}' ${sample_index}`

mkdir ${dir}/Bowtie2_Alignments/Scer/renamed_bams


#CREATE SYMBOLIC LINKS FOR BAMS RENAMED BASED ON SAMPLE INFORMATION LIST
parallel --header : ln -s ${dir}/Bowtie2_Alignments/Scer/bam/{SampleID}_nodups.bam \
	${dir}/Bowtie2_Alignments/Scer/renamed_bams/{SampleID}_{Treatment}_{ChIP}_{Type}_{Rep}.bam \
	::: ${SampleID} :::+ ${Treatment} :::+ ${ChIP} :::+ ${Type} :::+ ${Rep}

#CREATE SYMBOLIC LINKS FOR BAIS RENAMED BASED ON SAMPLE INFORMATION LIST
parallel --header : ln -s ${dir}/Bowtie2_Alignments/Scer/bam/{SampleID}_nodups.bam.bai \
	${dir}/Bowtie2_Alignments/Scer/renamed_bams/{SampleID}_{Treatment}_{ChIP}_{Type}_{Rep}.bam.bai \
	::: ${SampleID} :::+ ${Treatment} :::+ ${ChIP} :::+ ${Type} :::+ ${Rep}

#VERIFY SYMBOLIC LINKS
ls -l ${dir}/Bowtie2_Alignments/Scer/renamed_bams/*


