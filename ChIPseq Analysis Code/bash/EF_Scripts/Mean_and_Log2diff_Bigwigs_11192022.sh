#!/bin/bash
# see http://core.sam.pitt.edu/node/5678
#SBATCH -N 1
#SBATCH -t 0-02:00:00 #2hrs
#SBATCH -c 8 #8cores
#SBATCH --mem=64g
#SBATCH --job-name=BigWigs

##USE THIS TO OPEN AN ENVIRONMENT FOR TESTING
# srun -n1 -t02:00:00 -c 8 --mem=16g --pty bash


#LOAD MODULES
module load deeptools/3.3.0

#SET DIRECTORIES FOR CURRENT ANALYSIS
dir=/bgfs/karndt/Spt6_AID_ChIPSeq_analysis
cd ${dir}

#MAKE DIRECTORIES FOR BIGWIGS
mkdir ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged
mkdir ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc

######################
## AVERAGED BIGWIGS ##
######################

# 8WG16
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_8WG16_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_8WG16_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_8WG16.bw

bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_8WG16_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_8WG16_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_8WG16.bw

# 
# # HSV
# bigwigCompare \
# --bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_HSV_IP_Rep1_normalized.bw \
# --bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_HSV_IP_Rep2_normalized.bw \
# --binSize 1 \
# --skipNonCoveredRegions \
# --skipZeroOverZero \
# --operation mean \
# --numberOfProcessors "max" \
# --outFileFormat "bigwig" \
# --outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_HSV.bw
# 
# bigwigCompare \
# --bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_HSV_IP_Rep1_normalized.bw \
# --bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_HSV_IP_Rep2_normalized.bw \
# --binSize 1 \
# --skipNonCoveredRegions \
# --skipZeroOverZero \
# --operation mean \
# --numberOfProcessors "max" \
# --outFileFormat "bigwig" \
# --outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_HSV.bw
# 
# 
# # Spt5
# bigwigCompare \
# --bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_Spt5_IP_Rep1_normalized.bw \
# --bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_Spt5_IP_Rep2_normalized.bw \
# --binSize 1 \
# --skipNonCoveredRegions \
# --skipZeroOverZero \
# --operation mean \
# --numberOfProcessors "max" \
# --outFileFormat "bigwig" \
# --outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_Spt5.bw
# 
# bigwigCompare \
# --bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_Spt5_IP_Rep1_normalized.bw \
# --bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_Spt5_IP_Rep2_normalized.bw \
# --binSize 1 \
# --skipNonCoveredRegions \
# --skipZeroOverZero \
# --operation mean \
# --numberOfProcessors "max" \
# --outFileFormat "bigwig" \
# --outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_Spt5.bw
# 
# 
# # V5
# bigwigCompare \
# --bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_V5_IP_Rep1_normalized.bw \
# --bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_V5_IP_Rep2_normalized.bw \
# --binSize 1 \
# --skipNonCoveredRegions \
# --skipZeroOverZero \
# --operation mean \
# --numberOfProcessors "max" \
# --outFileFormat "bigwig" \
# --outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_V5.bw
# 
# bigwigCompare \
# --bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_V5_IP_Rep1_normalized.bw \
# --bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_V5_IP_Rep2_normalized.bw \
# --binSize 1 \
# --skipNonCoveredRegions \
# --skipZeroOverZero \
# --operation mean \
# --numberOfProcessors "max" \
# --outFileFormat "bigwig" \
# --outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_V5.bw


####################
## LOG2FC BIGWIGS ##
####################

# 8WG16
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_8WG16.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_8WG16.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation log2 \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc/8WG16_log2fc.bw

# # HSV
# bigwigCompare \
# --bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_HSV.bw \
# --bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_HSV.bw \
# --binSize 1 \
# --skipNonCoveredRegions \
# --skipZeroOverZero \
# --operation log2 \
# --numberOfProcessors "max" \
# --outFileFormat "bigwig" \
# --outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc/HSV_log2fc.bw
# 
# # Spt5
# bigwigCompare \
# --bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_Spt5.bw \
# --bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_Spt5.bw \
# --binSize 1 \
# --skipNonCoveredRegions \
# --skipZeroOverZero \
# --operation log2 \
# --numberOfProcessors "max" \
# --outFileFormat "bigwig" \
# --outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc/Spt5_log2fc.bw
# 
# # V5
# bigwigCompare \
# --bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_V5.bw \
# --bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_V5.bw \
# --binSize 1 \
# --skipNonCoveredRegions \
# --skipZeroOverZero \
# --operation log2 \
# --numberOfProcessors "max" \
# --outFileFormat "bigwig" \
# --outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc/V5_log2fc.bw
# 
# 
