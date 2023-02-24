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
dir=/bgfs/karndt/Spt6_AID_Histone_ChIP_analysis
cd ${dir}

#MAKE DIRECTORIES FOR BIGWIGS
mkdir ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged
mkdir ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc

######################
## AVERAGED BIGWIGS ##
######################


# H2Bub
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_H2Bub_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_H2Bub_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_H2Bub.bw

bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_H2Bub_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_H2Bub_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_H2Bub.bw



# H3
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_H3_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_H3_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_H3.bw

bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_H3_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_H3_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_H3.bw


# H3K4me2
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_H3K4me2_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_H3K4me2_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_H3K4me2.bw

bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_H3K4me2_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_H3K4me2_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_H3K4me2.bw


# H3K4me3
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_H3K4me3_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_H3K4me3_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_H3K4me3.bw

bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_H3K4me3_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_H3K4me3_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_H3K4me3.bw


# Myc
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_Myc_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/DMSO_Myc_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_Myc.bw

bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_Myc_IP_Rep1_normalized.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/IAA_Myc_IP_Rep2_normalized.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation mean \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_Myc.bw


####################
## LOG2FC BIGWIGS ##
####################

# H2Bub
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_H2Bub.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_H2Bub.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation log2 \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc/H2Bub_log2fc.bw

# H3
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_H3.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_H3.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation log2 \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc/H3_log2fc.bw

# H3K4me2
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_H3K4me2.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_H3K4me2.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation log2 \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc/H3K4me2_log2fc.bw

# H3K4me3
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_H3K4me3.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_H3K4me3.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation log2 \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc/H3K4me3_log2fc.bw

# Myc
bigwigCompare \
--bigwig1 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/IAA_Myc.bw \
--bigwig2 ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/averaged/DMSO_Myc.bw \
--binSize 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation log2 \
--numberOfProcessors "max" \
--outFileFormat "bigwig" \
--outFileName ${dir}/Bowtie2_Alignments/Scer/normalized_bigwigs/log2fc/Myc_log2fc.bw


