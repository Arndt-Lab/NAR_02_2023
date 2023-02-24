#!/bin/bash

#Set directory
dir=/Volumes/Mitch_Home_Backup/Arndt_Lab_Consulting/ChIPseq_Redo_Spt6_AID_112022/DeepTools_Analysis_ChIPseq_Redo_11192022/
cd ${dir}


##############
## FIGURE 5 ##
##############


# METAPLOTS

#8WG16
computeMatrix scale-regions \
-b 500 -a 500 \
-S ${dir}/averaged_bigwigs/DMSO_8WG16.bw \
${dir}/averaged_bigwigs/IAA_8WG16.bw  \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
--binSize 25 \
--missingDataAsZero \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/Spt6_AID_8wg16_matrix.gz

plotProfile -m ${dir}/analysis/Spt6_AID_8wg16_matrix.gz \
-out ${dir}/analysis/Spt6_AID_8wg16_matrix.png \
--dpi 300 \
--plotType=lines \
--perGroup \
--colors "#000000" "#990000" \
--samplesLabel "DMSO 8WG16" "IAA 8WG16" \
--legendLocation "upper-left" \
--averageType "mean" \
--plotHeight 8 \
--plotWidth 12 \
--yMin 0 \
--yMax 100 \
--regionsLabel "8WG16" \
--plotFileFormat "png"



#Spt5
computeMatrix scale-regions \
-b 500 -a 500 \
-S ${dir}/averaged_bigwigs/DMSO_Spt5.bw \
${dir}/averaged_bigwigs/IAA_Spt5.bw  \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
--binSize 25 \
--missingDataAsZero \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/Spt6_AID_Spt5_matrix.gz

plotProfile -m ${dir}/analysis/Spt6_AID_Spt5_matrix.gz \
-out ${dir}/analysis/Spt6_AID_Spt5_matrix.png \
--dpi 300 \
--plotType=lines \
--perGroup \
--colors "#000000" "#990000" \
--samplesLabel "DMSO Spt5" "IAA Spt5" \
--legendLocation "upper-left" \
--averageType "mean" \
--plotHeight 8 \
--plotWidth 12 \
--yMin 0 \
--yMax 600 \
--regionsLabel "Spt5" \
--plotFileFormat "png"


#V5
computeMatrix scale-regions \
-b 500 -a 500 \
-S ${dir}/averaged_bigwigs/DMSO_V5.bw \
${dir}/averaged_bigwigs/IAA_V5.bw  \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
--binSize 25 \
--missingDataAsZero \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/Spt6_AID_V5_matrix.gz

plotProfile -m ${dir}/analysis/Spt6_AID_V5_matrix.gz \
-out ${dir}/analysis/Spt6_AID_V5_matrix.png \
--dpi 300 \
--plotType=lines \
--perGroup \
--colors "#000000" "#990000" \
--samplesLabel "DMSO V5" "IAA V5" \
--legendLocation "upper-left" \
--averageType "mean" \
--plotHeight 8 \
--plotWidth 12 \
--yMin 0 \
--yMax 50 \
--regionsLabel "V5" \
--plotFileFormat "png"


#HSV
computeMatrix scale-regions \
-b 500 -a 500 \
-S ${dir}/averaged_bigwigs/DMSO_HSV.bw \
${dir}/averaged_bigwigs/IAA_HSV.bw  \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
--binSize 25 \
--missingDataAsZero \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/Spt6_AID_HSV_matrix.gz

plotProfile -m ${dir}/analysis/Spt6_AID_HSV_matrix.gz \
-out ${dir}/analysis/Spt6_AID_HSV_matrix.png \
--dpi 300 \
--plotType=lines \
--perGroup \
--colors "#000000" "#990000" \
--samplesLabel "DMSO HSV" "IAA HSV" \
--legendLocation "upper-left" \
--averageType "mean" \
--plotHeight 8 \
--plotWidth 12 \
--yMin 0 \
--yMax 200 \
--regionsLabel "HSV" \
--plotFileFormat "png"


#Myc
computeMatrix scale-regions \
-b 500 -a 500 \
-S ${dir}/averaged_bigwigs/DMSO_Myc.bw \
${dir}/averaged_bigwigs/IAA_Myc.bw  \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
--binSize 25 \
--missingDataAsZero \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/Spt6_AID_Myc_matrix.gz

plotProfile -m ${dir}/analysis/Spt6_AID_Myc_matrix.gz \
-out ${dir}/analysis/Spt6_AID_Myc_matrix.png \
--dpi 300 \
--plotType=lines \
--perGroup \
--colors "#000000" "#990000" \
--samplesLabel "DMSO Myc" "IAA Myc" \
--legendLocation "upper-left" \
--averageType "mean" \
--plotHeight 8 \
--plotWidth 12 \
--yMin 0 \
--yMax 200 \
--regionsLabel "Myc" \
--plotFileFormat "png"


# HEATMAP
file_length=`wc -l ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all.bed  | awk '{print $1}'`
heatmap_height=`expr "$file_length" / 100`

computeMatrix reference-point \
--referencePoint "TSS" -b 500 -a 5000 \
-S ${dir}/log2fc_bigwigs/8WG16_log2fc.bw \
${dir}/log2fc_bigwigs/Spt5_log2fc.bw \
${dir}/log2fc_bigwigs/V5_log2fc.bw \
${dir}/log2fc_bigwigs/HSV_log2fc.bw \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all.bed \
--binSize 25 \
--missingDataAsZero \
--sortUsing region_length \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/EF_Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all_matrix.gz


plotHeatmap -m ${dir}/analysis/EF_Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all_matrix.gz \
-out ${dir}/analysis/EF_Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all_log2fc2.png \
--dpi 300 \
--sortUsing region_length \
--missingDataColor 1 \
--sortRegions ascend \
--colorMap 'seismic' \
--zMin -5 \
--zMax 5 \
--whatToShow 'plot, heatmap and colorbar' \
--heatmapHeight "$heatmap_height" \
--heatmapWidth 12


# HEATMAP
file_length=`wc -l ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all.bed  | awk '{print $1}'`
heatmap_height=`expr "$file_length" / 100`

computeMatrix reference-point \
--referencePoint "TSS" -b 500 -a 5000 \
-S ${dir}/averaged_bigwigs/DMSO_8WG16.bw \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all.bed \
--binSize 25 \
--missingDataAsZero \
--sortUsing region_length \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/EF_Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all_matrix_test.gz


plotHeatmap -m ${dir}/analysis/EF_Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all_matrix_test.gz \
-out ${dir}/analysis/EF_Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all_log2fc2_test.png \
--dpi 300 \
--sortUsing region_length \
--missingDataColor 1 \
--sortRegions ascend \
--colorMap 'Greys' \
--zMin 0 \
--zMax 50 \
--whatToShow 'plot, heatmap and colorbar' \
--heatmapHeight "$heatmap_height" \
--heatmapWidth 5


# CORRELATION HEATMAP


multiBigwigSummary BED-file \
	--BED ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all.bed \
	--bwfiles ${dir}/log2fc_bigwigs/8WG16_log2fc.bw \
	${dir}/log2fc_bigwigs/Spt5_log2fc.bw \
	${dir}/log2fc_bigwigs/V5_log2fc.bw \
	${dir}/log2fc_bigwigs/HSV_log2fc.bw \
	--numberOfProcessors "max" \
	--outRawCounts ${dir}/analysis/multiBigwigSummary_log2fc_AID.txt \
	--outFileName ${dir}/analysis/multiBigwigSummary_log2fc_AID.npz 


plotCorrelation \
	--corData ${dir}/analysis/multiBigwigSummary_log2fc_AID.npz \
	--corMethod 'spearman' \
	--whatToPlot 'heatmap' \
	--plotNumbers \
	--skipZeros \
	--plotHeight 8 \
	--plotWidth 8 \
	--outFileCorMatrix ${dir}/analysis/multiBigwigSummary_log2fc_spearmans_correlation_heatmap_AID.tsv \
	--plotFile ${dir}/analysis/multiBigwigSummary_log2fc_spearmans_correlation_heatmap_AID.png



##############
## FIGURE 6 ##
##############

# METAPLOTS

#H3
computeMatrix scale-regions \
-b 500 -a 500 \
-S ${dir}/averaged_bigwigs/DMSO_H3.bw \
${dir}/averaged_bigwigs/IAA_H3.bw  \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
--binSize 25 \
--missingDataAsZero \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/Spt6_AID_H3_matrix.gz

plotProfile -m ${dir}/analysis/Spt6_AID_H3_matrix.gz \
-out ${dir}/analysis/Spt6_AID_H3_matrix.png \
--dpi 300 \
--plotType=lines \
--perGroup \
--colors "#000000" "#990000" \
--samplesLabel "DMSO H3" "IAA H3" \
--legendLocation "upper-left" \
--averageType "mean" \
--plotHeight 8 \
--plotWidth 12 \
--yMin 0 \
--yMax 30 \
--regionsLabel "H3" \
--plotFileFormat "png"



#H2Bub
computeMatrix scale-regions \
-b 500 -a 500 \
-S ${dir}/averaged_bigwigs/DMSO_H2Bub.bw \
${dir}/averaged_bigwigs/IAA_H2Bub.bw  \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
--binSize 25 \
--missingDataAsZero \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/Spt6_AID_H2Bub_matrix.gz

plotProfile -m ${dir}/analysis/Spt6_AID_H2Bub_matrix.gz \
-out ${dir}/analysis/Spt6_AID_H2Bub_matrix.png \
--dpi 300 \
--plotType=lines \
--perGroup \
--colors "#000000" "#990000" \
--samplesLabel "DMSO H2Bub" "IAA H2Bub" \
--legendLocation "upper-left" \
--averageType "mean" \
--plotHeight 8 \
--plotWidth 12 \
--yMin 0 \
--yMax 100 \
--regionsLabel "H2Bub" \
--plotFileFormat "png"


#H3K4me2
computeMatrix scale-regions \
-b 500 -a 500 \
-S ${dir}/averaged_bigwigs/DMSO_H3K4me2.bw \
${dir}/averaged_bigwigs/IAA_H3K4me2.bw \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
--binSize 25 \
--missingDataAsZero \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/Spt6_AID_H3K4me2_matrix.gz

plotProfile -m ${dir}/analysis/Spt6_AID_H3K4me2_matrix.gz \
-out ${dir}/analysis/Spt6_AID_H3K4me2_matrix.png \
--dpi 300 \
--plotType=lines \
--perGroup \
--colors "#000000" "#990000" \
--samplesLabel "DMSO H3K4me2" "IAA H3K4me2" \
--legendLocation "upper-left" \
--averageType "mean" \
--plotHeight 8 \
--plotWidth 12 \
--yMin 0 \
--yMax 150 \
--regionsLabel "H3K4me2" \
--plotFileFormat "png"


#H3K4me3
computeMatrix scale-regions \
-b 500 -a 500 \
-S ${dir}/averaged_bigwigs/DMSO_H3K4me3.bw \
${dir}/averaged_bigwigs/IAA_H3K4me3.bw  \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
--binSize 25 \
--missingDataAsZero \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/Spt6_AID_H3K4me3_matrix.gz

plotProfile -m ${dir}/analysis/Spt6_AID_H3K4me3_matrix.gz \
-out ${dir}/analysis/Spt6_AID_H3K4me3_matrix.png \
--dpi 300 \
--plotType=lines \
--perGroup \
--colors "#000000" "#990000" \
--samplesLabel "DMSO H3K4me3" "IAA H3K4me3" \
--legendLocation "upper-left" \
--averageType "mean" \
--plotHeight 8 \
--plotWidth 12 \
--yMin 0 \
--yMax 200 \
--regionsLabel "H3K4me3" \
--plotFileFormat "png"


# HEATMAP
file_length=`wc -l ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all.bed  | awk '{print $1}'`
heatmap_height=`expr "$file_length" / 100`

computeMatrix reference-point \
--referencePoint "TSS" -b 500 -a 5000 \
-S ${dir}/log2fc_bigwigs/H2Bub_log2fc.bw \
${dir}/log2fc_bigwigs/H3_log2fc.bw \
${dir}/log2fc_bigwigs/H3K4me2_log2fc.bw \
${dir}/log2fc_bigwigs/H3K4me3_log2fc.bw \
-R ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all.bed \
--binSize 25 \
--missingDataAsZero \
--sortUsing region_length \
--averageTypeBins mean \
--numberOfProcessors "max" \
-out ${dir}/analysis/HMOD_Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all_matrix.gz


plotHeatmap -m ${dir}/analysis/HMOD_Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all_matrix.gz \
-out ${dir}/analysis/HMOD_Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all_log2fc2.png \
--dpi 300 \
--sortUsing region_length \
--missingDataColor 1 \
--sortRegions ascend \
--colorMap 'seismic' \
--zMin -5 \
--zMax 5 \
--whatToShow 'plot, heatmap and colorbar' \
--heatmapHeight "$heatmap_height" \
--heatmapWidth 12









###############
## FIGURE S8 ##
###############

# CORRELATION HEATMAP


multiBigwigSummary BED-file \
	--BED ${dir}/analysis/Scer_transcripts_w_verifiedORFs-nonoverlapping_decending_quartiles_all.bed \
	--bwfiles ${dir}/log2fc_bigwigs/8WG16_log2fc.bw \
	${dir}/log2fc_bigwigs/Spt5_log2fc.bw \
	${dir}/log2fc_bigwigs/V5_log2fc.bw \
	${dir}/log2fc_bigwigs/HSV_log2fc.bw \
	${dir}/log2fc_bigwigs/H3_log2fc.bw \
	${dir}/log2fc_bigwigs/H2Bub_log2fc.bw \
	${dir}/log2fc_bigwigs/H3K4me2_log2fc.bw \
	${dir}/log2fc_bigwigs/H3K4me3_log2fc.bw \
	--numberOfProcessors "max" \
	--outRawCounts ${dir}/analysis/multiBigwigSummary_log2fc_AID_all.txt \
	--outFileName ${dir}/analysis/multiBigwigSummary_log2fc_AID_all.npz 


plotCorrelation \
	--corData ${dir}/analysis/multiBigwigSummary_log2fc_AID_all.npz \
	--corMethod 'spearman' \
	--whatToPlot 'heatmap' \
	--plotNumbers \
	--skipZeros \
	--plotHeight 8 \
	--plotWidth 8 \
	--outFileCorMatrix ${dir}/analysis/multiBigwigSummary_log2fc_spearmans_correlation_heatmap_AID_all.tsv \
	--plotFile ${dir}/analysis/multiBigwigSummary_log2fc_spearmans_correlation_heatmap_AID_all.png




