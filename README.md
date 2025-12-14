# GBM-Surprisal
This is the repository for the Surprisal analysis and cell-cell signaling in GBM for the manuscript ""

## Data Preprocessing Part
Please refer to https://github.com/Shuozhen/zUMIs_for_Patho-DBiT for the raw fastq to the gene count matrix processes and zUMIs setup.

## Downstream Processing
1. We run each of the sample 
2. For the NICHES analysis, we run NICHES based on log normalized assay "Spatial" with k=9 as neighboring spots, but put the label as "NeighborhoodToCell_SCT_9"
3. The surprisal analysis is based on raw expression matrix and the scripts can be found here @https://github.com/ShouryoGhose/Surprisal-Analysis-of-Human-GBM-samples

