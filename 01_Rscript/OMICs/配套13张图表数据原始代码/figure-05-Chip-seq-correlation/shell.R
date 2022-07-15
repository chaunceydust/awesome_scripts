# https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html
# https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html
multiBigwigSummary  bins -b  *WT*bw  -o wt_results.npz -p 8

plotCorrelation -in wt_results.npz  \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
--plotFileFormat pdf \
-o heatmap_SpearmanCorr_readCounts.pdf   \
--outFileCorMatrix SpearmanCorr_readCounts.tab


multiBigwigSummary  bins -b    *WT*bw  -o wt_merge_results.npz -p 8

plotCorrelation -in wt_merge_results.npz  \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
--plotFileFormat pdf \
-o merge_heatmap_SpearmanCorr_readCounts.pdf   \
--outFileCorMatrix merge_SpearmanCorr_readCounts.tab

