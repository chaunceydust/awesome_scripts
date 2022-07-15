# https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html
# https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html
multiBigwigSummary  bins -b  *.bw -o results.npz -p 4
plotCorrelation -in results.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
--plotFileFormat pdf \
-o heatmap_SpearmanCorr_readCounts.pdf   \
--outFileCorMatrix SpearmanCorr_readCounts.tab