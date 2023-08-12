# Integrative analysis

[Integrative analysis](https://www.notion.so/Integrative-analysis-2016fdf0e5cf445db8709a776adbffd5?pvs=21)

The structural variant data is in `data/sv_sclc_ucologne_2015.tsv`.

1. For each gene involved in a structural variant (SV), determine the expression level of the gene in the sample that harbours the SV.

1. **Identify SVs that satisfy the following critiera:**
    - **The involved pair of genes both have elevated expression levels in samples with the SV compared to samples without the SV.**
    - **The second gene in the pair is in frame.**

```r
### The second gene in the pair is in frame.
identify_SVs <- structural_variant_data[,c(1,3,5,7,10)]
identify_SVs <- identify_SVs[identify_SVs$site2_effect_on_frame == "in-frame",]
genes_affected <- data.frame()
genes_affected <- cbind(identify_SVs$site1_hugo_symbol, identify_SVs$site2_hugo_symbol)
genes_affected <- as.data.frame(genes_affected)
colnames(genes_affected)[1] <- "gene_symbol"

transcriptomic_dgList <- DGEList(counts = transcriptomic_data)
transcriptomic_dgList <- calcNormFactors(transcriptomic_dgList, method="TMM")
transcriptomic_cpm <- cpm(transcriptomic_dgList, log=TRUE, prior.count=1, normalized.lib.sizes = TRUE)

### site 1 genes
transcriptomic_SV <- cbind(gene_symbol = rownames(transcriptomic_cpm), transcriptomic_cpm)
merged_site1_transcriptomic <- merge(genes_affected, transcriptomic_SV, by = "gene_symbol")
rownames(merged_site1_transcriptomic) <- merged_site1_transcriptomic$gene_symbol
merged_site1_transcriptomic <- merged_site1_transcriptomic[,-c(1:2)]
### site 2 genes
colnames(genes_affected)[2] <- "gene_symbol_2"
colnames(transcriptomic_SV)[1] <- "gene_symbol_2"
merged_site2_transcriptomic <- merge(genes_affected, transcriptomic_SV, by = "gene_symbol_2")
rownames(merged_site2_transcriptomic) <- merged_site2_transcriptomic$gene_symbol_2
merged_site2_transcriptomic <- merged_site2_transcriptomic[,-c(1:2)]

### TXNRD3 and LOC728588 not exist in transcriptomic data.
grep("TXNRD3", colnames(merged_site1_transcriptomic))
####integer(0)
grep("LOC728588", colnames(merged_site1_transcriptomic))
####integer(0)
merged_site2_transcriptomic <- merged_site2_transcriptomic[-c(7,12),]
identify_SVs <- identify_SVs[-c(6,7),]

### The involved pair of genes both have elevated expression levels
merged_transcriptomic <- rbind(merged_site1_transcriptomic, merged_site2_transcriptomic)
merged_transcriptomic_numeric <- as.data.frame(sapply(merged_transcriptomic, as.numeric))
rownames(merged_transcriptomic_numeric) <- rownames(merged_transcriptomic)
merged_transcriptomic_numeric <- merged_transcriptomic_numeric[,c(64,74,1:63,65:73,75:81)]
#### naming the samples as with SV and without SV
sample_col <- data.frame(sample = rep(c("with SV", "without SV"), c(2,79)))
row.names(sample_col) <- colnames(merged_transcriptomic_numeric)
library(pheatmap)
heatmap_SV <- pheatmap(merged_transcriptomic_numeric, scale = "row",
                       main = "Heatmap for each gene involved in a structural variant (SV)",
                       xlab = "Samples",
                       ylab = "Genes",
                       fontsize_col = 7, fontsize_row = 7, annotation_col = sample_col)
heatmap_SV
```

![sv heatmap](https://github.com/yerimeeei/exercise-r/assets/134043926/2267723e-bcb0-4e5f-9393-d3f3b8ae457f)

**EFHB-TYR Fusion**

- From the heatmap for each gene involved in a SV, pair of both EFHB and TYR genes both have elevated expression levels in sample with the SV compared to samples without the SV.
