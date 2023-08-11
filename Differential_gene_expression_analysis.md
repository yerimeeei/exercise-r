# Differential gene expression analysis

The clinical data is in `data/pheno_sclc_ucologne_2015.tsv`.

1. **Define two groups of tumours as early stage (stages I-II) vs. advanced stage tumours (stages III-IV), while excluding samples missing stage information.**

```r
tumor_stages <- data.frame(pheno_sclc_ucologne_2015[c(1,7)])
unique(tumor_stages$uicc_tumor_stage)

tumor_stages$stages <- NA
tumor_stages$stages[which(str_detect(tumor_stages$uicc_tumor_stage, "^I"))] <- "early stage"
tumor_stages$stages[which(str_detect(tumor_stages$uicc_tumor_stage, "^II"))] <- "early stage"
tumor_stages$stages[which(str_detect(tumor_stages$uicc_tumor_stage, "^III"))] <- "advanced stage"
tumor_stages$stages[which(str_detect(tumor_stages$uicc_tumor_stage, "^IV"))] <- "advanced stage"
tumor_stages <- drop_na(tumor_stages)
```

2. **Identify genes that differentially expressed in early vs. advanced stage tumours using an appropriate R package**

```r
library(DESeq2)
### when using the DESeq2 package, the count table should be unnormalized data.
rownames(tumor_stages) <- tumor_stages$patient_id
ncol(transcriptomic_data)
nrow(tumor_stages)
tumor_stages <- tumor_stages[rownames(tumor_stages) %in% colnames(transcriptomic_data),]
count_table <- transcriptomic_data[,colnames(transcriptomic_data) %in% rownames(tumor_stages)]
count_table <- as.matrix(count_table)
all(colnames(count_table) %in% rownames(tumor_stages))
all(colnames(count_table) == rownames(tumor_stages))

dds <- DESeqDataSetFromMatrix(countData = round(count_table),colData = tumor_stages, design = ~ stages)
dds <- DESeq(dds)
res <- results(dds, contrast = c("stages", "early stage", "advanced stage"))
res <- res[which(abs(res$log2FoldChange) > 0),]
res <- res[order(res$padj),]
## to find cut-off point
hist(res$padj, breaks=10)
```

![histogram_cut_1](https://github.com/yerimeeei/exercise-r/assets/134043926/dbc1dba1-2263-48ab-be3f-4ef1fb0255bd)

```r
sig_genes <- subset(res, padj < 0.8 & abs(log2FoldChange) > 1)
hist(sig_genes$padj, breaks=30)
sig_genes <- subset(res, padj < 0.7 & abs(log2FoldChange) > 1)
```

![histogram_cut_2](https://github.com/yerimeeei/exercise-r/assets/134043926/55639ec5-2bc2-4617-9a63-0dedbde7d706)

### G**enes that differentially expressed in early vs. advanced stage tumours**

: ANXA1, C8orf84, ECM1, AMIGO2
- Differentially expressed in early stage tumours: ECM1, AMIGO2
- Differentially expressed in advanced stage tumors: ANXA1, C8orf84

```r
library(EnhancedVolcano)
#selected = c("ANXA1", "C8orf84", "ECM1", "AMIGO2")
volcano_plot <- EnhancedVolcano(res, x = "log2FoldChange",y = "padj",
                                ylim = c(0,1), pCutoff = 0.7,
                                lab = row.names(res), #selectLab = selected,
                                title = 'Volcano Plot of early vs. advanced stage tumours')
volcano_plot
```

![volplot](https://github.com/yerimeeei/exercise-r/assets/134043926/a37ca5c6-aa56-4f1b-b3a9-82ec272076a0)

