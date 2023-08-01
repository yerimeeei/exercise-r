# Programming exercise for R

## Mutation data analysis

[mutation data analysis](https://www.notion.so/mutation-data-analysis-6bbe1846a5f744c292163a47c02a2f44?pvs=21)

The mutation data is in `data/mutations_sclc_ucologne_2015.rds`.

1. **Identify the top 10 most frequently mutated genes. Identify samples whose mutation count is in the 80 to 90 percentile.**

```r
# Mutation Data Analysis
mutation_data <- readRDS("/Users/rimi/bioinformatics/DV/exercise-r/data/mutations_sclc_ucologne_2015.rds")
## 1.1 Identify the top 10 most frequently mutated genes.
library(dplyr)
frequency_gene <- mutation_data %>% group_by(mutation_data$gene) %>% summarise(freq=n())
frequency_gene$ranks <- rank(frequency_gene$freq)
frequency_gene <- frequency_gene[order(frequency_gene$ranks, decreasing = TRUE),]
frequency_gene <- frequency_gene[-c(3)]
### TTN, TP53, MUC16, RB1, RYR2, LRP1B, ZFHX4, CSMD3, USH2A, NAV3 = SYNE1

## 1.2 Identify samples whose mutation count is in the 80 to 90 percentile.
library(dplyr)
length(unique(mutation_data$sample_id))
frequency_sample <- mutation_data %>% group_by(mutation_data$sample_id) %>% summarise(freq=n())
frequency_sample$ranks <- rank(frequency_sample$freq)
frequency_sample <- frequency_sample[order(frequency_sample$ranks, decreasing = TRUE),]
colnames(frequency_sample)[1] <- "sample_id"frequency_sample <- frequency_sample[-c(3)]
quantile(frequency_sample$freq, c(0.8,0.9))
samples <- frequency_sample[frequency_sample$freq >= quantile(frequency_sample$freq, c(0.8)) & frequency_sample$freq <= quantile(frequency_sample$freq, c(0.9)),]
samples$sample_id
```

<img width="248" alt="Screenshot_2023-07-30_at_2 44 56_PM" src="https://github.com/yerimeeei/exercise-r/assets/134043926/dda0b39a-c057-4132-ab50-afcb40aab3f5">

- TTN, TP53, MUC16, RB1, RYR2, LRP1B, ZFHX4, CSMD3, USH2A, NAV3 = SYNE1
- "sclc_ucologne_2015_S01020" "sclc_ucologne_2015_S02344" "sclc_ucologne_2015_S02285"
"sclc_ucologne_2015_S01022" "sclc_ucologne_2015_S02376" "sclc_ucologne_2015_S00841"
"sclc_ucologne_2015_S01023" "sclc_ucologne_2015_S02248" "sclc_ucologne_2015_S02328"
"sclc_ucologne_2015_S01861" "sclc_ucologne_2015_S02242" "sclc_ucologne_2015_S02384”
2. **Categorize variants based on their expected effects using the `data/mutation_effects.tsv` table. Generate a count matrix containing the numbers of loss-of-function and neutral mutations for each gene.**

```r
mutation_effects <- read.delim("~/bioinformatics/DV/exercise-r/data/mutation_effects.tsv")
unique(mutation_data$variant_class)
### "Missense_Mutation"      "Silent"                 "Translation_Start_Site"
### "Nonsense_Mutation"      "Splice_Site"            "Frame_Shift_Del"       
### "Intron"                 "In_Frame_Del"           "Splice_Region"         
### "Frame_Shift_Ins"        "RNA"                    "Nonstop_Mutation"      
### "3'Flank"                "5'UTR"                  "5'Flank"               
### "In_Frame_Ins"           "3'UTR"                  "IGR"                   
### "frame_shift_del"
variants <- mutation_data[c(1,6)]
variants_effect <- merge(variants, mutation_effects, by = "variant_class")
variants_effect <- variants_effect[c(2,1,3)]

colnames(frequency_gene)[1] <- "gene"
counts <- merge(frequency_gene, variants_effect, by = "gene")
counts <- counts[!counts$effect == "uncertain",]
counts <- counts[-c(2,3)]
counts <- counts %>% group_by(gene, effect) %>% summarise(freq=n())
count_neutral <- counts[counts$effect == "neutral",]
count_loss_of_function <- counts[counts$effect == "loss_of_function",]
colnames(count_neutral)[3] <- "neutral"
colnames(count_loss_of_function)[3] <- "loss_of_function"
count_neutral <- count_neutral[-c(2)]
count_loss_of_function <- count_loss_of_function[-c(2)]
count_matrix <- merge(count_neutral, count_loss_of_function, by = "gene", all=TRUE)
count_matrix[is.na(count_matrix)] <- 0
View(count_matrix)
```

3. **Implement a statistical test that determines whether a gene has a significantly higher proportion of loss-of-function mutations (excluding mutations with uncertain effects), compared to other genes.**

```r
df_stat_test <- merge(count_matrix, frequency_gene, by = "gene", all = TRUE)
df_stat_test <- df_stat_test[-c(2)]
df_stat_test[is.na(df_stat_test)] <- 0
df_stat_test$proportion <- df_stat_test$loss_of_function/df_stat_test$freq

### z test only for n >= 30 
z_test <- df_stat_test[df_stat_test$freq >= 30,]
z_test <- z_test[order(z_test$proportion, decreasing = TRUE),]
stat_test <- function(z_test, gene_row){  
	p_a <- z_test[gene_row,"proportion"]  
	n_1 <- z_test[gene_row,"freq"]  
	p_b <- mean(z_test[-gene_row, "proportion"])  
	n_2 <- sum(z_test[-gene_row, "freq"])  
	p <- sum(z_test[, "loss_of_function"]) / sum(z_test[, "freq"])  
	q <- 1-p  z_score <- p_a*p_b/(sqrt(p*q/n_1+p*q/n_2))  
	return(c(z_test[gene_row,"gene"], z_score, p, q))
}
```

4. **Identify candidate tumour suppressor genes using this statistical test, adjusting for multiple hypothesis testing. The output table should contain:**
    - **gene symbol**
    - **an estimate of effect size (Z score)**
    - **p value**
    - **q value**

```r
candidates <- data.frame()
for (gene_row in 1:38){  
	candidates <- rbind(candidates, stat_test(z_test, gene_row))
}
colnames(candidates)[1] <- "gene_symbol"
colnames(candidates)[2] <- "sestimate_of_effect_size"
colnames(candidates)[3] <- "p_value"
colnames(candidates)[4] <- "q_value"
```

<img width="787" alt="Screenshot_2023-07-30_at_6 02 25_PM" src="https://github.com/yerimeeei/exercise-r/assets/134043926/2ee71402-5b19-4794-82dc-2bafd561e654">

- z score (estimate of effect size) ≥ 1.96, the difference is significant at 5%.
- Thus, only RB1 gene is significant at 5% resulted by z-test.
5. **Perform a literature search and explain the function of each candidate gene in the context of cancer, as well as specifically in small cell lung cancer.**

### **Retinoblastoma 1 (RB1) gene**

- The RB1 gene provides instructions for making a protein called pRB. This protein acts as a tumor suppressor, which **regulates cell growth and keeps cells from dividing too fast or in an uncontrolled way**.
- **[Epidemiology and *Rb*1 gene of retinoblastoma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3340672/#:~:text=In%20children%20with%20the%20heritable%20genetic%20form%20of%20retinoblastoma%20there,)%2C%20a%20cancer%20may%20develop.)**
    - In children with the heritable genetic form of retinoblastoma there is a mutation on chromosome 13, called the ***Rb*1 gene**. The genetic codes found in chromosomes control the way in which **cells grow and develop within the body. If a portion of the code is missing or altered (mutation), a cancer may develop**.
- [RB Tumor Suppressor in Small Cell Lung Cancer](https://www.notion.so/Programming-exercise-for-R-0a0780faec4949419e402ee8a2383ea6?pvs=21)
    - Small cell lung cancer (SCLC) is an exceptionally aggressive malignancy with limited therapeutic options. **The molecular hallmarks of SCLC include inactivating genomic alterations of the retinoblastoma gene (*RB1*), leading to the loss of Rb protein expression, and concomitant *TP53* alterations.**

### MLL2 gene

- ****[mutation of MLL2 in non–small cell lung carcinoma](https://www.nature.com/articles/srep06036)****
    - Recent studies showed that **MLL2-inactivating mutations** were frequently observed in several cancers, establishing it as a **novel tumor suppressor**. Frequent mutation of MLL2 was also observed in COSMIC lung cancer cohort (64 out of 431, 15%), we thus suggest that **MLL2 might be a novel cancer gene associated with NSCLC.**
    - Somatic mutations of MLL2 were identified only in three tumors in the discovery cohort. **The study found that MLL2 expression was either silenced or significantly reduced in all tumor tissues compared with adjacent non-tumor lung tissues**, **regardless of the mutation status**. All together, **they found that MLL2 was frequently mutated and repressed in NSCLC, supporting its role as a critical tumor suppressor**.

### TP53 gene

- [Role of TP53 mutations in *EGFR*-Mutated Non-Small-Cell Lung Cancer](https://www.notion.so/Programming-exercise-for-R-0a0780faec4949419e402ee8a2383ea6?pvs=21)
    - Tumor protein 53 (*TP53*) gene is the most frequently mutated gene in cancer, including NSCLC. ***TP53* mutations are able to induce carcinogenesis, tumor development and resistance to therapy, influencing patient prognosis and responsiveness to therapy**.
    - The paper assume so far, as observed in other tumors, ***TP53* is a prognostic factor;** however, after many reports concerning its effect on the efficacy of *EGFR*-TKIs, now need to understand some of the molecular processes that link such mutations and resistance to *EGFR*-targeted therapy, as this could guide resistance to therapy.
