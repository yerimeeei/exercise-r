# Programming Exercises for R

## Mutation data analysis

The mutation data is in `data/mutations_sclc_ucologne_2015.rds`.

1. Identify the top 10 most frequently mutated genes.
   Identify samples whose mutation count is in the 80 to 90 percentile.
2. Categorize variants based on their expected effects using the
   `data/mutation_effects.tsv` table.
   Generate a count matrix containing the numbers of loss-of-function
   and neutral mutations for each gene.
3. Implement a statistical test that determines whether a gene has a
   significantly higher proportion of loss-of-function mutations
   (excluding mutations with uncertain effects),
   compared to other genes.
4. Identify candidate tumour suppressor genes using this statistical test,
   adjusting for multiple hypothesis testing.
   The output table should contain:
     - gene symbol
     - an estimate of effect size
     - p value
     - q value
5. Perform a literature search and explain the function of each candidate
   gene in the context of cancer, as well as specifically in small cell 
   lung cancer.

## Transcriptomic data normalization

The transcriptomic data is in `data/expr_sclc_ucologne_2015.rds`.

1. Perform an appropriate log transformation on the data.
2. Implement a median polish algorithm from scartch.
3. Compare the residuals of your algorithm and `stats::medpolish`.
4. Plot heatmaps of the results before and after median polish.
5. Output the median polished residual matrix as the normalized transcriptomic data.

## Differential gene expression analysis

The clinical data is in `data/pheno_sclc_ucologne_2015.tsv`.

1. Define two groups of tumours as early stage (stages I-II) vs. advanced stage
   tumours (stages III-IV), while excluding samples missing stage information.
2. Identify genes that differentially expressed in early vs. advanced stage
   tumours using an appropriate R package.

## Integrative analysis

The structural variant data is in `data/sv_sclc_ucologne_2015.tsv`.

1. For each gene involved in a structural variant (SV), determine the expression 
   level of the gene in the sample that harbours the SV.
2. Identify SVs that satisfy the following critiera:
      - The involved pair of genes both have elevated expression levels
        in samples with the SV compared to samples without the SV.
      - The second gene in the pair is in frame.

