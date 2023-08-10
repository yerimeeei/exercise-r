# Transcriptomic data normalization

[Transcriptomic Data Normalization](https://www.notion.so/Transcriptomic-Data-Normalization-13b30a5acfa048faa533aa0837925af1?pvs=21)

The transcriptomic data is in `data/expr_sclc_ucologne_2015.rds`

1. **Perform an appropriate log transformation on the data.**

```r
# Transcriptomic data normalization
transcriptomic_data <- readRDS("/Users/rimi/bioinformatics/DV/exercise-r/data/expr_sclc_ucologne_2015.rds")

## 1. Perform an appropriate log transformation on the data.
library(edgeR)transcriptomic_dgList <- DGEList(counts = transcriptomic_data)
transcriptomic_dgList <- calcNormFactors(transcriptomic_dgList, method="TMM")
transcriptomic_cpm <- cpm(transcriptomic_dgList, log=TRUE, prior.count=1, normalized.lib.sizes = TRUE)
```

2. **Implement a median polish algorithm from scartch.**
- [Writing R package from scratch](https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/)
- [The median polish](http://mgimond.github.io/ES218/two_way.html)
- [Youtube: median polish](https://www.youtube.com/watch?v=RtC9ZMOYgk8)

```r
install.packages("devtools")
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)setwd('/Users/rimi/')
create("medianpolish")

# Function 'update_col'
update_col <- function(A, n_row, n_col){  
	# Find median by col  
	for (i in seq(1, n_col)){
    col <- A[,i]
    col <- col[1:n_row-1]
    # Replace the last row
    med <- median(col)
    A[n_row,i] <- A[n_row,i]+med
    # Update the (nxm) cells
    for(j in seq(1,n_row-1)){
      A[j,i] <- A[j,i]-med
    }
  }
  return(A)}
### Function 'update_row'
update_row <- function(A, n_row, n_col){
  for (i in seq(1,n_row)){
    # Find median by row
    row <- A[i,]
    row <- row[1:n_col-1]
    # Replace the last column
    med <- median(row)
    A[i,n_col] <- A[i,n_col]+med
    # Update the (nxm) cells
    for(j in seq(1,n_col-1)){
      A[i,j] <- A[i,j]-med
    }
  }
  return(A)}

medianpolish <- function(A){
  # padding - (nxm) matrix to (n+1)x(m+1) matrix
  A <- rbind(A, numeric(ncol(A))) 
  A <- cbind(A, numeric(nrow(A))) 
  n_row <- nrow(A)
  n_col <- ncol(A)
  A_prev <- A  

  A <- update_row(A, n_row, n_col)
  A <- update_col(A, n_row, n_col)  

  while( !identical(A, A_prev) ){
    A_prev <- A
    A <- update_row(A, n_row, n_col)
    A <- update_col(A, n_row, n_col)
  }  
  row_effect <- A[-n_row, n_col]
  col_effect <- A[n_row, -n_col]
  m <- A[n_row,n_col]
  output <- c( m, row_effect, col_effect)
  output <- list(
    'Overall' = m,
    'Row Effect' = row_effect,
    'Col Effect' = col_effect,
    'Residuals' = A[-n_row,-n_col]
  )
  return(output)
}
```

3. **Compare the residuals of your algorithm and `stats::medpolish`.**

```r
A <- transcriptomic_cpm
my_output <- medianpolish(A)
r_output <- medpolish(A)
```

- redisuals of my_output

![my_output](https://github.com/yerimeeei/exercise-r/assets/134043926/7beaa7c9-626c-4ff9-b01e-120513c16bcb)

- residuals of r_output

![r_output](https://github.com/yerimeeei/exercise-r/assets/134043926/6d8d8dbc-44fb-4a7f-ad9b-55f7738151a6)

4. **Plot heatmaps of the results before and after median polish.**

```r
library(tidyverse)
library(magrittr)
library(pheatmap)
boxplot(transcriptomic_cpm, las = 2)
heatmap_before <- pheatmap(transcriptomic_cpm, scale = "row",
          main = "Transcriptomic data Heatmap before medianpolish",
         xlab = "Samples",
         ylab = "Genes",
         fontsize_col = 7,
 fontsize_row = 7)
heatmap_before
```

- heatmap before median polish

![pheatmap_transcriptomic_cpm_(scale)](https://github.com/yerimeeei/exercise-r/assets/134043926/9ab4864f-9c87-4d95-b91f-68f5e848add0)

```r
heatmap_after <- pheatmap(transcriptomic_medpolish, scale = "row",
                          main = "Transcriptomic data Heatmap after medianpolish",
                          xlab = "Samples",
                          ylab = "Genes",
                          fontsize_col = 7, fontsize_row = 7)
heatmap_after
```

- heatmap after median polish

![pheatmap_transcriptomic_cpm_after_medianpolish](https://github.com/yerimeeei/exercise-r/assets/134043926/5b27c7c3-605d-4152-8c00-29c33438cd27)

5. **Output the median polished residual matrix as the normalized transcriptomic data.**

```r
# update medianpolish function
medianpolish <- function(A){
  # padding - (nxm) matrix to (n+1)x(m+1) matrix
  A <- rbind(A, numeric(ncol(A))) 
  A <- cbind(A, numeric(nrow(A))) 
  n_row <- nrow(A)
  n_col <- ncol(A)
  A_prev <- A
  A <- update_row(A, n_row, n_col)
  A <- update_col(A, n_row, n_col)  

  while( !identical(A, A_prev) ){
    A_prev <- A
    A <- update_row(A, n_row, n_col)
    A <- update_col(A, n_row, n_col)
  }  
  row_effect <- A[-n_row, n_col]
  col_effect <- A[n_row, -n_col]
  m <- A[n_row,n_col]
  output <- c( m, row_effect, col_effect)
  output <- list(
    'Overall' = m,
    'Row Effect' = row_effect,
    'Col Effect' = col_effect,
    'Residuals' = A[-n_row,-n_col])
  residual_matrix <- data.frame(A[-n_row, -n_col]
  )  
return(residual_matrix)
}
## median polished residual matrix (normalized transcriptomic data)
transcriptomic_medpolish <- medianpolish(A)
```

![Screenshot_2023-08-09_at_6 48 43_PM](https://github.com/yerimeeei/exercise-r/assets/134043926/4dc4d9d2-d639-44e6-8aa1-6495fbeb44ac)
