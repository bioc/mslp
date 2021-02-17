# slGene
slGene provides a comprehensive pipeline to predict potential pipeline to predict the potential synthetic lethality partners of 
tumour mutations, based on data of gene expression, mutation profiling and genetic screens in cell lines.

## Installation
``` r
# From GitHub:
devtools::install_github("tidyverse/ggplot2")
```

! More on the vignette.
## Basic usage
Generally the pipleline involves four steps.

**Prepare data**
The mutation, CNA, expression, and zscore data are available in cBioPortal.
``` r
#- Path to input files.
P_mut  <- "data_mutations_extended.txt"
P_cna  <- "data_CNA.txt"
P_expr <- "data_RNA_Seq_v2_expression_median.txt"
P_z    <- "data_RNA_Seq_v2_mRNA_median_Zscores.txt"

res <- pp_tcga(P_mut, P_cna, P_expr, P_z)

saveRDS(res$mut_data, "mut_data.rds")
saveRDS(res$expr_data, "expr_data.rds")
saveRDS(res$zscore_data, "zscore_data.rds")
```

**SLPs from compensationModule**
We identify SLPs that are compensated for loss of funciton of mutations.
``` r
zscore_data <- readRDS("../../01_data/brca/zscore_data.rds")
mutcna_data <- readRDS("../../01_data/brca/mut_data.rds")

res <- compSLP(zscore_data, mutcna_data, ncore = 4, fast = TRUE, huge = TRUE)
saveRDS(res, file = "compSLP_res.rds")
```

**SLPs from correlatoinModule**
We identify SLPs that are correlated with mutations in WT patients
``` r
expr_data   <- readRDS("expr_data.rds")
mutcna_data <- readRDS("mut_data.rds")

#- Filter results by im threshold
res       <- corr_slp(expr_data, mutcna_data, mutgene = mutgene, ncore = 4)
im_thresh <- 0.0016
res       <- res[im >= im_thresh]
saveRDS(res, "corrSLP_res.rds")
```

**Identify consensus SLP**
Among the above predicted SLPs, we identify SLPs are screen hits of same mutations
among pairs of related cancer cell lines.

``` r
#- First we merged the predicted primary SLPs.
comp_res <- readRDS("compSLP_res.rds")
corr_res <- readRDS("corrSLP_res.rds")
comb_res <- merge_slp(comp_res, corr_res)

#- Read mutation and genetic screen data of related cancer cell lines. Note the cells 
# line id should match.
cell_mut    <- readRDS("cell_mut.rds")
cell_screen <- readRDS("cell_screen.rds")

#- First, identify screen hits are they are SLPs for each cell line.
allcell <- intersect(cell_mut$Tumor_Sample_Barcode, cell_screen$cell_line)
scr_res <- lapply(allcell, scr_slp, comb_res, cell_mut, cell_screen)

#- Identifiy consensus SLPs.
k_res <- cons_slp(scr_res, tumour_slp = comb_res, method = "Kappa_test", ncore = 4)
```

See details in the vignette.
