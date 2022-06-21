#' Process tumour genomic data
#'
#' Preprocess mutation, cna, expression and zscore datsets in TCGA PanCancer Atlas by cBioPortal.
#'
#' @param p_mut path of muation data, like "data_mutations_uniprot.txt" provided by cBioPortal.
#' @param p_cna path of copy number variation data, like "data_CNA.txt".
#' @param p_exprs path of normalized RNAseq expression data, like "data_RNA_Seq_v2_expression_median.txt".
#' @param p_score path of zscore data, like "data_RNA_Seq_v2_mRNA_median_Zscores.txt".
#' @param freq_thresh threshold to select recurrent mutations.
#' @param expr_thresh threshold to remove low expression genes.
#' @param hypermut_thresh threshold for hpyermutations.
#' @details
#'   It is designed to process the TCGA data provided by cBioPortal. In mutation data, "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
#'   "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation" are selected for the downstream analysis,
#'   In CNA data, genes with GISTIC value equal to -2 are used. Patients with hypermutations are removed.
#'   Low expression genes, or genes that are not detected in any patient are filtered out.
#' @return
#'   Return a list of mut_data, expr_data and zscore_data, while expr_data and zscore_data are matrix (entrez_id by patients),
#'   mut_data is a data.table with two columns of "patientid" and "mut_entrez".
#' @references
#'   Cerami et al. The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data. Cancer Discovery. May 2012 2; 401.
#'   Gao et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Sci. Signal. 6, pl1 (2013).
#' @examples
#' #- See vignette for more details.
#' if (FALSE) {
#' P_mut  <- "data_mutations_extended.txt"
#' P_cna  <- "data_CNA.txt"
#' P_expr <- "data_RNA_Seq_v2_expression_median.txt"
#' P_z    <- "data_RNA_Seq_v2_mRNA_median_Zscores.txt"
#' res    <- pp_tcga(P_mut, P_cna, P_expr, P_z)
#' saveRDS(res$mut_data, "mut_data.rds")
#' saveRDS(res$expr_data, "expr_data.rds")
#' saveRDS(res$zscore_data, "zscore_data.rds")
#' }
#' @export
pp_tcga <- function(p_mut,
    p_cna,
    p_exprs,
    p_score,
    freq_thresh     = 0.02,
    expr_thresh     = 10,
    hypermut_thresh = 300) {
  ensembl_id <- Gene <- Variant_Classification <- Entrez_Gene_Id <- V1 <- patientid <- mut_entrez <- NULL

  message("(II) Pre-process mutation data.")
  mut_type <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation")
  anno     <- as.data.table(org.Hs.eg.db::org.Hs.egENSEMBL)[match(unique(ensembl_id), ensembl_id)]

  mut_data <- fread(p_mut, select = c("Tumor_Sample_Barcode", "Gene", "Variant_Classification")) %>%
    .[!is.na(Gene)] %>%
    .[Variant_Classification %in% mut_type] %>%
    unique

  #- Add gene symbol
  mut_data <- merge(mut_data, anno, by.x = "Gene", by.y = "ensembl_id") %>%
    .[, -c("Gene", "Variant_Classification")] %>%
    setnames(c("gene_id", "Tumor_Sample_Barcode"), c("mut_entrez", "patientid")) %>%
    unique

  message("(II) Pre-process CNA data.")
  cna_data <- fread(p_cna) %>%
    .[!is.na(Entrez_Gene_Id)] %>%
    unique

  cna_mtx  <- as.matrix(cna_data[, -c(1:2)])
  cna_idx  <- as.data.table(which(cna_mtx == -2, arr.ind = TRUE))

  cna_idx[, `:=` (patientid   = colnames(cna_mtx)[cna_idx$col],
                  mut_entrez  = as.character(cna_data$Entrez_Gene_Id[cna_idx$row]))]

  #- Combine the cna and mutation data to get a broad sense of mutation data, and remove low frequency mutations
  mutcna_data <- unique(rbind(mut_data[, c(1:2)], cna_idx[, c(3:4)]))
  recu_thresh <- mutcna_data$patientid %>% unique %>% length %>% multiply_by(freq_thresh) %>% round
  mutcna_data <- mutcna_data[mutcna_data[, .I[.N >= recu_thresh], by = mut_entrez][, V1]]

  message("(II) Pre-process RNAseq expression data.")
  expr_data <- fread(p_exprs) %>%
    .[!is.na(Entrez_Gene_Id)] %>%
    unique

  expr_data <- expr_data[, -c(1:2)] %>%
    as.matrix %>%
    set_rownames(as.character(expr_data$Entrez_Gene_Id)) %>%
    .[rowSums(!is.finite(.)) != ncol(.), ] %>%
    .[rowMeans(., na.rm = TRUE) >= expr_thresh, ]

  message("(II) Pre-process z-score data.")
  zscore_data <- fread(p_score) %>%
    .[!is.na(Entrez_Gene_Id)] %>%
    unique

  zscore_data <- zscore_data[, -c(1:2)] %>%
    as.matrix %>%
    set_rownames(as.character(zscore_data$Entrez_Gene_Id)) %>%
    .[rowSums(!is.finite(.)) != ncol(.), ] %>%
    .[intersect(rownames(.), rownames(expr_data)), ]

  #- Same genes set in zscore dand expr.
  expr_data <- expr_data[rownames(zscore_data), ]

  #- Remove patietns with hypermutations
  mutcna_data <- Reduce(intersect, list(colnames(expr_data), colnames(zscore_data), unique(mutcna_data$patientid))) %>%
    table(mutcna_data$patientid)[.] %>%
    .[. <= hypermut_thresh] %>%
    names %>%
    {mutcna_data[patientid %in% .]}

  return(list(mut_data = mutcna_data, expr_data = expr_data, zscore_data = zscore_data))
}

#' Merge SLPs
#'
#' Merge predcted SLPs from comp_slp and corr_slp.
#'
#' @param comp_data predicted SLPs from \code{\link{comp_slp}}.
#' @param corr_data predicted SLPs from \code{\link{corr_slp}}.
#' @return A data.table.
#'   \describe{
#'     \item{mut_entrez}{Entrez ids of mutations.}
#'     \item{mut_symbol}{Gene symbols of mutations.}
#'     \item{slp_entrez}{Entrez ids of SLPs.}
#'     \item{slp_symbol}{Gene symbols of SLPs.}
#'     \item{pvalue}{p_value from \code{\link[RankProd]{RankProducts}}.}
#'     \item{fdr}{"BH" adjusted pvalue via \code{\link[stats]{p.adjust}}.}
#'     \item{im}{The importance value returned by \code{\link{genie3}}.}
#'     \item{dualhit}{Whether the slp is identified by \code{\link{corr_slp}} and \code{\link{comp_slp}}.}
#' }
#' @examples
#' data("example_z")
#' data("comp_mut")
#' comp_res <- comp_slp(example_z, comp_mut)
#'
#' data("example_expr")
#' data("corr_mut")
#' corr_res <- corr_slp(example_expr, corr_mut)
#'
#' res <- merge_slp(comp_res, corr_res)
#' @export
merge_slp <- function(comp_data, corr_data) {
  mut_entrez <- idx <- dualhit <- slp_entrez <- pvalue <- im <- NULL

  comp_data <- comp_data[, idx := paste(mut_entrez, slp_entrez, sep = "_")]
  corr_data <- corr_data[, idx := paste(mut_entrez, slp_entrez, sep = "_")]

  comm_data <- merge(comp_data, corr_data[, -c("mut_entrez", "slp_entrez", "mut_symbol", "slp_symbol")], by.x = "idx", by.y = "idx") %>%
    setcolorder(c("idx", "mut_entrez", "mut_symbol", "slp_entrez", "slp_symbol", "pvalue", "fdr", "im"))

  comp_data <- comp_data[!(idx %in% comm_data$idx)] %>%
    .[, `:=` (im = NA)] %>%
    setcolorder(c("idx", "mut_entrez", "mut_symbol", "slp_entrez", "slp_symbol", "pvalue", "fdr", "im"))

  corr_data <- corr_data[!(idx %in% comm_data$idx)] %>%
    .[, `:=` (pvalue = NA, fdr = NA)] %>%
    setcolorder(c("idx", "mut_entrez", "mut_symbol", "slp_entrez", "slp_symbol", "pvalue", "fdr", "im"))

  res <- rbindlist(list(comm_data, comp_data, corr_data)) %>%
    .[, dualhit := "NO"] %>%
    .[!is.na(im) & !is.na(pvalue), dualhit := "YES"] %>%
    .[, idx := NULL]

  return(res)
}

#' Estimate the importance threshold for GENIE3
#'
#' Estimate the importance threshold based on repetition GENIE3 results via ROC.
#'
#' @param permu_data permuated \code{\link{corr_slp}} results.
#' @param fdr_thresh fdr threshold to selected "TRUE" SLPs.
#' @details
#'   We first generate a SLPs by repetition matrix from repetition GENIE3 results.
#'   SLPs with high im value in repetitions are selected and condsidered as "TRUE" SLPs via the rank product algorithm.
#'   Then for each repetion, we perform receiver operating characteristic curve analysis and select an optimal threshold by "youden" approach.
#'   The optimal thresholds are averaged to get the final threshold.
#' @return A data.table with mut_entrez (mutation entrez_id) and roc_thresh (estimated im threshold).
#' @examples
#' #- Toy examples.
#' require(future)
#' require(doFuture)
#' plan(multisession, workers = 2)
#' data(example_expr)
#' data(corr_mut)
#' mutgene    <- sample(intersect(corr_mut$mut_entrez, rownames(example_expr)), 2)
#' nperm      <- 5
#' res        <- lapply(seq_len(nperm), function(x) corr_slp(example_expr, corr_mut, mutgene = mutgene))
#' roc_thresh <- est_im(res)
#' plan(sequential)
#' @export
est_im <- function(permu_data, fdr_thresh = 0.001) {
  permu_id <- gene <- im <- slp_entrez <- fdr <- mut_entrez <- NULL

  mutgene <- lapply(permu_data, function(x) x$mut_entrez) %>% unlist %>% unique

  suppressPackageStartupMessages(
    im_res <- foreach(gene = mutgene) %dopar% {
      #- Get the slp_entrez for a mutation.
      mtx <- lapply(seq_along(permu_data), function(x) permu_data[[x]][mut_entrez == gene][, permu_id := paste0("perm_", x)][, mut_entrez := NULL])

      mtx[lengths(mtx) == 0] <- NULL

      if (length(mtx) > 0) {
        mtx <- rbindlist(mtx)

        mtx <- dcast(mtx, permu_id ~ slp_entrez, value.var = "im", fill = NA) %>%
          .[, permu_id := NULL] %>%
          as.matrix %>%
          t

        #- Top im.
        invisible(capture.output(rp_res <- RankProd::RankProducts(mtx,
                                         rep(1, ncol(mtx)),
                                         logged           = FALSE,
                                         gene.names       = rownames(mtx),
                                         calculateProduct = FALSE)))

        #- SLPs with high importances are selected.
        true_slp_res <- data.table(slp_entrez = rownames(rp_res$pval), fdr = p.adjust(rp_res$pval[, 2], method = "BH")) %>%
          .[fdr <= fdr_thresh] %>%
          .[order(fdr)]

        #- ROC analysis, for each run of permutation data.
        roc_thresh  <- apply(mtx, 2, fn_ROC, true_slp_res$slp_entrez) %>%
          .[. != -1] %>%
          mean %>%
          {data.table(mut_entrez = gene, roc_thresh = .)}

        return(roc_thresh)
      }
    }
  )

  im_res[lengths(im_res) == 0] <- NULL
  if (length(im_res) > 0) {
    im_res <- rbindlist(im_res)
  } else {
    message("(II) im estimation is NULL for all mutations.")
  }

  return(im_res)
}

fn_ROC <- function(x, true_slp) {
  slp_entrez <- tr_slp <- NULL

  true_slp <- intersect(true_slp, na.omit(x) %>% names)
  roc_data <- data.table(slp_entrez = names(x), im = x, tr_slp = rep("NO", length(x))) %>%
    .[slp_entrez %in% true_slp, tr_slp := "YES"] %>%
    .[, tr_slp := factor(tr_slp)] %>%
    as.data.frame %>%
    set_rownames(names(x))

  roc_res <- try(pROC::roc(tr_slp ~ im, roc_data, levels = c("NO", "YES"), direction = "<"), silent = TRUE)

  if (inherits(roc_res, "try-error")) {
    roc_thresh <- -1
  } else {
    roc_thresh <- pROC::coords(roc_res, "best", best.method = "youden", transpose = TRUE)[1]
  }

  return(roc_thresh)
}
