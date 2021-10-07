#' Identify SLPs via compensation
#'
#' Identify SLPs compensating for the loss of function of mutations. The up-regulated SLPs are selected via the rank prodcuts algorithm, with option calculateProduct = FALSE for a robust results and capacity on large datasets.
#'
#' @param zscore_data a matrix (genes by patients) reflecting gene expression related to wide type samples. For example, generated from \code{\link{pp_tcga}}.
#' @param mut_data a data.table with columns "patientid" and "mut_entrez".
#' @param ncore number of cores for parallel computing.
#' @param mutgene identify SLPs for sepecific muatation (gene symbols). If NULL (by default), the intersection genes between zscore_data and mut_data are used.
#' @param positive_perc keep genes with postive zscore in at least positive_perc * number of mutation patients.
#' @param p_thresh pvalue threshold to filter out results.
#' @param ... additional parameters to \code{\link[RankProd]{RankProducts}}.
#' @return A data.table with predicted SLPs.
#'   \describe{
#'     \item{mut_entrez}{Entrez ids of mutations.}
#'     \item{mut_symbol}{Gene symbols of mutations.}
#'     \item{slp_entrez}{Entrez ids of SLPs.}
#'     \item{slp_symbol}{Gene symbols of SLPs.}
#'     \item{pvalue}{p_value from \code{\link[RankProd]{RankProducts}}.}
#'     \item{fdr}{"BH" adjusted pvalue via \code{\link[stats]{p.adjust}}.}
#' }
#' @export
comp_slp <- function(zscore_data,
    mut_data,
    ncore         = 2,
    mutgene       = NULL,
    positive_perc = 0.5,
    p_thresh      = 0.01,
    ...) {
  i <- symbol <- fdr <- mut_entrez <- V1 <- NULL

  #- Mutations found in at least two patients.
  mut_lite <- mut_data[mut_data[, .I[.N >= 2], by = mut_entrez][, V1]]

  if (is.null(mutgene)) {
    mutgene <- intersect(mut_lite$mut_entrez, rownames(zscore_data))
  } else {
    mutgene <- intersect(mutgene, rownames(zscore_data))
  }

  message("(==) Number of mutations: ", length(mutgene), ".")

  if(ncore > 1) {
    doFuture::registerDoFuture()
    future::plan(future::multisession, workers = ncore)

    suppressPackageStartupMessages(
      res <- foreach (i = mutgene) %dopar% {
        fn_sub_comp_slp(i, zscore_data, mut_lite, positive_perc = positive_perc, p_thresh = p_thresh, ...)
      }
    )
  } else {
    suppressPackageStartupMessages(
      res <- foreach (i = mutgene) %do% {
        fn_sub_comp_slp(i, zscore_data, mut_lite, positive_perc = positive_perc, p_thresh = p_thresh, ...)
      }
    )
  }

  res[lengths(res) == 0] <- NULL

  if (length(res) > 0) {
    res <- rbindlist(res)

    anno <- as.data.table(org.Hs.eg.db::org.Hs.egSYMBOL)[match(unique(symbol), symbol)]
    res  <- merge(res, anno, by.x = "mut_entrez", by.y = "gene_id", all.x = TRUE) %>%
      setnames("symbol", "mut_symbol") %>%
      merge(anno, by.x = "slp_entrez", by.y = "gene_id", all.x = TRUE) %>%
      setnames("symbol", "slp_symbol") %>%
      setcolorder(c("mut_entrez", "mut_symbol", "slp_entrez", "slp_symbol", "pvalue", "fdr")) %>%
      .[order(fdr)]

    return(res)
  } else {
    message("(II) No potential SLPs from compensationModule.")
  }
}

#- Internal function running rankprod.
fn_sub_comp_slp <- function(gene,
    zscore_data,
    mut_data,
    positive_perc,
    p_thresh,
    ...) {
  patientid <- mut_entrez <- pvalue <- fdr <- slp_entrez <- NULL

  mut_patient <- mut_data[mut_entrez == gene, unique(patientid)]

  if (length(mut_patient) >= 2) {
    comut <- mut_data[patientid %in% mut_patient, unique(mut_entrez)]

    #- Remove co-mutated genes, and keep a clean data.
    mtx <- zscore_data[!(rownames(zscore_data) %in% comut), colnames(zscore_data) %in% mut_patient] %>%
      na.omit %>%
      .[rowSums(. > 0) >= (positive_perc * ncol(.)), ]

    if (nrow(mtx) > 0) {
      wd_mtx <- zscore_data[, !(colnames(zscore_data) %in% mut_patient)] %>% extract(rownames(mtx), )
      mtx    <- mtx[rowMeans(mtx) > rowMeans(wd_mtx, na.rm = TRUE), ]

      invisible(capture.output(res <- RankProd::RankProducts(mtx, cl = rep(1, ncol(mtx)), logged = TRUE, gene.names = rownames(mtx), calculateProduct = FALSE, ...)))
      compres <- data.table(mut_entrez = gene, slp_entrez = rownames(res$pval), pvalue = res$pval[, 2], fdr = p.adjust(res$pval[, 2], method = "BH")) %>%
                            .[pvalue <= p_thresh] %>%
                            .[order(fdr)]

      return(compres)
    } else {
      message("(II) Empty matrix.")
    }

  } else {
    message("(II) Less than two patients have mutations in: ", gene, ".")
  }
}
