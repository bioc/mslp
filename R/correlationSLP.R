#' Identify SLPs via correlation
#'
#' Identify SLPs of mutations based on co-expression. GENIE3 is employed to
#' find genes highly correlated with mutations in wide type patients.
#'
#' @param expr_data an expression matrix, genes by patients.
#' @param mut_data a data.table with columns "patientid" and "mut_entrez".
#' @param ncore number of cores for parallel computing.
#' @param mutgene identify SLPs for sepecific muatation (gene symbols). If NULL (by default), the intersection genes between expr_data and mut_data are used.
#' @param im_thresh minimum importance threshold.
#' @param topgene top N genes above the \code{im_thresh}.
#' @param ... further parameters to \code{\link{genie3}}.
#' @return A data.table with predicted SLPs.
#'   \describe{
#'     \item{mut_entrez}{Entrez ids of mutations.}
#'     \item{mut_symbol}{Gene symbols of mutations.}
#'     \item{slp_entrez}{Entrez ids of SLPs.}
#'     \item{slp_symbol}{Gene symbols of SLPs.}
#'     \item{fdr}{"BH" adjusted pvalue via \code{\link[stats]{p.adjust}}.}
#'     \item{im}{The importance value returned by \code{\link{genie3}}.}
#' }
#' @export
corr_slp <- function(expr_data, mut_data, ncore = 2, mutgene = NULL, im_thresh = 0.001, topgene = 2000, ...) {
  i <- symbol <- im <- mut_entrez <- NULL

  if (is.null(mutgene)) {
    #- Mutations detected in expression data.
    mutgene <- intersect(unique(mut_data$mut_entrez), rownames(expr_data))
  } else {
    #- Mutations given.
    mutgene <- intersect(mutgene, rownames(expr_data))
  }

  message("(II) Number of mutations: ", length(mutgene), ".")

  doFuture::registerDoFuture()
  future::plan(future::multisession, workers = ncore)

  suppressPackageStartupMessages(
    genie3_res <- foreach(i = mutgene) %dorng% {
      fn_sub_corr_slp(i, expr_data, mut_data, im_thresh, topgene, ...)
    }
  )

  genie3_res[lengths(genie3_res) == 0] <- NULL

  if(length(genie3_res) > 0) {
    genie3_res <- rbindlist(genie3_res)
    anno       <- as.data.table(org.Hs.eg.db::org.Hs.egSYMBOL)[match(unique(symbol), symbol)]
    genie3_res <- merge(genie3_res, anno, by.x = "mut_entrez", by.y = "gene_id", all.x = TRUE) %>%
      setnames("symbol", "mut_symbol") %>%
      merge(anno, by.x = "slp_entrez", by.y = "gene_id", all.x = TRUE) %>%
      setnames("symbol", "slp_symbol") %>%
      setcolorder(c("mut_entrez", "mut_symbol", "slp_entrez", "slp_symbol", "im")) %>%
      .[order(-im)]

    return(genie3_res)
  } else {
    message("(II) No potential SLPs from correlationModule.")
  }
}

#- Internal function running GENIE3.
fn_sub_corr_slp <- function(gene, expr_data, mut_data, im_thresh, topgene, ...) {
  mut_entrez <- im <- slp_entrez <- patientid <- NULL

  mut_patient <- mut_data[mut_entrez == gene, unique(patientid)]

  if (length(mut_patient) > 0) {
    comut    <- setdiff(mut_data[patientid %in% mut_patient, unique(mut_entrez)], gene)
    sub_expr <- expr_data[!(rownames(expr_data) %in% comut), !(colnames(expr_data) %in% mut_patient)]
  } else {
    sub_expr <- expr_data
  }

  rm(expr_data)
  gc()

  genecorr <- cor(sub_expr[gene, ], t(sub_expr[rownames(sub_expr) != gene, ])) %>% extract(1, ) %>% na.omit %>% .[. > 0]

  if (length(genecorr) > 0) {
    pred_gene  <- names(genecorr)
    genie3_mtx <- sub_expr[c(gene, pred_gene), ]
    genie3_res <- genie3(genie3_mtx, ngene = 1, ...)

    rm(genie3_mtx, sub_expr)
    gc()

    im_res <- getlink(genie3_res[, 1, drop = FALSE]) %>%
      setnames(c("from.gene", "to.gene"), c("slp_entrez", "mut_entrez")) %>%
      setcolorder("mut_entrez") %>%
      .[im >= im_thresh] %>%
      .[1:topgene] %>%
      na.omit
    return(im_res)
  }
}
