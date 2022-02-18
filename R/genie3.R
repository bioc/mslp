#' Run GENIE3
#'
#' Calculate the weight matrix between genes via randomForest, modified from original codes by Huynh-Thu, V.A.
#'
#' @param expr.matrix  exrepssion matrix (genes by samples).
#' @param ngene an integer, only up to ngene (included) targets (responsible variables).
#' @param K choice of number of input genes randomly, must be one of "sqrt", "all", an integar.
#' @param nb.trees number of trees in ensemble for each target gene (default 1000).
#' @param input.idx subset of genes used as input genes (default all genes). A vector of indices or gene names is accepted.
#' @param importance.measure type of variable importance measure, "IncNodePurity" or "\%IncMSE".
#' @param seed random number generator seed for replication of analyses.
#' @param trace index of currently computed gene is reported (default FALSE).
#' @param ... parameter to randomForest.
#' @return A weighted adjacency matrix of inferred network, element w_ij (row i, column j) gives the importance of the link from regulatory gene i to target gene i.
#' @references
#'   Huynh-Thu, V.A., Irrthum, A., Wehenkel, L., and Geurts, P. (2010). Inferring Regulatory Networks from Expression Data Using Tree-Based Methods. PLoS ONE 5, e12776.
#' @examples
#' \dontrun{
#' mtx <- matrix(sample(1000, 100), nrow = 5)
#' mtx <- rbind(mtx[1, ] * 2 + rnorm(20), mtx)
#' res <- genie3(mtx, K = 1, nb.trees = 100)
#' }
#' @export
genie3 <- function(expr.matrix,
    ngene              = NULL,
    K                  = "sqrt",
    nb.trees           = 1000,
    input.idx          = NULL,
    importance.measure = "IncNodePurity",
    # seed               = NULL,
    trace              = FALSE, ...) {
  # set random number generator seed if seed is given
  # if (!is.null(seed)) set.seed(seed)

  # to be nice, report when parameter importance.measure is not correctly spelled
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
  }

  # transpose expression matrix to (samples x genes)
  expr.matrix   <- t(expr.matrix)
  # normalize expression matrix
  expr.matrix   <- apply(expr.matrix, 2, function(x) {(x - mean(x)) / sd(x)})
  # setup weight matrix
  num.samples   <- dim(expr.matrix)[1]
  num.genes     <- dim(expr.matrix)[2]
  gene.names    <- colnames(expr.matrix)
  weight.matrix <- matrix(0.0, nrow = num.genes, ncol = num.genes)
  rownames(weight.matrix) <- gene.names
  colnames(weight.matrix) <- gene.names

  # get number of input genes, names of input genes
  if (is.null(input.idx)) {
    num.input.genes  <- num.genes
    input.gene.names <- gene.names
  } else {
    num.input.genes <- length(input.idx)
    # input gene indices given as integers
    if (is.numeric(input.idx)) {
      input.gene.names <- gene.names[input.idx]
    # input gene indices given as names
    } else {
      input.gene.names  <- input.idx
      # for security, abort if some input gene name is not in gene names
      missing.gene.names <- setdiff(input.gene.names, gene.names)
      if (length(missing.gene.names) != 0) {
        for (missing.gene.name in missing.gene.names) {
          message(paste("Gene ", missing.gene.name, " was not in the expression matrix\n", sep = ""))
        }
        stop("Aborting computation")
      }
    }
  }

  if (is(k, "numeric")) {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.input.genes))
  } else if (K == "all") {
    mtry <- num.input.genes - 1
  } else {
    stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
  }

  if (trace) {
    message(paste("Starting RF computations with ", nb.trees,
                  " trees/target gene,\nand ", mtry,
                  " candidate input genes/tree node\n",
                  sep = ""))
  }

  if (is.null(ngene)) ngene <- num.genes
  for (target.gene.idx in seq_len(ngene)) {
    if (trace) message(paste("Computing gene ", target.gene.idx, "/", ngene, "\n", sep = ""))

    target.gene.name       <- gene.names[target.gene.idx]
    these.input.gene.names <- setdiff(input.gene.names, target.gene.name)

    x  <- expr.matrix[, these.input.gene.names]
    y  <- expr.matrix[, target.gene.name]
    rf <- randomForest::randomForest(x, y, mtry = mtry, ntree = nb.trees, importance = TRUE, ...)
    im <- randomForest::importance(rf)[, importance.measure]
    im.names <- names(im)
    #- The columns are the target genes, i.e., the response variable.
    weight.matrix[im.names, target.gene.name] <- im
  }

  weight.matrix <- weight.matrix / num.samples
  return(weight.matrix)
}

#' Get sorted list of regulatory links in GENIE3 results
#'
#' Take genie3 output and sort the links.
#'
#' @param weight.matrix a weighted adjacency matrix as returned by get.weight.matrix.targ.
#' @param report.max maximum number of links to report (default all links).
#' @return A data.table of links with columns "from.gene", "to.gene", "im".
#' @export

getlink <- function(weight.matrix, report.max = NULL) {
  im <- NULL

  res <- data.table(from.gene = rownames(weight.matrix)[c(row(weight.matrix))],
                    to.gene   = colnames(weight.matrix)[c(col(weight.matrix))],
                    im        = c(weight.matrix)) %>%
    .[im > 0] %>%
    .[order(-im)]

  if (!is.null(report.max) && nrow(res) > report.max) res <- res[seq_len(report.max)]

  return(res)
}
