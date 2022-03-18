#' Project Drive screen data for two brca cell lines
#'
#' Screen hits are selected with RSA <= -2.
#'
#' @format A data.table with three columns: "screen_entrez", "screen_symbol" and "cell_line".
"brca_screen"

#' Mutations in two brca cell lines.
#'
#' Mutation gene ids and related cell lines.
#'
#' @format A data.table.
"brca_mut"

#' Patients mutations to be use in the comp_slp
#'
#' Mutations and related TCGA ids.
#'
#' @format A data.table.
"comp_mut"

#' Patients mutations to be use in the corr_slp
#'
#' Mutations and related TCGA ids.
#'
#' @format A data.table.
"corr_mut"

#' Expression data to be used in comp_slp
#'
#' Expresion matrix, genes by samples.
#'
#' @format A matrix.
"example_expr"

#' Expression data to be used in corr_slp
#'
#' Z score matrix, genes by samples.
#'
#' @format A matrix.
"example_z"
