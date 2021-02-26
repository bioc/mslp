#' @import data.table
#' @importFrom stats cor na.omit p.adjust sd
#' @importFrom utils combn capture.output
#' @importFrom magrittr %>% set_rownames multiply_by extract extract2
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG %dorng%
NULL

#- deal with . in magrittr
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
