#' AFatlasqtl: Adaptive-focus implementation of atlasQTL
#'
#' @keywords internal
"_PACKAGE"
#'
#' @section atlasqtl functions: atlasqtl, assign_bFDR, map_hyperprior_elicitation, print.atlasqtl, set_hyper, set_init, summary.atlasqtl.
#'
#' @name AFatlasqtl
#' @useDynLib AFatlasqtl, .registration = TRUE
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor dnorm median pnorm qnorm rbeta rbinom rgamma rnorm setNames uniroot var
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline legend matplot points
#' @import tictoc
NULL