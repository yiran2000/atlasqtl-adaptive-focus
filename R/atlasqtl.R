# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#

#' Flexible sparse regression for hierarchically-related responses using 
#' annealed variational inference. 
#'
#' @param Y Response data matrix of dimension n x q, where n is the number of
#'   samples and q is the number of response variables.
#' @param X Input matrix of dimension n x p, where p is the number of candidate
#'   predictors. \code{X} cannot contain NAs. No intercept must be supplied.
#' @param p0 Vector of size 2 whose entries are the prior expectation and 
#'   variance of the number of predictors associated with each response.
#'   Must be \code{NULL} if \code{list_init} and \code{list_hyper} are both 
#'   non-\code{NULL}.
#' @param anneal Parameters for annealing scheme. Must be a vector whose first
#'   entry is the type of schedule: 1 = geometric spacing (default), 
#'   2 = harmonic spacing or 3 = linear spacing, the second entry is the initial 
#'   temperature (default is 2), and the third entry is the temperature grid 
#'   size (default is 10). If \code{NULL}, no annealing is performed.
#' @param tol Tolerance for the stopping criterion (default is 0.1).
#' @param maxit Maximum number of iterations allowed (default is 1000).
#' @param user_seed Seed set for reproducible default choices of hyperparameters
#'   (if \code{list_hyper} is \code{NULL}) and initial variational parameters
#'   (if \code{list_init} is \code{NULL}). Default is \code{NULL}, no
#'   seed set.
#' @param verbose Integer specifying the level of verbosity during execution: 0, 
#'   for no message, 1 for standard verbosity (default), 2 for detailed 
#'   verbosity.
#' @param list_hyper An object of class "\code{hyper}" containing the model
#'   hyperparameters. Must be specified using the \code{\link{set_hyper}}
#'   function or set to \code{NULL} for default hyperparameters.
#' @param list_init An object of class "\code{init}" containing the initial
#'   variational parameters. Must be specified using the \code{\link{set_init}} 
#'   function or set to \code{NULL} for a default initialization.
#' @param save_hyper If \code{TRUE}, the hyperparameters used for the model are
#'   returned.
#' @param save_init If \code{TRUE}, the initial variational parameters used for
#'   the inference are returned (note that the size of the resulting objects is
#'   likely to be large). Default is \code{FALSE}.
#' @param full_output If \code{TRUE}, the inferred variational parameters for 
#'   all parameters are returned.
#' @param thinned_elbo_eval If \code{TRUE}, the convergence criterion 
#'   (variational lower bound on the marginal log-likelihood, or ELBO) is not 
#'   evaluated at each iteration, in order to save computational resources, but 
#'   after a given batch of iterations where the batch size is set in an 
#'   adaptive fashion accounting for the number of iterations roughly expected 
#'   until convergence. 
#' @param checkpoint_path Path where to save temporary checkpoint outputs. 
#'   Default is \code{NULL}, for no checkpointing.
#' @param trace_path Path where to save trace plot for the variance of hotspot
#'   propensities. Default is \code{NULL}, for no trace saved.
#' @param add_collinear_back Boolean for re-adding the collinear variables X 
#'   (removed as part of \code{atlasqtl} to ensure the stability of inference) 
#'   to the final posterior summaries. Default is \code{FALSE}.
#' @param CAVI String for CAVI scheme, one of "full", "random", "adaptive". 
#' default is "full"
#' @param start_partial for number of iterations before partial update in random 
#' adaptive CAVI. Default is \code{NULL} for CAVI = "full" and 50 for CAVI = "random" 
#' or "adaptive"
#' @param epsilon_scheme for the choice of epsilon (level of randomness) in the 
#' selection probability. One of "minimum", "iteration" and "ELBO". Only for 
#' adaptive CAVI. If not specified, default is "iteration" with geom_alpha = 0.95 
#' and min_epsilon = 0.
#' @param min_epsilon for the minimum epsilon value. Must be specified between 0-1
#' for CAVI = "adaptive", epsilon_scheme = "minimum". For other cases, default as 0
#' @param geom_alpha for the geometric decreasing speed for CAVI = "adaptive", 
#' epsilon_scheme = "iteration". Default 0.95
#' @param fix_size for the fixed porportion of local factors to update if specified. 
#' Only for the CAVI = "adaptive". Default is \code{NULL}, i.e., choice is purely
#' based on Bernoulli probablity defined by the epsilon_scheme, not fixed size. 
#' @details \code{atlasqtl} implements a flexible hierarchical regression 
#'   framework that allows information-sharing across responses and predictors, 
#'   thereby enhancing the detection of weak effects. \code{atlasqtl} is 
#'   tailored to the detection of hotspot predictors, i.e., predictors 
#'   associated with many responses: it is based on a fully Bayesian model 
#'   whereby the hotspot propensity is modelled using the horseshoe shrinkage 
#'   prior: its global-local formulation shrinks noise globally and hence can 
#'   accommodate highly sparse settings, while being robust to individual 
#'   signals, thus leaving the effects of hotspots unshrunk. Inference is 
#'   carried out using a scalable variational algorithm coupled with a novel 
#'   simulated annealing procedure, which is applicable to large dimensions, 
#'   i.e., thousands of response variables, and allows efficient exploration of 
#'   multimodal distributions.
#' 
#'   The columns of the response matrix \code{Y} are centered within the 
#'   \code{atlasqtl} call, and the columns of the candidate predictor matrix 
#'   \code{X} are standardized.
#'
#' @return An object of class "\code{atlasqtl}" containing the following 
#'   variational estimates and settings:
#'  \item{beta_vb}{Estimated effect size matrix of dimension p x q. Entry (s, t) 
#'                 corresponds to the variational posterior mean 
#'                 (mu_beta_vb_st x gam_vb_st) of the regression effect between 
#'                 candidate predictor s and response t.}
#'  \item{gam_vb}{Posterior inclusion probability matrix of dimension p x q.
#'                Entry (s, t) corresponds to the variational posterior 
#'                probability of association between candidate predictor s 
#'                and response t.}
#'  \item{theta_vb}{Vector of length p containing the variational posterior 
#'                  mean of the hotspot propensity. Entry s corresponds to the 
#'                  propensity of candidate predictor s to be associated with 
#'                  many responses.}
#'  \item{X_beta_vb}{Matrix of dimension n x q equal to the matrix product 
#'                   X beta_vb, where X is the rescaled predictor matrix with 
#'                   possible collinear predictor(s) excluded.}
#'  \item{zeta_vb}{Vector of length q containing the variational posterior 
#'                mean of the response importance. Entry t corresponds to the 
#'                propensity of response t have associated predictors.}
#'  \item{converged}{A boolean indicating whether the algorithm has converged
#'                   before reaching \code{maxit} iterations.}
#'  \item{it}{Final number of iterations.}
#'  \item{lb_opt}{Optimized variational lower bound (ELBO) on the marginal
#'                log-likelihood.}
#'  \item{diff_lb}{Difference in ELBO between the last and penultimate
#'                 iterations (to be used as a convergence diagnostic when 
#'                 \code{maxit} has been reached).}
#'  \item{rmvd_cst_x}{Vectors containing the indices of constant variables in 
#'                    \code{X} removed prior to the analysis.}
#'  \item{rmvd_coll_x}{Vectors containing the indices of variables in \code{X} 
#'                     removed prior to the analysis because collinear to other
#'                     variables. The entry name indicates the corresponding 
#'                     variable kept in the analysis.}
#'  \item{list_hyper, list_init}{If \code{save_hyper}, resp. \code{save_init},
#'                               are \code{TRUE}, the hyperparameters, resp. 
#'                               initial variational parameters, used for 
#'                               inference are saved as output.}
#'  \item{time_init} Time spent on initialization
#'  \item{perform_df} {If \code{eval_perform} is \code{TRUE} Details of timing, 
#'  epsilon, subsample size, ELBO, ELBO_diff 
#'  on each iteration if is 
#'  \item{...}{If \code{full_output} is \code{TRUE} all inferred variational 
#'             parameters are returned.}
#'  
#' @examples
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' # Example with small problem sizes:
#' #
#' n <- 200; p <- 50; p_act <- 10; q <- 100; q_act <- 50
#'
#' # Candidate predictors (subject to selection)
#' #
#' # Here example with common genetic variants under Hardy-Weinberg equilibrium
#' #
#' X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.25), nrow = n)
#'
#' # shuffle indices 
#' shuff_x_ind <- sample(p)
#' shuff_y_ind <- sample(q)
#' 
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' # Association pattern and effect sizes
#' #
#' pat <- matrix(FALSE, ncol = q, nrow = p)
#' bool_x <- shuff_x_ind <= p_act
#' bool_y <- shuff_y_ind <= q_act
#' 
#' pat_act <- beta_act <- matrix(0, nrow = p_act, ncol = q_act)
#' pat_act[sample(p_act * q_act, floor(p_act * q_act / 5))] <- 1
#' beta_act[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#' 
#' pat[bool_x, bool_y] <- pat_act
#'
#' # Gaussian responses
#' #
#' Y_act <- matrix(rnorm(n * q_act, mean = X_act %*% beta_act), nrow = n)
#' Y_inact <- matrix(rnorm(n * (q - q_act)), nrow = n)
#'
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]
#'
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' # Expectation and variance for the prior number of predictors associated with
#' # each response
#' #
#' p0 <- c(mean(colSums(pat)), 10)
#' 
#' #default, full CAVI
#' res_atlas <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "full")
#' 
#' #random CAVI
#' vb <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "random")
#' 
#' #AF CAVI, default is epsilon_scheme = iteration, start_partial = 50, geom_alpha = 0.95, min_epsilon = 0
#' vb <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "adaptive")
#' 
#' #other epsilon_scheme of AF CAVI: minimum, i.e., always set epsilon to a minimum value
#' vb <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "adaptive", epsilon_scheme = "iteration", min_epsilon = 0.01)
#' 
#' #other epsilon_scheme of AF CAVI: ELBO, i.e., epsilon is a logistic function of log ELBO difference
#' vb <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "adaptive", epsilon_scheme = "ELBO")
#'
#' @references
#' Helene Ruffieux, Anthony C. Davison, Jorg Hager, Jamie Inshaw, Benjamin P. 
#' Fairfax, Sylvia Richardson, Leonardo Bottolo, A global-local approach for 
#' detecting hotspots in multiple-response regression, arXiv:1811.03334, 2018.
#'
#' @seealso \code{\link{set_hyper}}, \code{\link{set_init}}.
#'
#' @export
#'
atlasqtl <- function(Y, X, p0, anneal = c(1, 2, 10), maxit = 1000, 
                     user_seed = NULL, verbose = 1, list_hyper = NULL, 
                     list_init = NULL, save_hyper = FALSE, save_init = FALSE, 
                     full_output = FALSE, thinned_elbo_eval = TRUE, 
                     checkpoint_path = NULL, trace_path = NULL, 
                     add_collinear_back = FALSE,
                     tol = 0.01,
                     CAVI = "full", start_partial=NULL,epsilon_scheme=NULL,min_epsilon = NULL,
                     geom_alpha = NULL,epsilon_param=NULL,fix_size=NULL,
                     eval_perform = F) {
  
  if (verbose != 0){
    cat(paste0("\n======================= \n",
               "== PREPROCESSING ... == \n",
               "======================= \n\n"))
  }
  
  check_verbose_(verbose)
  
  #check annealing setups
  if (verbose != 0) cat("== Checking the annealing schedule ... \n\n")
  
  check_annealing_(anneal)
  
  if (verbose != 0) cat("... done. == \n\n")
  
  #check partial-update scheme setups
  if (verbose != 0) cat("== Prepare the CAVI scheme parameters ... \n\n")
  
  config_CAVI <- parpare_partial_(CAVI,start_partial,
                                  epsilon_scheme,
                                  min_epsilon,
                                  geom_alpha,fix_size)
  
  if (verbose != 0) cat("... done. == \n\n")
  
  #prepare data
  if (verbose != 0) cat("== Preparing the data ... \n\n")
  
  dat <- prepare_data_(Y, X, tol, maxit, user_seed, verbose, checkpoint_path, 
                       trace_path)
  
  bool_rmvd_x <- dat$bool_rmvd_x
  
  X <- dat$X
  Y <- dat$Y
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  names_x <- colnames(X)
  names_y <- colnames(Y)
  
  shr_fac_inv <- q # = 1 / shrinkage_factor for global variance
  
  if (verbose != 0) cat("... done. == \n\n")
  
  
  if (verbose != 0) cat("== Preparing the hyperparameters ... \n\n")
  
  if (is.null(list_hyper) | is.null(list_init)) {
    
    check_structure_(p0, "vector", "numeric", 2)
    check_positive_(p0)
    
  } else {
    
    if (!is.null(p0))
      warning(paste0("Provided argument p0 not used, as both list_hyper ",
                     "and list_init were provided."))
    
  }
  
  list_hyper <- prepare_list_hyper_(list_hyper, Y, p, p0, bool_rmvd_x, 
                                    names_x, names_y, verbose)
  
  if (verbose != 0) cat("... done. == \n\n")
  
  
  if (verbose != 0) cat("== Preparing the parameter initialization ... \n\n")
  
  list_init <- prepare_list_init_(list_init, Y, p, p0, bool_rmvd_x, shr_fac_inv,
                                  user_seed, verbose)
  
  if (verbose != 0) cat("... done. == \n\n")
  
  
  if (verbose != 0){
    
    
    cat(paste0("**************************************************** \n",
               "Number of samples: ", n, "\n",
               "Number of (non-redundant) candidate predictors: ", p, "\n",
               "Number of responses: ", q, "\n",
               "**************************************************** \n\n"))
    
    cat(paste0("======================================================= \n",
               "== ATLASQTL: fast global-local hotspot QTL detection == \n",
               "======================================================= \n\n"))
  }
  
  
  hs <- TRUE
  debug <- TRUE
  
  if (hs) {
    
    df <- 1
    batch = "y" #use the C++ for the for loop
    
    res_atlas <- atlasqtl_global_local_core_(Y, X, shr_fac_inv, anneal, df,
                                             maxit, verbose, list_hyper,
                                             list_init, checkpoint_path,
                                             trace_path, full_output,
                                             thinned_elbo_eval, debug,
                                             batch,
                                             tol,
                                             config_CAVI,
                                             eval_perform)
    
  } else {
    
    if (!is.null(trace_path)) 
      warning(paste0("Provided argument trace_path not used, when using the ",
                     "global-scale-only model."))
    
    res_atlas <- atlasqtl_global_core_(Y, X, shr_fac_inv, anneal, df, tol, 
                                       maxit, verbose, list_hyper, list_init, 
                                       checkpoint_path, full_output, 
                                       thinned_elbo_eval, debug, 
                                       batch)
    
  }
  
  res_atlas$p0 <- p0
  
  res_atlas$rmvd_cst_x <- dat$rmvd_cst_x
  res_atlas$rmvd_coll_x <- dat$rmvd_coll_x
  
  if (add_collinear_back & length(res_atlas$rmvd_coll_x) > 0) {
    
    list_coll_back <- add_collinear_back_(res_atlas$beta_vb, 
                                          res_atlas$gam_vb, 
                                          res_atlas$theta_vb, 
                                          dat$initial_colnames_X, 
                                          res_atlas$rmvd_coll_x,
                                          verbose)
    
    res_atlas$beta_vb <- list_coll_back$beta_vb
    res_atlas$gam_vb <- list_coll_back$gam_vb
    res_atlas$theta_vb <- list_coll_back$theta_vb
    
  }
  
  if (save_hyper) res_atlas$list_hyper <- list_hyper
  if (save_init) res_atlas$list_init <- list_init
  
  class(res_atlas) <- "atlasqtl"
  
  if (verbose != 0) cat("... done. == \n\n")
  
  res_atlas
  
}
