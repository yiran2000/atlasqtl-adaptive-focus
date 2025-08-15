rm(list = ls())

set.seed(123)

############################
## simulate basic dataset ##
############################

n <- 100; p <- 75; q <- 20; p_act <- 10

# candidate predictors (subject to selection)
X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.2), nrow = n)
X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.2), nrow = n)
X <- cbind(X_act, X_inact)[, sample(p)]

beta <-  matrix(rnorm(p_act * q), nrow = p_act)

# Gaussian outcomes
Y <- matrix(rnorm(n * q, mean = X_act %*% beta, sd = 1), nrow = n)

# remove constant variables (needed for checking dimension consistency)
X <- scale(X)
rm_cst <- function(mat_sc) mat_sc[, !is.nan(colSums(mat_sc))]
rm_coll <- function(mat_sc) mat_sc[, !duplicated(mat_sc, MARGIN = 2)]

X <- rm_cst(X)
X <- rm_coll(X)

p <- ncol(X)


########################
## atlasqtl inference ##
########################

p0 <- c(5, 25)

# Continuous outcomes, no covariates
#
#default, full CAVI
res_atlas <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "full")

#random CAVI
vb <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "random")

#AF CAVI, default is epsilon_scheme = iteration, start_partial = 50, geom_alpha = 0.95, min_epsilon = 0
vb <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "adaptive")

#other epsilon_scheme of AF CAVI: minimum, i.e., always set epsilon to a minimum value
vb <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "adaptive", epsilon_scheme = "iteration", min_epsilon = 0.01)

#other epsilon_scheme of AF CAVI: ELBO, i.e., epsilon is a logistic function of log ELBO difference
vb <- atlasqtl(Y = Y, X = X, p0 = p0, CAVI = "adaptive", epsilon_scheme = "ELBO")

