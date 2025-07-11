# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#
# Internal core function to call the variational algorithm for global-local 
# hotspot propensity modelling. 
# See help of `atlasqtl` function for details.
#
library(tictoc)

# lognormal_cdf <- function(x, mu, sigma, m) {
#   return(m + (1 - m) * pnorm((log(x) - mu) / sigma))
# }
# emax <- function(x, e0=0, emax = 1, ec50, n) {
#   e0 + ((emax-e0) * x^n) / (ec50^n +x^n)
# }

# update e:
# e = lognormal_cdf(diff_lb, mu=epsilon[1], sigma=epsilon[2], m=epsilon[3])
# e = 1- tanh_function(log(it)/2)

logistic_function = function(x) {
  1/ (1 + exp(-x))
}

# tanh_function <- function(x) {
#   (exp(x) - exp(-x)) / (exp(x) + exp(-x))
# }



atlasqtl_global_local_core_ <- function(Y, X, shr_fac_inv, anneal, df, 
                                        maxit, verbose, list_hyper, list_init, 
                                        checkpoint_path = NULL, 
                                        trace_path = NULL, full_output = FALSE, 
                                        thinned_elbo_eval = TRUE,
                                        debug = FALSE, 
                                        batch, 
                                        tol,
                                        config_CAVI,
                                        eval_perform) {
  
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Y)
  
  #devel specification of partial update that we don't use anymore
  max_partial = NULL #maximum number of iterations for partial-update, if one wants to add
  end_full = F #whether the algorithm should go back to full in the end
  ifZ = F #whether Z is also partial 
  ifTau = F #whether tau is also partial 
  
  #specifications of partial update
  
  scheme = config_CAVI$scheme
  start_partial = config_CAVI$start_partial
  epsilon_scheme = config_CAVI$epsilon_scheme
  min_epsilon = config_CAVI$min_epsilon
  geom_alpha = config_CAVI$geom_alpha
  fix_size = config_CAVI$fix_size
  
  
  # Pre-computing large matrices
  # 
  tic("precompute")
  if (any(is.na(Y))) {
    
    mis_pat <- ifelse(is.na(Y), 0, 1)
    Y[is.na(Y)] <- 0
    X_norm_sq <- crossprod(X^2, mis_pat)
    
    cp_X_rm <- lapply(1:q, function(k) {
      if (any(mis_pat[,k] == 0)) {
        ind <- which(mis_pat[,k] == 0)
        crossprod(X[ind,, drop = FALSE])
      } else {
        matrix(0, nrow = p, ncol = p)
      }
    })
    
  } else {
    mis_pat <- X_norm_sq <- cp_X_rm <- NULL
  }
  
  Y_norm_sq <- colSums(Y^2) # must be after the if as some entries of y set to 0 when missing values
  cp_X <- crossprod(X)
  cp_Y_X <- crossprod(Y, X)
  t_init = toc(quiet = TRUE)
  time_init = t_init$toc - t_init$tic
  

  # Gathering initial variational parameters. Do it explicitly.
  # with() function not used here as the objects will be modified later.
  #
  gam_vb <- list_init$gam_vb
  mu_beta_vb <- list_init$mu_beta_vb # mu_beta_vb[s, t] = variational posterior mean of beta_st | gamma_st = 1 
  # (beta_vb[s, t] = varational posterior mean of beta_st)
  sig02_inv_vb <- list_init$sig02_inv_vb  # horseshoe global precision parameter
  sig2_beta_vb <- list_init$sig2_beta_vb  # sig2_beta_vb[t]  = variational variance of beta_st | gamma_st = 1 
  # (same for all s as X is standardized)
  sig2_theta_vb <- list_init$sig2_theta_vb 
  tau_vb <-list_init$tau_vb
  theta_vb <- list_init$theta_vb
  zeta_vb <- list_init$zeta_vb
  
  rm(list_init)
  
  theta_plus_zeta_vb <- sweep(tcrossprod(theta_vb, rep(1, q)), 2, zeta_vb, `+`)
  log_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE)
  log_1_min_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE, lower.tail = FALSE)
  
  
  # Preparing trace saving if any
  #
  trace_ind_max <- trace_var_max <- NULL
  
  # Preparing annealing if any
  #
  anneal_scale <- TRUE
  
  if (is.null(anneal)) {
    annealing <- FALSE
    c_cst <- c_s <- 1 # c_s for scale parameters
    it_endAnneal <- 1 # first non-annealed iteration 
  } else {
    annealing <- TRUE
    ladder <- get_annealing_ladder_(anneal, verbose)
    c_cst <- ladder[1]
    c_s <- ifelse(anneal_scale, c_cst, 1) 
    it_endAnneal <- anneal[3] # first non-annealed iteration 
  }
  
  eps <- .Machine$double.eps^0.5
  
  if(thinned_elbo_eval){
    times_conv_sched <- c(1, 2, 5, 10, 100)
    batch_conv_sched <- c(1, 5, 10, 25, 50)
  }else{
    times_conv_sched <- 1
    batch_conv_sched <- 1
  }
 
  ind_batch_conv <- length(batch_conv_sched) + 1 # so that, the first time, it enters in the loop below 
  batch_conv <- 1 
  

  with(list_hyper, { # list_init not used with the with() function to avoid
    # copy-on-write for large objects
    # Response-specific parameters: objects derived from t02
    #
    t02_inv <- 1 / t02
    sig2_zeta_vb <- update_sig2_c0_vb_(p, t02, c = c_cst) # stands for a diagonal matrix of size d with this value on the (constant) diagonal
    
    vec_sum_log_det_zeta <- - q * (log(t02) + log(p + t02_inv))
    
    
    # Stored/precomputed objects
    #
    beta_vb <- update_beta_vb_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE) # first time keep sweep = TRUE even when missing data, since uses the initial parameter sig2_beta_vb which is a vector.
    

    cp_X_Xbeta <- update_cp_X_Xbeta_(cp_X, beta_vb, cp_X_rm)
    

    
    # Fixed VB parameter
    #
    nu_xi_inv_vb <- 1 # no change with annealing 
    
    # Initialize ELBO, iterations, error terms
    #
    lb_new <- -Inf #the latest ELBO
    it <- 0 #total iteration index
    it_0 = 0 #partial iteration index
    diff_lb = Inf #current difference of ELBO
    #define the error term in response selection
    e = 1
    
    # Initialize algorithm status
    #
    converged <- FALSE #marks whether converged
    # init = T #marker whether we are in the initial full update stage
    partial = F #marks whether in the partial-update stage or convergence evaluation stage
    
    # Initialize empty vectors to keep track of the algorithm
    #
    if(eval_perform){
      partial_ls = list() #list of T/F, whether in the partial-update stage of not
      subsample_ls = list() #list of T/F, whether subsampling in the convergence evaluation stage
      ELBO_ls = list() #ELBO
      ELBO_diff_ls = list() #difference of ELBO in two consecutive iterations
      it_ls = list() #iteration number
      e_ls = list() #error term in response selection\
      annealing_ls = list()
      
      time_loop_ls = list() #runtime of CoreDualLoop
      time_total_ls = list() #total runtime of one iteration
      time_sigma_ls = list()
      time_tau_ls = list()
      time_sig2_beta_ls = list()
      time_m2_beta_ls = list()
      time_Z_ls = list()
      time_theta_ls = list()
      time_zeta_ls = list()
      time_thetaPzeta_ls = list()
      time_ELBO_ls = list()
      
      subsample_size_ls = list() #number of responses selected in the iteration
      
      # it_eval_ls = list() #iterations where ELBO is evaluated
      # zeta_ls = list() #zeta
      # select_prob_ls = list()
      # r_vc_ls = list()
      
    }
    
    # if(ifZ){
    #   Z <- update_Z_(gam_vb, theta_plus_zeta_vb, log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, c = c)
    # }
    
    # CAVI starts
    #
    while ((!converged) & (it < maxit)) {
      

      t0 = Sys.time()
      
      lb_old <- lb_new
      it <- it + 1
      
      start_partial_iter = T
      
      if(start_partial_iter){ #starting partial update based on iterations, the optimal choice
        if(start_partial != Inf){
          if(it < start_partial & partial == F){
            partial = F
          }else{
            partial = T
            it_0 = it_0 + 1
          }
        }
        
      }else{ #starting partial update based on ELBO, do not use it anymore
        if(start_partial != Inf){
          if(diff_lb > start_partial*tol & partial == F){
            partial = F
          }else{
            partial = T
            it_0 = it_0 + 1
          }
        }
      }
      

      
      if (verbose != 0 &  (it == 1 | it %% max(5, batch_conv) == 0)) 
        cat(paste0("Iteration ", format(it), "... \n"))
      
      
      # generate subsample
      if(partial){
        
        #calculate the selection probability, with or without error
        if(!is.null(epsilon_scheme)){
          #in the adaptive CAVI scheme, calculate the error term epsilon for selection
          if(epsilon_scheme == "ELBO"){
            e = max(min_epsilon, logistic_function(log(diff_lb)))
          }else if(epsilon_scheme == "iteration"){
            if(geom_alpha == 0){
              e = 0
            }else{
              e <- max(min_epsilon, geom_alpha^(it_0-1))
            }
          }else if (epsilon_scheme == "minimum"){ #the minimum error scheme
            e = min_epsilon
          }
          
          r_vc = 1- apply((1 - gam_vb), 2, prod) #PPI
          select_prob =  (1 - e)*r_vc + e #probability of selecting each response by adding the error
          
        }else{
          #in the random selection scheme, select_prob is just 0.5
          select_prob = rep(0.5, q)
        }

        
        if(batch == "y"){
          if(is.null(fix_size)){ #fix_size is only for the adaptive scheme, whether selecting a fixed number or based on Bernoulli 
            sample_q = c(0:(q-1))[rbinom(q, size = 1, prob = select_prob) == 1]
          }else{
            sample_q = sample(c(0:(q-1)), size = fix_size*q, replace = FALSE, prob = select_prob)
          }
        }else{
          sample_q = c(1:q)[rbinom(q, size = 1, prob = select_prob) == 1]
        }
        
        # select_prob_ls = append(select_prob_ls, list(select_prob))
        
      }else{
        if(batch == "y"){
          sample_q = sample(0:(q-1))
        }else{
          sample_q = sample(1:q)
        }
        
        # select_prob_ls = append(select_prob_ls, list(rep(1, q)))
      }
      
      # r_vc_ls = append(r_vc_ls, list(1- apply((1 - gam_vb), 2, prod)))
      
      #record partial and subsample_q
      if(eval_perform){
        partial_ls = c(partial_ls, partial)
        annealing_ls = c(annealing_ls, annealing)
      }
      
      # update VB parameters
      tic("sigma")
      # % # sigma^2, global 
      nu_vb <- update_nu_vb_(nu, sum(gam_vb), c = c_cst)
      rho_vb <- update_rho_vb_(rho, m2_beta, tau_vb, c = c_cst)
      
      sig2_inv_vb <- nu_vb / rho_vb
      # % #
      t_sigma = toc(quiet = TRUE)
      time_sigma = t_sigma$toc - t_sigma$tic
      
      
      tic("tau")
      # % # tau, local
      if(ifTau & partial & it > 1){
        eta_vb[sample_q] <- update_eta_vb_(n, eta[sample_q], gam_vb[,sample_q], mis_pat, c = c_cst)
        kappa_vb[sample_q] <- update_kappa_vb_partial_(n, Y_norm_sq, cp_Y_X, cp_X_Xbeta, kappa, 
                                                       beta_vb, m2_beta, sig2_inv_vb, X_norm_sq, sample_q, c = c_cst)
        
        tau_vb[sample_q] <- eta_vb[sample_q] / kappa_vb[sample_q]
        log_tau_vb[sample_q] <- update_log_tau_vb_(eta_vb[sample_q], kappa_vb[sample_q])
      }else{
        eta_vb <- update_eta_vb_(n, eta, gam_vb, mis_pat, c = c_cst)
        kappa_vb <- update_kappa_vb_(n, Y_norm_sq, cp_Y_X, cp_X_Xbeta, kappa, 
                                     beta_vb, m2_beta, sig2_inv_vb, X_norm_sq, c = c_cst)
        
        tau_vb <- eta_vb / kappa_vb
        log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
        log_sig2_inv_vb <- update_log_sig2_inv_vb_(nu_vb, rho_vb)
      }
      
      t_tau = toc(quiet = TRUE)
      time_tau = t_tau$toc - t_tau$tic
      
      
      #beta
      tic("sig2_beta")
      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb, tau_vb, X_norm_sq, c = c_cst) 
      t_sig2_beta = toc(quiet = TRUE)
      time_sig2_beta = t_sig2_beta$toc - t_sig2_beta$tic
      

      
      # different possible batch-coordinate ascent schemes:
      #
      
      if (batch == "y") { # optimal scheme
        
        # C++ Eigen call for expensive updates 
        #
        # Shuffle updates order at each iteration
        #
        shuffled_ind <- as.integer(sample(0:(p-1))) # Zero-based index in C++
        
        if (is.null(mis_pat)) {
          
          tic("for loop")

          coreDualLoop(cp_X, cp_Y_X, gam_vb, log_Phi_theta_plus_zeta,
                       log_1_min_Phi_theta_plus_zeta, log_sig2_inv_vb, log_tau_vb,
                       beta_vb, cp_X_Xbeta, mu_beta_vb, sig2_beta_vb, tau_vb,
                       shuffled_ind, sample_q, c = c_cst)


          time = toc(quiet = TRUE)
          t = time$toc - time$tic
          
        } else {
          tic("for loop")
          coreDualMisLoop(cp_X, cp_X_rm, cp_Y_X, gam_vb, log_Phi_theta_plus_zeta, 
                          log_1_min_Phi_theta_plus_zeta, log_sig2_inv_vb, log_tau_vb, 
                          beta_vb, cp_X_Xbeta, mu_beta_vb, sig2_beta_vb, tau_vb, 
                          shuffled_ind, sample_q, c = c_cst)
          time = toc(quiet = TRUE)
          t = time$toc - time$tic
          
          
        }
        
        
      } else if (batch == "0"){ # no batch, used only internally (slower)
        # schemes "x" of "x-y" are not batch concave
        # hence not implemented as they may diverge
        
        for (k in sample_q) {
          
          
          if (is.null(mis_pat)) {
            
            for (j in sample(1:p)) {
              
              cp_X_Xbeta[k, ] <- cp_X_Xbeta[k, ] - beta_vb[j, k] * cp_X[j, ]
              
              mu_beta_vb[j, k] <- c * sig2_beta_vb[k] * tau_vb[k] * (cp_Y_X[k, j] - cp_X_Xbeta[k, j])
              
              gam_vb[j, k] <- exp(-log_one_plus_exp_(c * (pnorm(theta_vb[j] + zeta_vb[k], lower.tail = FALSE, log.p = TRUE) -
                                                            pnorm(theta_vb[j] + zeta_vb[k], log.p = TRUE) -
                                                            log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                            mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                            log(sig2_beta_vb[k]) / 2)))
              
              beta_vb[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]

              cp_X_Xbeta[k, ] <- cp_X_Xbeta[k, ] + beta_vb[j, k] * cp_X[j, ]

              
            }
            
          } else {
            
            for (j in sample(1:p)) {

              cp_X_Xbeta[k, ] <- cp_X_Xbeta[k, ] - beta_vb[j, k] * (cp_X[j, ] - cp_X_rm[[k]][j, ])

              mu_beta_vb[j, k] <- c * sig2_beta_vb[j, k] * tau_vb[k] * (cp_Y_X[k, j] - cp_X_Xbeta[k, j])

              gam_vb[j, k] <- exp(-log_one_plus_exp_(c * (pnorm(theta_vb[j] + zeta_vb[k], lower.tail = FALSE, log.p = TRUE) -
                                                            pnorm(theta_vb[j] + zeta_vb[k], log.p = TRUE) -
                                                            log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                            mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[j, k]) -
                                                            log(sig2_beta_vb[j, k]) / 2)))

              beta_vb[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]

              cp_X_Xbeta[k, ] <- cp_X_Xbeta[k, ] + beta_vb[j, k] * (cp_X[j, ] - cp_X_rm[[k]][j, ])
            }
            
          }
        }
        
      } else {
        
        stop ("Batch scheme not defined. Exit.")
        
      }
      
      tic("beta2")
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, mis_pat = mis_pat)
      t_m2_beta = toc(quiet = TRUE)
      time_m2_beta = t_m2_beta$toc - t_m2_beta$tic

      tic("Z")
      if(ifZ==T & it > 1 & partial){
        Z <- update_Z_partial_(Z, gam_vb, theta_plus_zeta_vb, log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, sample_q, c = c_cst)
      }else{
        Z <- update_Z_(gam_vb, theta_plus_zeta_vb, log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, c = c_cst)
      }
      t_Z = toc(quiet = TRUE)
      time_Z = t_Z$toc - t_Z$tic
      
      
      tic("theta")
      # keep this order!
      #  
      L_vb <- c_s * sig02_inv_vb * shr_fac_inv * (theta_vb^2 + sig2_theta_vb - 2 * theta_vb * m0 + m0^2) / 2 / df 
      rho_xi_inv_vb <- c_s * (A2_inv + sig02_inv_vb) 
      
      if (annealing & anneal_scale) {
        
        lam2_inv_vb <- update_annealed_lam2_inv_vb_(L_vb, c_s, df)
        
      } else {
        
        Q_app <- Q_approx_vec(L_vb)
        
        if (df == 1) {
          
          lam2_inv_vb <- 1 / (Q_app * L_vb) - 1
          
        } else if (df == 3) {
          
          lam2_inv_vb <- exp(-log(3) - log(L_vb) + log(1 - L_vb * Q_app) - log(Q_app * (1 + L_vb) - 1)) - 1 / 3
          
        } else {
          # also works for df = 3 but might be slightly less efficient than the above
          
          exponent <- (df + 1) / 2
          
          lam2_inv_vb <- sapply(1:p, function(j) {
            
            exp(log(compute_integral_hs_(df, L_vb[j] * df, m = exponent, n = exponent, Q_ab = Q_app[j])) -
                  log(compute_integral_hs_(df, L_vb[j] * df, m = exponent, n = exponent - 1, Q_ab = Q_app[j])))
            
          })
          
        }
        
      }
      
      xi_inv_vb <- nu_xi_inv_vb / rho_xi_inv_vb
      
      sig2_theta_vb <- update_sig2_c0_vb_(q, 1 / (sig02_inv_vb * lam2_inv_vb * shr_fac_inv), c = c_cst)
      
      theta_vb <- update_theta_vb_(Z, m0, sig02_inv_vb * lam2_inv_vb * shr_fac_inv, sig2_theta_vb,
                                   vec_fac_st = NULL, zeta_vb, is_mat = FALSE, c = c_cst)
      
      t_theta = toc(quiet = TRUE)
      time_theta = t_theta$toc - t_theta$tic
      
      tic("zeta")
      nu_s0_vb <- update_nu_vb_(1 / 2, p, c = c_s)
      
      rho_s0_vb <- c_s * (xi_inv_vb + 
                            sum(lam2_inv_vb * shr_fac_inv * (theta_vb^2 + sig2_theta_vb - 2 * theta_vb * m0 + m0^2)) / 2) 
      
      sig02_inv_vb <- as.numeric(nu_s0_vb / rho_s0_vb)
      
      zeta_vb <- update_zeta_vb_(Z, theta_vb, n0, sig2_zeta_vb, t02_inv,
                                 is_mat = FALSE, c = c_cst) 

      t_zeta = toc(quiet = TRUE)
      time_zeta = t_zeta$toc - t_zeta$tic
      
      # #save all zeta_vb to visualize
      # zeta_ls = append(zeta_ls, list(zeta_vb))
      
      tic("theta+zeta")
      
      theta_plus_zeta_vb <- sweep(tcrossprod(theta_vb, rep(1, q)), 2, zeta_vb, `+`)
      log_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE)
      log_1_min_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE, lower.tail = FALSE)  
      
      t_thetaPzeta = toc(quiet = TRUE)
      time_thetaPzeta = t_thetaPzeta$toc - t_thetaPzeta$tic
      
      
      
      if (verbose == 2 && (it == 1 | it %% max(5, batch_conv) == 0)) {
        
        cat(paste0("Variational hotspot propensity global scale: ", 
                   format(sqrt(rho_s0_vb / (nu_s0_vb - 1) / shr_fac_inv), digits = 3), ".\n"))
        cat("Approximate variational hotspot propensity local scale: \n")
        # print(summary(sqrt(1 / lam2_inv_vb)))
        cat("\n")
        
      }
      
      if (!is.null(trace_path) && (it == 1 | it %% 25 == 0)) {
        
        list_traces <- plot_trace_var_hs_(lam2_inv_vb, sig02_inv_vb, q, it, 
                                          trace_ind_max, trace_var_max, trace_path)
        
        trace_ind_max <- list_traces$trace_ind_max
        trace_var_max <- list_traces$trace_var_max
        
      }
      
      tic("ELBO")
      # Update annealing parameters
      #
      if (annealing) {
        
        if (verbose != 0 & (it == 1 | it %% 5 == 0))
          cat(paste0("Temperature = ", format(1 / c_cst, digits = 4), "\n\n"))
        
        sig2_zeta_vb <- c_cst * sig2_zeta_vb
        
        c_cst <- ifelse(it < length(ladder), ladder[it + 1], 1)
        c_s <- ifelse(anneal_scale, c_cst, 1)
        
        sig2_zeta_vb <- sig2_zeta_vb / c_cst
        
        if (isTRUE(all.equal(c_cst, 1))) {
          
          annealing <- FALSE
          
          if (verbose != 0)
            cat("** Exiting annealing mode. **\n\n")
        }
        
        
      } else {
        
      # Evaluate ELBO
      if(partial){it_e = it_0}else{it_e = it}
      if (it <= it_endAnneal + 1 | it_e %% batch_conv == 0 | it_e %% batch_conv == 1 ){ 
        # it <= it_endAnneal + 1 evaluate the ELBO for the first two non-annealed iterations
        # to (also) evaluate convergence between two consecutive iterations

        lb_new <- elbo_global_local_full_(Y, sample_q, partial, ifTau, A2_inv, beta_vb, df, eta, eta_vb, gam_vb,
                                          kappa, kappa_vb, L_vb, lam2_inv_vb, log_tau_vb,
                                          log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, 
                                          m0, m2_beta, n0, nu, nu_s0_vb, nu_vb, nu_xi_inv_vb,
                                          Q_app, rho, rho_s0_vb, rho_vb, rho_xi_inv_vb,
                                          shr_fac_inv, sig02_inv_vb, sig2_beta_vb,
                                          sig2_inv_vb,  sig2_theta_vb, sig2_zeta_vb,
                                          t02_inv, tau_vb, theta_vb, vec_sum_log_det_zeta,
                                          xi_inv_vb, zeta_vb, X_norm_sq, Y_norm_sq, cp_Y_X, cp_X_Xbeta, mis_pat)



        
        if (verbose != 0 & (it == it_endAnneal | it %% max(5, batch_conv) == 0))
          cat(paste0("ELBO = ", format(lb_new), "\n\n"))
        
        if (debug && lb_new + eps < lb_old){
          stop("ELBO not increasing monotonically. Exit. ")
        }

        # diff_lb = abs(lb_new - lb_old)
        diff_lb = lb_new - lb_old
        
        # # Record the iteration where ELBO is evaluated
        # if(eval_perform){
        #   it_eval_ls = c(it_eval_ls, it)
        # }


        sum_exceed <- sum(diff_lb > (times_conv_sched * tol))
        # times_conv_sched*tol sets how many more iterations should be conducted before the next ELBO evaluaion, 
        # which is defined by batch_conv_shed
        
        if (sum_exceed == 0){
          
          if (end_full & partial){ #switch back to full update if end_full
              partial = F
            } else if (end_full == F & partial){
              converged = T
            } else{
              converged = T
            }
          
        }else if (ind_batch_conv > sum_exceed) {
          
          if(partial){
            # If ELBO diff is too large, set the next time where ELBO is evaluated depending on how big the difference is
            ind_batch_conv <- sum_exceed
            batch_conv <- batch_conv_sched[ind_batch_conv]
          }else{
            ind_batch_conv = length(batch_conv_sched) + 1 
            batch_conv = 1
          }


        }
      
      }
        

        
        checkpoint_(it, checkpoint_path, beta_vb, gam_vb, theta_vb, zeta_vb, 
                    converged, lb_new, lb_old,
                    lam2_inv_vb = lam2_inv_vb, sig02_inv_vb = sig02_inv_vb,
                    names_x = colnames(X), names_y = colnames(Y))
          
          
      }
      t_ELBO = toc(quiet = TRUE)
      time_ELBO = t_ELBO$toc - t_ELBO$tic


      if(!is.null(max_partial)){
        if(end_full & partial & it_0 >= max_partial){
          partial = F
        }
      }
      
      
      #run time of the entire iteration
      t1 = Sys.time()-t0
      
      if(eval_perform){
        it_ls = c(it_ls, it)
        ELBO_ls = c(ELBO_ls, lb_new)
        e_ls = c(e_ls, e)
        ELBO_diff_ls = c(ELBO_diff_ls, diff_lb)
        

        time_loop_ls = c(time_loop_ls, t)
        time_total_ls = c(time_total_ls, t1)
        time_sigma_ls = c(time_sigma_ls, time_sigma)
        time_tau_ls = c(time_tau_ls, time_tau)
        time_sig2_beta_ls = c(time_sig2_beta_ls, time_sig2_beta)
        time_m2_beta_ls = c(time_m2_beta_ls, time_m2_beta)
        time_Z_ls = c(time_Z_ls, time_Z)
        time_theta_ls = c(time_theta_ls, time_theta)
        time_zeta_ls = c(time_zeta_ls, time_zeta)
        time_thetaPzeta_ls = c(time_thetaPzeta_ls, time_thetaPzeta)
        time_ELBO_ls = c(time_ELBO_ls, time_ELBO)
        
        subsample_size_ls = c(subsample_size_ls, length(sample_q))
      }
      
    }
    
    
    checkpoint_clean_up_(checkpoint_path)
    
    
    if (verbose != 0) {
      if (converged) {
        cat(paste0("Convergence obtained after ", format(it), " iterations. \n",
                   "Optimal marginal log-likelihood variational lower bound ",
                   "(ELBO) = ", format(lb_new), ". \n\n"))
      } else {
        warning("Maximal number of iterations reached before convergence. Exit.")
      }
    }
    
    lb_opt <- lb_new
    s02_vb <- rho_s0_vb / (nu_s0_vb - 1) / shr_fac_inv # inverse gamma mean, with embedded shrinkage
    
    if (full_output) { # for internal use only
      
      create_named_list_(beta_vb, eta_vb, gam_vb, kappa_vb, lam2_inv_vb, 
                         nu_s0_vb, nu_vb, nu_xi_inv_vb, rho_s0_vb, rho_vb, 
                         rho_xi_inv_vb, shr_fac_inv, sig02_inv_vb, sig2_beta_vb, 
                         sig2_inv_vb,  sig2_theta_vb, sig2_zeta_vb, tau_vb, 
                         theta_vb, cp_Y_X, cp_X, cp_X_Xbeta, xi_inv_vb, zeta_v)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      names_n <- rownames(Y)
      
      rownames(gam_vb) <- rownames(beta_vb) <- names_x
      colnames(gam_vb) <- colnames(beta_vb) <- names_y
      names(theta_vb) <- names_x
      names(zeta_vb) <- names_y
      names(lam2_inv_vb) <- names_x
      
      # diff_lb <- abs(lb_opt - lb_old)
      diff_lb <- lb_opt - lb_old
      
      if(eval_perform){
        perform_df = data.frame(
          iter = it_ls %>% unlist,
          partial = partial_ls %>% unlist,
          ELBO = ELBO_ls %>% unlist,
          ELBO_diff = ELBO_diff_ls %>% unlist,
          e = e_ls %>% unlist,
          
          time_loop = time_loop_ls %>% unlist,
          time_total = time_total_ls %>% unlist,
          time_sigma = time_sigma_ls %>% unlist,
          time_tau = time_tau_ls %>% unlist,
          time_sig2_beta =time_sig2_beta_ls %>% unlist,
          time_m2_beta = time_m2_beta_ls %>% unlist,
          time_Z = time_Z_ls %>% unlist,
          time_theta =time_theta_ls %>% unlist,
          time_zeta = time_zeta_ls %>% unlist,
          time_thetaPzeta = time_thetaPzeta_ls %>% unlist,
          time_ELBO = time_ELBO_ls %>% unlist,
          
          subsample_size = subsample_size_ls %>% unlist,
          annealing = annealing_ls %>% unlist
        )
        
        create_named_list_(beta_vb, gam_vb, theta_vb, zeta_vb, 
                           n, p, q, anneal, converged, it, maxit, tol, lb_opt, 
                           diff_lb,
                           perform_df,time_init)
                           # it_eval_ls,
                           # zeta_ls,
                           # select_prob_ls,
                           # r_vc_ls,

        
      } else {
        create_named_list_(beta_vb, gam_vb, theta_vb, zeta_vb, 
                           n, p, q, anneal, converged, it, maxit, tol, lb_opt, 
                           diff_lb)
        
      }
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `atlasqtl_struct_core` algorithm.
#
elbo_global_local_full_ <- function(Y, sample_q, partial, ifTau, A2_inv, beta_vb, df, eta, eta_vb, gam_vb, 
                               kappa, kappa_vb, L_vb, lam2_inv_vb, log_tau_vb,
                               log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, m0, m2_beta, 
                               n0, nu, nu_s0_vb, nu_vb, nu_xi_inv_vb, 
                               Q_app, rho, rho_s0_vb, rho_vb, rho_xi_inv_vb, 
                               shr_fac_inv, sig02_inv_vb, sig2_beta_vb, 
                               sig2_inv_vb,  sig2_theta_vb, sig2_zeta_vb, 
                               t02_inv, tau_vb, theta_vb, vec_sum_log_det_zeta, 
                               xi_inv_vb, zeta_vb, X_norm_sq, Y_norm_sq, cp_Y_X, 
                               cp_X_Xbeta, mis_pat) {
  
  n <- nrow(Y)
  p <- length(L_vb)
  
  # needed for monotonically increasing elbo.
  #
  # if(ifTau & partial){
  #   eta_vb[sample_q] <- update_eta_vb_(n, eta[sample_q], gam_vb[,sample_q], mis_pat)
  #   kappa_vb[sample_q] <- update_kappa_vb_partial_(n, Y_norm_sq, cp_Y_X, cp_X_Xbeta, kappa, 
  #                                                  beta_vb, m2_beta, sig2_inv_vb, X_norm_sq, sample_q)
  #   
  #   log_tau_vb[sample_q] <- update_log_tau_vb_(eta_vb[sample_q], kappa_vb[sample_q])
  # }else{
    eta_vb <- update_eta_vb_(n, eta, gam_vb, mis_pat)
    kappa_vb <- update_kappa_vb_(n, Y_norm_sq, cp_Y_X, cp_X_Xbeta, kappa, beta_vb, 
                                 m2_beta, sig2_inv_vb, X_norm_sq)
    log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  # }
  
  
  nu_vb <- update_nu_vb_(nu, sum(gam_vb))
  rho_vb <- update_rho_vb_(rho, m2_beta, tau_vb)
  
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(nu_vb, rho_vb)
  
  log_sig02_inv_vb <- update_log_sig2_inv_vb_(nu_s0_vb, rho_s0_vb)
  log_xi_inv_vb <- update_log_sig2_inv_vb_(nu_xi_inv_vb, rho_xi_inv_vb)
  
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb, mis_pat)
  
  elbo_B <- e_beta_gamma_(gam_vb, log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, log_sig2_inv_vb, 
                          log_tau_vb, zeta_vb, 
                          theta_vb, m2_beta, sig2_beta_vb, sig2_zeta_vb,
                          sig2_theta_vb, sig2_inv_vb, tau_vb)
  
  elbo_C <- e_theta_hs_(lam2_inv_vb, L_vb, log_sig02_inv_vb + log(shr_fac_inv), 
                        m0, theta_vb, Q_app, sig02_inv_vb * shr_fac_inv, 
                        sig2_theta_vb, df)
  
  elbo_D <- e_zeta_(zeta_vb, n0, sig2_zeta_vb, t02_inv, vec_sum_log_det_zeta)
  
  elbo_E <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_F <- e_sig2_inv_hs_(xi_inv_vb, nu_s0_vb, log_xi_inv_vb, log_sig02_inv_vb, 
                           rho_s0_vb, sig02_inv_vb) # sig02_inv_vb
  
  elbo_G <- e_sig2_inv_(1 / 2, nu_xi_inv_vb, log_xi_inv_vb, A2_inv, 
                        rho_xi_inv_vb, xi_inv_vb) # xi_inv_vb
  
  elbo_H <- e_sig2_inv_(nu, nu_vb, log_sig2_inv_vb, rho, rho_vb, sig2_inv_vb)
  
  as.numeric(elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G + elbo_H)
  
}
