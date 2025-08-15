/*
 * This file is part of the `AFatlasqtl` R package:
 *      https://github.com/yiran2000/atlasqtl-adaptive-focus
 *      
 * Functions for computationally expensive updates in algorithms without
 * external information.
 *
 * These functions use Eigen::Map to pass large matrices by reference from R.
 * Given dimensionalities involved in some applications, copying such matrices
 * would imply a prohibitive RAM overconsumption.
 *
 */

#include "utils.h"
/*
 double crossprod(const Eigen::VectorXd &one, const Eigen::VectorXd &two) {
 // R Crossproduct of vectors = sum of element-wise products
 // Assume one.size() == two.size() without checking
 double sum = 0;
 for (int i = 0; i < one.size(); i++) {
 sum += one[i] * two[i];
 }
 return sum;
 }
 */

double logOnePlusExp(double x) {
  double m = x;
  if (x < 0)
    m = 0;
  return log(exp(x-m) + exp(-m)) + m;
}

// for atlasqtl_core function
// [[Rcpp::export]]
void coreDualLoop(const MapMat cp_X,
                  const MapMat cp_Y_X,
                  MapArr2D gam_vb,
                  const MapArr2D log_Phi_theta_plus_zeta,
                  const MapArr2D log_1_min_Phi_theta_plus_zeta,
                  const double log_sig2_inv_vb,
                  const MapArr1D log_tau_vb,
                  MapMat m1_beta,
                  MapMat cp_betaX_X,
                  MapArr2D mu_beta_vb,
                  const MapArr1D sig2_beta_vb,
                  const MapArr1D tau_vb,
                  const Eigen::VectorXi shuffled_ind,
                  const Eigen::VectorXi sample_q,
                  const double c = 1) {
  
  // Rcout << "one\n";
  
  const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;
  
  for (int a = 0; a < sample_q.size(); a++) {
    int k = sample_q[a];
    
    for (int b = 0; b < shuffled_ind.size(); b++) {
      int j = shuffled_ind[b];
      
      
      double m1_beta_jk = m1_beta(j, k);
      //double cp_betaX_X_kj = cp_betaX_X(k, j) - m1_beta_jk * cp_X(j,j);
      double cp_betaX_X_jk = cp_betaX_X(j, k) - m1_beta_jk * cp_X(j,j);
      
      mu_beta_vb(j, k) = c * sig2_beta_vb[k] * tau_vb[k] * (cp_Y_X(k, j) - cp_betaX_X_jk);
      
      gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) - log_Phi_theta_plus_zeta(j, k)
                                               - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb[k])
                                               + cst[k])));
                                               
                                               m1_beta(j, k) = gam_vb(j, k) * mu_beta_vb(j, k);
                                               
                                               cp_betaX_X.col(k) += (m1_beta(j, k) - m1_beta_jk) * cp_X.col(j);
                                               
                                               
    } 
  }
}

// // [[Rcpp::export]]
// void coreDualLoopZ(const MapMat cp_X,
//                   const MapMat cp_Y_X,
//                   MapMat Z,
//                   MapArr2D gam_vb,
//                   const MapArr2D theta_plus_zeta_vb,
//                   const MapArr2D log_Phi_theta_plus_zeta,
//                   const MapArr2D log_1_min_Phi_theta_plus_zeta,
//                   const double log_sig2_inv_vb,
//                   const MapArr1D log_tau_vb,
//                   MapMat m1_beta,
//                   MapMat cp_betaX_X,
//                   MapArr2D mu_beta_vb,
//                   const MapArr1D sig2_beta_vb,
//                   const MapArr1D tau_vb,
//                   const Eigen::VectorXi shuffled_ind,
//                   const Eigen::VectorXi sample_q,
//                   const double c = 1) {
//   
//   // Rcout << "one\n";
//   
//   const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;
//   
//   const bool use_given_log_pnorm = std::abs(c - 1.0) < 1e-8;
//   const double sqrt_c = std::sqrt(c);
//   
//   for (int a = 0; a < sample_q.size(); a++) {
//     int k = sample_q[a];
//     
//     // Rcout << "two" << a << " " << k <<  "\n";
//     // Rcout << sample_q << "\n";
//     
//     for (int b = 0; b < shuffled_ind.size(); b++) {
//       int j = shuffled_ind[b];
//       
//       
//       double m1_beta_jk = m1_beta(j, k);
//       //double cp_betaX_X_kj = cp_betaX_X(k, j) - m1_beta_jk * cp_X(j,j);
//       double cp_betaX_X_jk = cp_betaX_X(j, k) - m1_beta_jk * cp_X(j,j);
//       
//       mu_beta_vb(j, k) = c * sig2_beta_vb[k] * tau_vb[k] * (cp_Y_X(k, j) - cp_betaX_X_jk);
//       
//       gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) - log_Phi_theta_plus_zeta(j, k)
//                                                - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb[k])
//                                                + cst[k])));
//                                                
//        m1_beta(j, k) = gam_vb(j, k) * mu_beta_vb(j, k);
//        
//        cp_betaX_X.col(k) += (m1_beta(j, k) - m1_beta_jk) * cp_X.col(j);
//                                                
//        // ==== Inline update of Z ====
//        double u = sqrt_c * theta_plus_zeta_vb(j, k);
//        
//        double lp = log_Phi_theta_plus_zeta(j, k);
//        double l1p = log_1_min_Phi_theta_plus_zeta(j, k);
//        double log_p = use_given_log_pnorm ? lp : R::pnorm(u, 0.0, 1.0, true, true);
//        double log_1_p = use_given_log_pnorm ? l1p : R::pnorm(u, 0.0, 1.0, false, true);
//        
//        double imr0 = inv_mills_ratio_scalar(0, u, log_1_p, log_p);
//        double imr1 = inv_mills_ratio_scalar(1, u, log_1_p, log_p);
//        
//        Z(j, k) = (gam_vb(j, k) * (imr1 - imr0) + imr0) / sqrt_c + theta_plus_zeta_vb(j, k);
//     }
//   } 
//   
//   // //vectorized update for Z
//   // for (int a = 0; a < sample_q.size(); a++) {
//   //   int k = sample_q[a];
//   //   
//   //   // extract full column k
//   //   Eigen::ArrayXd theta_col = theta_plus_zeta_vb.col(k);
//   //   Eigen::ArrayXd gam_col = gam_vb.col(k);
//   //   Eigen::ArrayXd log_phi_col = log_Phi_theta_plus_zeta.col(k);
//   //   Eigen::ArrayXd log_1_min_phi_col = log_1_min_Phi_theta_plus_zeta.col(k);
//   //   
//   //   Eigen::ArrayXd u = sqrt_c * theta_col;
//   //   
//   //   Eigen::ArrayXd log_p, log_1_p;
//   //   if (use_given_log_pnorm) {
//   //     log_p = log_phi_col;
//   //     log_1_p = log_1_min_phi_col;
//   //   } else {
//   //     log_p = u.unaryExpr([](double x) { return R::pnorm(x, 0.0, 1.0, true, true); });
//   //     log_1_p = u.unaryExpr([](double x) { return R::pnorm(x, 0.0, 1.0, false, true); });
//   //   }
//   //   
//   //   Eigen::ArrayXd imr0 = (-(0.5 * u.square()) - 0.5 * std::log(2.0 * M_PI) - log_1_p).exp();
//   //   imr0 = (-imr0).min(-u);
//   //   
//   //   Eigen::ArrayXd imr1 = (-(0.5 * u.square()) - 0.5 * std::log(2.0 * M_PI) - log_p).exp();
//   //   imr1 = imr1.min(-u);
//   //   
//   //   Eigen::ArrayXd Z_col = (gam_col * (imr1 - imr0) + imr0) / sqrt_c + theta_col;
//   //   
//   //   // store back column k
//   //   Z.col(k) = Z_col;
//   // }
//   
// }
// 
// 
// 
// // // [[Rcpp::export]]
// // List update_log_pnorm_partial(const Eigen::ArrayXXd& theta_plus_zeta_vb, 
// //                               const IntegerVector& sample_q) {
// //   
// //   int n_rows = theta_plus_zeta_vb.rows();
// //   int n_cols = theta_plus_zeta_vb.cols();
// //   
// //   // initialize output matrices (full shape)
// //   Eigen::ArrayXXd log_Phi_theta_plus_zeta = theta_plus_zeta_vb;
// //   Eigen::ArrayXXd log_1_min_Phi_theta_plus_zeta = theta_plus_zeta_vb;
// //   
// //   for (int idx = 0; idx < sample_q.size(); idx++) {
// //     int k = sample_q[idx];
// //     
// //     Eigen::ArrayXd theta_col = theta_plus_zeta_vb.col(k);
// //     
// //     Eigen::ArrayXd log_p = theta_col.unaryExpr([](double x) { 
// //       return R::pnorm(x, 0.0, 1.0, true, true); 
// //     });
// //     
// //     Eigen::ArrayXd log_1_p = theta_col.unaryExpr([](double x) { 
// //       return R::pnorm(x, 0.0, 1.0, false, true); 
// //     });
// //     
// //     log_Phi_theta_plus_zeta.col(k) = log_p;
// //     log_1_min_Phi_theta_plus_zeta.col(k) = log_1_p;
// //   }
// // 
// // }



// for atlasqtl_core function handling missing values in Y
// [[Rcpp::export]]
void coreDualMisLoop(const MapMat cp_X,
                     const List cp_X_rm,
                     const MapMat cp_Y_X,
                     MapArr2D gam_vb,
                     const MapArr2D log_Phi_theta_plus_zeta,
                     const MapArr2D log_1_min_Phi_theta_plus_zeta,
                     const double log_sig2_inv_vb,
                     const MapArr1D log_tau_vb,
                     MapMat m1_beta,
                     MapMat cp_betaX_X, //MapArr2D X_beta_vb,
                     MapArr2D mu_beta_vb,
                     const MapArr2D sig2_beta_vb,
                     const MapArr1D tau_vb,
                     const Eigen::VectorXi shuffled_ind,
                     const Eigen::VectorXi sample_q,
                     const double c = 1) {
  
  const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb)/ 2;
  // int q = cp_Y_X.rows();
  
  for (int a = 0; a < sample_q.size(); a++) {
    int k = sample_q[a];
    MapMat cp_X_rm_k = as<MapMat>(cp_X_rm[k]);
    
    for (int b = 0; b < shuffled_ind.size(); b++) {
      int j = shuffled_ind[b];
      
      
      double m1_beta_jk = m1_beta(j, k);
      //double cp_betaX_X_kj = cp_betaX_X(k, j) - m1_beta(j, k)*(cp_X(j, j) - cp_X_rm_k(j, j));
      double cp_betaX_X_jk = cp_betaX_X(j, k) - m1_beta_jk*(cp_X(j, j) - cp_X_rm_k(j, j));
      
      // cp_betaX_X.row(k) += -m1_beta(j, k) * (cp_X.row(j) - cp_X_rm_k.row(j));
      
      mu_beta_vb(j, k) = c * sig2_beta_vb(j,k) * tau_vb[k] * (cp_Y_X(k, j) - cp_betaX_X_jk);
      
      gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) -
        log_Phi_theta_plus_zeta(j, k) - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb(j, k)) -
        log(sig2_beta_vb(j, k)) / 2 + cst[k])));
      
      m1_beta(j, k) = gam_vb(j, k) * mu_beta_vb(j, k);
      cp_betaX_X.col(k) += (m1_beta(j, k) - m1_beta_jk) * (cp_X.col(j) - cp_X_rm_k.col(j)); 
      
      
    } 
  }
  
}
