// C++ standard library headers
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>

// All R headers must be included within an extern "C" block
// when compiling with a C++ compiler.
extern "C" {
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
}

// Helper function for column-major indexing (as used by R and FORTRAN)
inline int idx(int row, int col, int num_rows) {
  return row + col * num_rows;
}

// Fisher-Yates shuffle using R's random number generator for reproducibility
void shuffle_indices(std::vector<int>& vec) {
  GetRNGstate();
  const int n = static_cast<int>(vec.size());
  for (int i = n - 1; i > 0; i--) {
    int j = static_cast<int>(unif_rand() * (i + 1.0));
    std::swap(vec[i], vec[j]);
  }
  PutRNGstate();
}

extern "C" {
  
  SEXP PEGS_lapack(SEXP Y_R, SEXP X_list_R, SEXP maxit_R, SEXP logtol_R, 
                   SEXP NonNegativeCorr_R, SEXP covbend_R, SEXP covMinEv_R, SEXP XFA_R) {
    
    // --- 1. Get Input Dimensions and Parameters ---
    const int n0 = nrows(Y_R);
    const int k  = ncols(Y_R);
    const double *Y_ptr = REAL(Y_R);
    
    const int    maxit            = asInteger(maxit_R);
    const double logtol           = asReal(logtol_R);
    const bool   NonNegativeCorr  = asLogical(NonNegativeCorr_R);
    const double covbend          = asReal(covbend_R);
    const double covMinEv         = asReal(covMinEv_R);
    int          XFA              = asInteger(XFA_R);
    if (XFA < 0) XFA = k; // If XFA is negative, use all components (no reduction)
    
    const int one = 1; // Used for BLAS calls
    
    if (!isVector(X_list_R)) {
      error("X_list_R must be a list of numeric matrices.");
    }
    const int n_effects = LENGTH(X_list_R);
    std::vector<const double*> X_ptr_list(n_effects);
    std::vector<int> p_vec(n_effects);
    std::vector<int> p_cumsum(n_effects + 1, 0);
    int p_total = 0;
    
    for(int i = 0; i < n_effects; ++i) {
      SEXP X_i_R = VECTOR_ELT(X_list_R, i);
      if (!isMatrix(X_i_R) || !isNumeric(X_i_R)) {
        error("Element %d of X_list_R is not a numeric matrix.", i + 1);
      }
      if (nrows(X_i_R) != n0) {
        error("Matrix %d in X list has %d rows, but Y has %d rows.", i + 1, nrows(X_i_R), n0);
      }
      X_ptr_list[i] = REAL(X_i_R);
      p_vec[i] = ncols(X_i_R);
      p_total += p_vec[i];
      p_cumsum[i+1] = p_total;
    }
    
    // --- 2. Data Preparation and Memory Allocation ---
    std::vector<double> Y(n0 * k);
    std::vector<double> Z(n0 * k, 0.0);
    
    for (int i = 0; i < n0 * k; ++i) {
      if (ISNAN(Y_ptr[i])) {
        Y[i] = 0.0;
      } else {
        Y[i] = Y_ptr[i];
        Z[i] = 1.0;
      }
    }
    
    std::vector<double> n(k, 0.0);
    for (int j = 0; j < k; ++j) {
      for (int i = 0; i < n0; ++i) n[j] += Z[idx(i, j, n0)];
    }
    
    std::vector<double> iN_orig(k);
    for (int i = 0; i < k; ++i) iN_orig[i] = (n[i] > 0) ? 1.0 / n[i] : 0.0;
    
    std::vector<double> mu(k, 0.0);
    for (int j = 0; j < k; ++j) {
      for (int i = 0; i < n0; ++i) mu[j] += Y[idx(i, j, n0)];
      mu[j] *= iN_orig[j];
    }
    
    std::vector<double> y(n0 * k);
    for (int j = 0; j < k; ++j) {
      for (int i = 0; i < n0; ++i) {
        y[idx(i, j, n0)] = (Y[idx(i, j, n0)] - mu[j]) * Z[idx(i, j, n0)];
      }
    }
    
    // --- 3. Pre-computation for Loop ---
    std::vector<std::vector<double>> XX_list(n_effects);
    std::vector<std::vector<double>> XSX_list(n_effects);
    std::vector<double> MSx_list(n_effects * k, 0.0);
    std::vector<double> TrXSX_list(n_effects * k, 0.0);
    std::vector<std::vector<double>> tilde_list(n_effects);
    
    for (int i_effect = 0; i_effect < n_effects; ++i_effect) {
      const int p_i = p_vec[i_effect];
      const double* X_ptr_i = X_ptr_list[i_effect];
      
      XX_list[i_effect].resize(p_i * k);
      XSX_list[i_effect].resize(p_i * k);
      tilde_list[i_effect].resize(p_i * k);
      
      for (int j = 0; j < k; ++j) {
        for (int i = 0; i < p_i; ++i) {
          double sum_sq = 0.0, sum_xz_in = 0.0;
          for (int r = 0; r < n0; ++r) {
            const double x_val = X_ptr_i[idx(r, i, n0)];
            const double z_val = Z[idx(r, j, n0)];
            sum_sq += x_val * x_val * z_val;
            sum_xz_in += x_val * z_val;
          }
          XX_list[i_effect][idx(i, j, p_i)] = sum_sq;
          sum_xz_in *= iN_orig[j];
          XSX_list[i_effect][idx(i, j, p_i)] = XX_list[i_effect][idx(i, j, p_i)] * iN_orig[j] - (sum_xz_in * sum_xz_in);
        }
      }
      
      for (int j = 0; j < k; ++j) {
        for (int i = 0; i < p_i; ++i) {
          MSx_list[idx(i_effect, j, n_effects)] += XSX_list[i_effect][idx(i, j, p_i)];
        }
      }
      
      for (int j = 0; j < k; ++j) {
        TrXSX_list[idx(i_effect, j, n_effects)] = n[j] * MSx_list[idx(i_effect, j, n_effects)];
      }
      
      char transa = 'T', transb = 'N';
      double alpha = 1.0, beta = 0.0;
      F77_CALL(dgemm)(&transa, &transb, &p_i, &k, &n0, &alpha, X_ptr_i, &n0, y.data(), &n0, &beta, tilde_list[i_effect].data(), &p_i FCONE FCONE);
    }
    
    std::vector<double> iN(k);
    for (int i = 0; i < k; ++i) iN[i] = (n[i] > 1) ? 1.0 / (n[i] - 1.0) : 0.0;
    
    std::vector<double> vy(k);
    for (int j = 0; j < k; ++j) {
      double y_norm_sq = F77_CALL(ddot)(&n0, &y[j * n0], &one, &y[j * n0], &one);
      vy[j] = y_norm_sq * iN[j];
    }
    
    std::vector<double> ve = vy;
    for (double &val : ve) val *= 0.5;
    
    std::vector<double> iVe(k);
    for (int i = 0; i < k; ++i) iVe[i] = (ve[i] > 0) ? 1.0 / ve[i] : 0.0;
    
    std::vector<std::vector<double>> vb_list(n_effects, std::vector<double>(k * k, 0.0));
    std::vector<std::vector<double>> iG_list(n_effects, std::vector<double>(k * k, 0.0));
    
    for (int i_effect = 0; i_effect < n_effects; ++i_effect) {
      for (int i = 0; i < k; ++i) {
        double msx_val = MSx_list[idx(i_effect, i, n_effects)];
        double& vb_diag = vb_list[i_effect][idx(i, i, k)];
        vb_diag = (msx_val > 0) ? ve[i] / msx_val : 0.0;
        if (vb_diag > 0) {
          iG_list[i_effect][idx(i, i, k)] = 1.0 / vb_diag;
        }
      }
    }
    
    // --- 4. Initialize Iteration Variables ---
    std::vector<std::vector<double>> b_list(n_effects);
    for(int i=0; i<n_effects; ++i) b_list[i].assign(p_vec[i] * k, 0.0);
    
    std::vector<double> e = y;
    std::vector<int> RGSvec(p_total);
    std::iota(RGSvec.begin(), RGSvec.end(), 0);
    
    std::vector<double> inflate_list(n_effects, 0.0);
    double cnv = 10.0;
    int numit = 0;
    
    // --- 5. Main Iteration Loop ---
    std::vector<double> LHS(k * k);
    
    while (numit < maxit) {
      auto beta0_list = b_list;
      shuffle_indices(RGSvec);
      
      for (int j_idx = 0; j_idx < p_total; ++j_idx) {
        const int J_global = RGSvec[j_idx];
        
        auto it = std::upper_bound(p_cumsum.begin(), p_cumsum.end(), J_global);
        const int effect_idx = std::distance(p_cumsum.begin(), it) - 1;
        const int local_J = J_global - p_cumsum[effect_idx];
        
        const int p_i = p_vec[effect_idx];
        const double* X_ptr_i = X_ptr_list[effect_idx];
        
        std::vector<double> b0(k);
        for (int i = 0; i < k; ++i) b0[i] = b_list[effect_idx][idx(local_J, i, p_i)];
        
        std::copy(iG_list[effect_idx].begin(), iG_list[effect_idx].end(), LHS.begin());
        for (int i = 0; i < k; ++i) {
          LHS[idx(i, i, k)] += XX_list[effect_idx][idx(local_J, i, p_i)] * iVe[i];
        }
        
        std::vector<double> RHS(k);
        for (int i = 0; i < k; ++i) {
          RHS[i] = F77_CALL(ddot)(&n0, &X_ptr_i[local_J * n0], &one, &e[i * n0], &one);
          RHS[i] += XX_list[effect_idx][idx(local_J, i, p_i)] * b0[i];
          RHS[i] *= iVe[i];
        }
        
        char uplo = 'U'; int info;
        F77_CALL(dpotrf)(&uplo, &k, LHS.data(), &k, &info FCONE);
        if (info != 0) {
          warning("Cholesky factorization failed for effect %d, marker %d. Skipping update.", effect_idx + 1, local_J + 1);
          continue;
        }
        int nrhs = 1;
        F77_CALL(dpotrs)(&uplo, &k, &nrhs, LHS.data(), &k, RHS.data(), &k, &info FCONE);
        const std::vector<double> b1 = RHS;
        
        for (int i = 0; i < k; ++i) {
          const double delta_b = b1[i] - b0[i];
          b_list[effect_idx][idx(local_J, i, p_i)] = b1[i];
          if (fabs(delta_b) > 1e-12) { 
            const double neg_delta_b = -delta_b;
            F77_CALL(daxpy)(&n0, &neg_delta_b, &X_ptr_i[local_J * n0], &one, &e[i * n0], &one);
          }
        }
      }
      
      for (int j = 0; j < k; ++j) {
        const double ve_sum = F77_CALL(ddot)(&n0, &e[j * n0], &one, &y[j * n0], &one);
        ve[j]  = ve_sum * iN[j];
      }
      for(int i=0; i<k; ++i) iVe[i] = (ve[i] > 0) ? 1.0 / ve[i] : 0.0;
      
      for (int i_effect = 0; i_effect < n_effects; ++i_effect) {
        const int p_i = p_vec[i_effect];
        std::vector<double> TildeHat(k * k);
        char transa = 'T', transb = 'N';
        double alpha = 1.0, beta = 0.0;
        F77_CALL(dgemm)(&transa, &transb, &k, &k, &p_i, &alpha,
                 b_list[i_effect].data(), &p_i,
                 tilde_list[i_effect].data(), &p_i,
                 &beta, TildeHat.data(), &k FCONE FCONE);
        
        for (int c = 0; c < k; ++c) {
          for (int r = c; r < k; ++r) {
            const double tr_r = TrXSX_list[idx(i_effect, r, n_effects)];
            const double tr_c = TrXSX_list[idx(i_effect, c, n_effects)];
            if (r == c) {
              vb_list[i_effect][idx(r, c, k)] = (tr_r > 0) ? TildeHat[idx(r, c, k)] / tr_r : 0.0;
            } else {
              const double denom = tr_r + tr_c;
              const double val = (denom > 0) ? (TildeHat[idx(r, c, k)] + TildeHat[idx(c, r, k)]) / denom : 0.0;
              vb_list[i_effect][idx(r, c, k)] = vb_list[i_effect][idx(c, r, k)] = val;
            }
          }
        }
        
        // --- START: New XFA and Bending Logic ---
        
        // XFA: Factor-analytic approximation
        if (XFA == 0) {
          // Make traits independent, keeping only diagonal variance
          std::vector<double> sd_diag(k);
          for(int i=0; i<k; ++i) sd_diag[i] = vb_list[i_effect][idx(i,i,k)];
          std::fill(vb_list[i_effect].begin(), vb_list[i_effect].end(), 0.0);
          for(int i=0; i<k; ++i) vb_list[i_effect][idx(i,i,k)] = sd_diag[i];
        } else if (XFA > 0 && XFA < k) {
          // Reconstruct vb from top XFA principal components
          std::vector<double> sd(k), inv_sd(k);
          for (int i = 0; i < k; ++i) {
            sd[i] = sqrt(std::max(vb_list[i_effect][idx(i, i, k)], 1e-24));
            inv_sd[i] = (sd[i] > 1e-12) ? 1.0 / sd[i] : 0.0;
          }
          
          // Convert vb to correlation matrix GC
          std::vector<double> GC = vb_list[i_effect];
          for (int c = 0; c < k; ++c) {
            for (int r = 0; r < k; ++r) {
              GC[idx(r, c, k)] *= inv_sd[r] * inv_sd[c];
            }
          }
          
          // Eigen decomposition of GC
          char jobz_eig = 'V', uplo_eig = 'U'; int info_eig;
          int lwork = std::max(1, 3 * k - 1);
          std::vector<double> work(lwork), eigvals(k);
          F77_CALL(dsyev)(&jobz_eig, &uplo_eig, &k, GC.data(), &k, eigvals.data(), work.data(), &lwork, &info_eig FCONE FCONE);
          
          // Reconstruct GC = V_reduced * D_reduced_diag * V_reduced'
          std::vector<double> V_reduced(k * XFA);
          std::vector<double> Temp(k * XFA);
          double alpha_gemm = 1.0, beta_gemm = 0.0;
          
          // V_reduced are the last XFA columns of GC (eigenvectors)
          // Temp = V_reduced * D_reduced_diag
          for (int j = 0; j < XFA; ++j) {
            int eig_idx = k - XFA + j;
            for (int i = 0; i < k; ++i) {
              Temp[idx(i, j, k)] = GC[idx(i, eig_idx, k)] * eigvals[eig_idx];
            }
          }
          // GC = Temp * V_reduced'
          char transa_gemm = 'N', transb_gemm = 'T';
          F77_CALL(dgemm)(&transa_gemm, &transb_gemm, &k, &k, &XFA, &alpha_gemm, Temp.data(), &k, &GC[idx(0, k - XFA, k)], &k, &beta_gemm, vb_list[i_effect].data(), &k FCONE FCONE);
          
          // Rescale back to covariance matrix
          for (int c = 0; c < k; ++c) {
            for (int r = c; r < k; ++r) {
              double val = vb_list[i_effect][idx(r, c, k)] * sd[r] * sd[c];
              vb_list[i_effect][idx(r, c, k)] = val;
              vb_list[i_effect][idx(c, r, k)] = val;
            }
            vb_list[i_effect][idx(c,c,k)] = sd[c]*sd[c]; // Enforce diagonal
          }
        }
        
        // Bending: Ensure positive semi-definiteness
        if (NonNegativeCorr) {
          for (double &val : vb_list[i_effect]) if (val < 0.0) val = 0.0;
        }
        
        std::vector<double> vb_copy_eig = vb_list[i_effect];
        std::vector<double> eigvals(k);
        char jobz = 'N', uplo_eig = 'U'; int info_eig;
        int lwork = std::max(1, 3 * k - 1); std::vector<double> work(lwork);
        F77_CALL(dsyev)(&jobz, &uplo_eig, &k, vb_copy_eig.data(), &k, eigvals.data(), work.data(), &lwork, &info_eig FCONE FCONE);
        
        double MinDVb = eigvals[0];
        for (int i = 1; i < k; ++i) if (eigvals[i] < MinDVb) MinDVb = eigvals[i];
        
        if (MinDVb < covMinEv) {
          const double new_inflate = fabs(MinDVb * covbend);
          if (new_inflate > inflate_list[i_effect]) inflate_list[i_effect] = new_inflate;
        }
        
        if(k >= 5 || MinDVb < covMinEv){
          for (int i = 0; i < k; ++i) vb_list[i_effect][idx(i, i, k)] += inflate_list[i_effect];
        }
        // --- END: New XFA and Bending Logic ---
        
        // Invert vb to get iG
        std::vector<double> vb_copy_inv = vb_list[i_effect];
        char uplo_inv = 'U'; int info_inv;
        F77_CALL(dpotrf)(&uplo_inv, &k, vb_copy_inv.data(), &k, &info_inv FCONE);
        if (info_inv == 0) {
          F77_CALL(dpotri)(&uplo_inv, &k, vb_copy_inv.data(), &k, &info_inv FCONE);
          if (info_inv == 0) {
            for (int c = 0; c < k; ++c)
              for (int r = c + 1; r < k; ++r)
                vb_copy_inv[idx(r, c, k)] = vb_copy_inv[idx(c, r, k)];
            iG_list[i_effect] = vb_copy_inv;
          } else {
            warning("dpotri failed for effect %d; iG not updated.", i_effect + 1);
          }
        } else {
          warning("dpotrf failed for bent vb in effect %d; iG not updated.", i_effect + 1);
        }
      } 
      
      std::vector<double> b0_e(k, 0.0);
      for (int j = 0; j < k; ++j) {
        for (int i = 0; i < n0; ++i) b0_e[j] += e[idx(i, j, n0)];
        b0_e[j] *= iN_orig[j]; 
        mu[j] += b0_e[j];
      }
      for (int j = 0; j < k; ++j) {
        for (int i = 0; i < n0; ++i) {
          e[idx(i, j, n0)] = (e[idx(i, j, n0)] - b0_e[j]) * Z[idx(i, j, n0)];
        }
      }
      
      double diff_sum_sq = 0.0;
      for (int i_effect = 0; i_effect < n_effects; ++i_effect) {
        for (size_t i = 0; i < b_list[i_effect].size(); ++i) {
          const double diff = beta0_list[i_effect][i] - b_list[i_effect][i];
          diff_sum_sq += diff * diff;
        }
      }
      cnv = (diff_sum_sq > 0) ? log10(diff_sum_sq) : -INFINITY;
      ++numit;
      if (numit % 100 == 0) { Rprintf("Iter: %d || Conv: %f\n", numit, cnv); }
      if (cnv < logtol) break;
    }
    
    // --- 6. Final Calculations and Output ---
    std::vector<double> h2(k);
    for (int i = 0; i < k; ++i) h2[i] = (vy[i] > 0) ? 1.0 - ve[i] / vy[i] : 0.0;
    
    std::vector<double> hat(n0 * k, 0.0);
    for (int i_effect = 0; i_effect < n_effects; ++i_effect) {
      char transa_final = 'N', transb_final = 'N';
      double beta = (i_effect == 0) ? 0.0 : 1.0;
      double alpha = 1.0;
      F77_CALL(dgemm)(&transa_final, &transb_final, &n0, &k, &p_vec[i_effect],
               &alpha, X_ptr_list[i_effect], &n0,
               b_list[i_effect].data(), &p_vec[i_effect],
                                              &beta, hat.data(), &n0 FCONE FCONE);
    }
    for (int j = 0; j < k; ++j) 
      for (int i = 0; i < n0; ++i) 
        hat[idx(i, j, n0)] += mu[j];
    
    std::vector<std::vector<double>> GC_list(n_effects, std::vector<double>(k*k, 0.0));
    for (int i_effect = 0; i_effect < n_effects; ++i_effect) {
      for (int c = 0; c < k; ++c) {
        for (int r = c; r < k; ++r) {
          const double sd_r = sqrt(std::max(0.0, vb_list[i_effect][idx(r, r, k)]));
          const double sd_c = sqrt(std::max(0.0, vb_list[i_effect][idx(c, c, k)]));
          if (sd_r > 1e-12 && sd_c > 1e-12) {
            const double val = vb_list[i_effect][idx(r, c, k)] / (sd_r * sd_c);
            GC_list[i_effect][idx(r, c, k)] = val;
            GC_list[i_effect][idx(c, r, k)] = val;
          }
        }
      }
    }
    
    // --- 7. Create R List for Output ---
    const char *names[] = {"mu", "b_list", "hat", "h2", "GC_list", "bend_list", "numit", "cnv", ""};
    SEXP res = PROTECT(mkNamed(VECSXP, names));
    
    SEXP mu_R = PROTECT(allocVector(REALSXP, k));
    if (k > 0) std::copy(mu.begin(), mu.end(), REAL(mu_R));
    SET_VECTOR_ELT(res, 0, mu_R);
    
    SEXP b_R_list = PROTECT(allocVector(VECSXP, n_effects));
    for (int i = 0; i < n_effects; ++i) {
      SEXP b_i_R = PROTECT(allocMatrix(REALSXP, p_vec[i], k));
      if (p_vec[i] * k > 0) std::copy(b_list[i].begin(), b_list[i].end(), REAL(b_i_R));
      SET_VECTOR_ELT(b_R_list, i, b_i_R);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(res, 1, b_R_list);
    
    SEXP hat_R = PROTECT(allocMatrix(REALSXP, n0, k));
    if (n0 * k > 0) std::copy(hat.begin(), hat.end(), REAL(hat_R));
    SET_VECTOR_ELT(res, 2, hat_R);
    
    SEXP h2_R = PROTECT(allocVector(REALSXP, k));
    if (k > 0) std::copy(h2.begin(), h2.end(), REAL(h2_R));
    SET_VECTOR_ELT(res, 3, h2_R);
    
    SEXP GC_R_list = PROTECT(allocVector(VECSXP, n_effects));
    for (int i = 0; i < n_effects; ++i) {
      SEXP GC_i_R = PROTECT(allocMatrix(REALSXP, k, k));
      if (k*k > 0) std::copy(GC_list[i].begin(), GC_list[i].end(), REAL(GC_i_R));
      SET_VECTOR_ELT(GC_R_list, i, GC_i_R);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(res, 4, GC_R_list);
    
    SEXP bend_R = PROTECT(allocVector(REALSXP, n_effects));
    if (n_effects > 0) std::copy(inflate_list.begin(), inflate_list.end(), REAL(bend_R));
    SET_VECTOR_ELT(res, 5, bend_R);
    
    SET_VECTOR_ELT(res, 6, ScalarInteger(numit));
    SET_VECTOR_ELT(res, 7, ScalarReal(cnv));
    
    UNPROTECT(7); 
    return res;
  }
  
} // extern "C"