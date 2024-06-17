functions{
  
  // Gaussian Copula Log Probability Density
  
  // Gaussian Copula Log Probability Density
  
  real multi_normal_cholesky_copula_lpdf(matrix U, matrix L) {
    int N = rows(U);
    int J = cols(U);
    matrix[J, J] Gammainv = chol2inv(L);
    return -N * sum(log(diagonal(L))) - 0.5 * sum(add_diag(Gammainv, -1.0) .* crossprod(U));
  }
  
  // Prepare data for LPDF
  
  real centered_gaussian_copula_cholesky_(array[,] matrix marginals, matrix L) {
    // Extract dimensions of final outcome matrix
    int N = rows(marginals[1][1]);
    int J = rows(L);
    matrix[N, J] U;
  
    // Iterate through marginal arrays, concatenating the outcome matrices by column
    // and aggregating the log-likelihoods (from continuous marginals) and jacobian
    // adjustments (from discrete marginals)
    real adj = 0;
    int pos = 1;
    for (m in 1 : size(marginals)) {
      int curr_cols = cols(marginals[m][1]);
    
      U[ : , pos : (pos + curr_cols - 1)] = marginals[m][1];
    
      adj += sum(marginals[m][2]);
      pos += curr_cols;
    }
  
    // Return the sum of the log-probability for copula outcomes and jacobian adjustments
    return multi_normal_cholesky_copula_lpdf(U | L) + adj;
  }
  

  // Continuous Marginal Distribution (Normal)
  
  array[] matrix normal_marginal(array[] real y, array[] real mu, array[] real sigma, int N) {
    array[2] matrix[N, 1] rtn; // empty 2D array to return
    // Initialise the jacobian adjustments to zero, as vectorised lpdf will be used
    rtn[2] = rep_matrix(0, N, 1);

    for (n in 1 : N) {
      rtn[1][n, 1] = (y[n] - mu[n]) / sigma[n]; // center RV
      rtn[2][n, 1] = normal_lpdf(y[n] | mu[n], sigma[n]); // "jacobian"
    }
  return rtn;
  } 

  array[] matrix probit_marginal(array[] int y, array[] real mu_glm, array[] real u_raw, vector cutpoints) {
    int N = num_elements(mu_glm); // number of observations
    array[2] matrix[N, 1] rtn; // empty 2D array to return
    
    for(n in 1:N){
      int C = num_elements(cutpoints) + 1; // total number of ord categories
      if(y[n] == 99){ // missing data
        rtn[1][n,1] = inv_Phi(u_raw[n]); // missing RV
        rtn[2][n,1] = 0;
        } else if(y[n] == 1){ // lowest bound
        real bound = Phi((cutpoints[1] - mu_glm[n])); // data augmentation
        rtn[1][n,1] = inv_Phi((bound*u_raw[n])); // latent RV
        rtn[2][n,1] = log(bound); // jacobian
      } else if (y[n] == C){ // highest bound
        real bound = Phi((cutpoints[C - 1] - mu_glm[n])); // data augmentation
        rtn[1][n,1] = inv_Phi(bound + (1-bound)*u_raw[n]); // latent RV
        rtn[2][n,1] = log1m(bound); // jacobian
      } else { // in between 
        real ub = Phi((cutpoints[y[n]] - mu_glm[n])); // data augmentation
        real lb = Phi((cutpoints[y[n] - 1] - mu_glm[n])); // data augmentation
        rtn[1][n,1] = inv_Phi((lb + (ub-lb)*u_raw[n])); // latent RV
        rtn[2][n,1] = log(ub-lb); // jacobian
      }
    }
    return rtn;
  }
}
data{
  
  // Global Variables
  
  int N; // total number of observations (rows)
  int M; // total number of response variables (columns)
  array[N] real x; // predictor variables (chronological age)
  
  // Continuous Variables
  
  ///FDL
  int<lower=0> FDL_complete; // # of complete observations
  int<lower=0> FDL_missing; // # of missing observations
  array[FDL_complete] real y_FDL; // observed response variable
  array[FDL_complete] int FDL_index_present; // index of present variable
  array[FDL_missing] int FDL_index_missing; // index of missing variable
  
  ///FMSB
  int<lower=0> FMSB_complete; // # of complete observations
  int<lower=0> FMSB_missing; // # of missing observations
  array[FMSB_complete] real y_FMSB; // observed response variable
  array[FMSB_complete] int FMSB_index_present; // index of present variable
  array[FMSB_missing] int FMSB_index_missing; // index of missing variable

  ///FDB
  int<lower=0> FDB_complete; // # of complete observations
  int<lower=0> FDB_missing; // # of missing observations
  array[FDB_complete] real y_FDB; // observed response variable
  array[FDB_complete] int FDB_index_present; // index of present variable
  array[FDB_missing] int FDB_index_missing; // index of missing variable

  ///TDL
  int<lower=0> TDL_complete; // # of complete observations
  int<lower=0> TDL_missing; // # of missing observations
  array[TDL_complete] real y_TDL; // observed response variable
  array[TDL_complete] int TDL_index_present; // index of present variable
  array[TDL_missing] int TDL_index_missing; // index of missing variable

  ///TPB
  int<lower=0> TPB_complete; // # of complete observations
  int<lower=0> TPB_missing; // # of missing observations
  array[TPB_complete] real y_TPB; // observed response variable
  array[TPB_complete] int TPB_index_present; // index of present variable
  array[TPB_missing] int TPB_index_missing; // index of missing variable

  ///TMSB
  int<lower=0> TMSB_complete; // # of complete observations
  int<lower=0> TMSB_missing; // # of missing observations
  array[TMSB_complete] real y_TMSB; // observed response variable
  array[TMSB_complete] int TMSB_index_present; // index of present variable
  array[TMSB_missing] int TMSB_index_missing; // index of missing variable

  ///TDB
  int<lower=0> TDB_complete; // # of complete observations
  int<lower=0> TDB_missing; // # of missing observations
  array[TDB_complete] real y_TDB; // observed response variable
  array[TDB_complete] int TDB_index_present; // index of present variable
  array[TDB_missing] int TDB_index_missing; // index of missing variable

  ///FBDL
  int<lower=0> FBDL_complete; // # of complete observations
  int<lower=0> FBDL_missing; // # of missing observations
  array[FBDL_complete] real y_FBDL; // observed response variable
  array[FBDL_complete] int FBDL_index_present; // index of present variable
  array[FBDL_missing] int FBDL_index_missing; // index of missing variable

  ///HDL
  int<lower=0> HDL_complete; // # of complete observations
  int<lower=0> HDL_missing; // # of missing observations
  array[HDL_complete] real y_HDL; // observed response variable
  array[HDL_complete] int HDL_index_present; // index of present variable
  array[HDL_missing] int HDL_index_missing; // index of missing variable

  ///HPB
  int<lower=0> HPB_complete; // # of complete observations
  int<lower=0> HPB_missing; // # of missing observations
  array[HPB_complete] real y_HPB; // observed response variable
  array[HPB_complete] int HPB_index_present; // index of present variable
  array[HPB_missing] int HPB_index_missing; // index of missing variable

  ///HMSB
  int<lower=0> HMSB_complete; // # of complete observations
  int<lower=0> HMSB_missing; // # of missing observations
  array[HMSB_complete] real y_HMSB; // observed response variable
  array[HMSB_complete] int HMSB_index_present; // index of present variable
  array[HMSB_missing] int HMSB_index_missing; // index of missing variable

  ///HDB
  int<lower=0> HDB_complete; // # of complete observations
  int<lower=0> HDB_missing; // # of missing observations
  array[HDB_complete] real y_HDB; // observed response variable
  array[HDB_complete] int HDB_index_present; // index of present variable
  array[HDB_missing] int HDB_index_missing; // index of missing variable

  ///RDL
  int<lower=0> RDL_complete; // # of complete observations
  int<lower=0> RDL_missing; // # of missing observations
  array[RDL_complete] real y_RDL; // observed response variable
  array[RDL_complete] int RDL_index_present; // index of present variable
  array[RDL_missing] int RDL_index_missing; // index of missing variable

  ///RPB
  int<lower=0> RPB_complete; // # of complete observations
  int<lower=0> RPB_missing; // # of missing observations
  array[RPB_complete] real y_RPB; // observed response variable
  array[RPB_complete] int RPB_index_present; // index of present variable
  array[RPB_missing] int RPB_index_missing; // index of missing variable
  
  ///RMSB
  int<lower=0> RMSB_complete; // # of complete observations
  int<lower=0> RMSB_missing; // # of missing observations
  array[RMSB_complete] real y_RMSB; // observed response variable
  array[RMSB_complete] int RMSB_index_present; // index of present variable
  array[RMSB_missing] int RMSB_index_missing; // index of missing variable
  
  ///RDB
  int<lower=0> RDB_complete; // # of complete observations
  int<lower=0> RDB_missing; // # of missing observations
  array[RDB_complete] real y_RDB; // observed response variable
  array[RDB_complete] int RDB_index_present; // index of present variable
  array[RDB_missing] int RDB_index_missing; // index of missing variable
  
  ///UDL
  int<lower=0> UDL_complete; // # of complete observations
  int<lower=0> UDL_missing; // # of missing observations
  array[UDL_complete] real y_UDL; // observed response variable
  array[UDL_complete] int UDL_index_present; // index of present variable
  array[UDL_missing] int UDL_index_missing; // index of missing variable
  
  ///UMSB
  int<lower=0> UMSB_complete; // # of complete observations
  int<lower=0> UMSB_missing; // # of missing observations
  array[UMSB_complete] real y_UMSB; // observed response variable
  array[UMSB_complete] int UMSB_index_present; // index of present variable
  array[UMSB_missing] int UMSB_index_missing; // index of missing variable

  // Ordinal Variables
  
  /// Dental Variables
  int C_dental; // # of dental categories
  
  ///max_m1
  array[N] int y_max_m1; //score

  ///max_m2
  array[N] int y_max_m2; //score

  ///max_m3
  array[N] int y_max_m3; //score

  ///max_pm1
  array[N] int y_max_pm1; //score

  ///max_pm2
  array[N] int y_max_pm2; //score

  ///max_c
  array[N] int y_max_c; //score


  ///max_i1
  array[N] int y_max_i1; //score

  ///max_i2
  array[N] int y_max_i2; //score

  ///man_m1
  array[N] int y_man_m1; //score

  ///man_m2
  array[N] int y_man_m2; //score

  ///man_m3
  array[N] int y_man_m3; //score

  ///man_pm1
  array[N] int y_man_pm1; //score

  ///man_pm2
  array[N] int y_man_pm2; //score

  ///max_c
  array[N] int y_man_c; //score

  ///max_i1
  array[N] int y_man_i1; //score
  
  ///max_i2
  array[N] int y_man_i2; //score
  
  /// Long Bone Epiphyses (and Calcaneus)
  int C_lb_ef;

  ///FH_EF
  array[N] int y_fh_ef; //score

  ///FGT_EF
  array[N] int y_fgt_ef; //score

  ///FLT_EF
  array[N] int y_flt_ef; //score 

  ///FDE_EF
  array[N] int y_fde_ef; //score 

  ///TPE_EF
  array[N] int y_tpe_ef; //score

  ///TDE_EF
  array[N] int y_tde_ef; //score

  ///FBPE_EF
  array[N] int y_fbpe_ef; //score
  
  ///FBDE_EF
  array[N] int y_fbde_ef; //score
  
  ///HPE_EF
  array[N] int y_hpe_ef; //score
  
  ///HDE_EF
  array[N] int y_hde_ef; //score
  
  ///HME_EF
  array[N] int y_hme_ef; //score
  
  ///RPE_EF
  array[N] int y_rpe_ef; //score
  
  ///RDE_EF
  array[N] int y_rde_ef; //score
  
  ///UPE_EF
  array[N] int y_upe_ef; //score
  
  ///UDE_EF
  array[N] int y_ude_ef; //score
  
  ///CT_EF
  array[N] int y_ct_ef; //score

  /// Pelvis Fusion
  int C_pelvis_ef;
  
  ///ISPR_EF
  array[N] int y_ispr_ef; // score
  
  ///ILIS_EF
  array[N] int y_ilis_ef; // score

  /// Carpal and Tarsal Ossification
  int C_carpal_oss;
  int C_tarsal_oss;
  
  ///CC_Oss
  array[N] int y_cc_oss; // score
  
  ///TC_Oss
  array[N] int y_tc_oss; // score
  
}
transformed data{
  
  int t_dental = C_dental - 1; // dental thresholds size
  int t_lb_ef = C_lb_ef - 1; // lb ef thresholds size
  int t_pelvis_ef = C_pelvis_ef - 1; // pelvis ef thresholds size
  int t_carpal = C_carpal_oss - 1; // carpal thresholds size
  int t_tarsal = C_tarsal_oss - 1; // carpal thresholds size
  
}
parameters{
  
  // Global Parameter
  cholesky_factor_corr[M] L; 
  
  // Continuous Parameters
  
  ///FDL
  real<lower=0> FDL_constant; // mean function
  real<lower=0> FDL_exponent; // mean function
  real<lower=0> FDL_offset; // mean function 
  real<lower=0> FDL_noise_intercept; // sd function
  real<lower=0> FDL_noise_slope; // sd function
  array[FDL_missing] real<lower=0> FDL_y_missing; // missing data parameter
  
  ///FMSB
  real<lower=0> FMSB_constant; // mean function
  real<lower=0> FMSB_exponent; // mean function
  real<lower=0> FMSB_offset; // mean function 
  real<lower=0> FMSB_noise_intercept; // sd function
  real<lower=0> FMSB_noise_slope; // sd function
  array[FMSB_missing] real<lower=0> FMSB_y_missing; // missing data parameter

  ///FDB
  real<lower=0> FDB_constant; // mean function
  real<lower=0> FDB_exponent; // mean function
  real<lower=0> FDB_offset; // mean function 
  real<lower=0> FDB_noise_intercept; // sd function
  real<lower=0> FDB_noise_slope; // sd function
  array[FDB_missing] real<lower=0> FDB_y_missing; // missing data parameter

  ///TDL
  real<lower=0> TDL_constant; // mean function
  real<lower=0> TDL_exponent; // mean function
  real<lower=0> TDL_offset; // mean function 
  real<lower=0> TDL_noise_intercept; // sd function
  real<lower=0> TDL_noise_slope; // sd function
  array[TDL_missing] real<lower=0> TDL_y_missing; // missing data parameter

  ///TPB
  real<lower=0> TPB_constant; // mean function
  real<lower=0> TPB_exponent; // mean function
  real<lower=0> TPB_offset; // mean function 
  real<lower=0> TPB_noise_intercept; // sd function
  real<lower=0> TPB_noise_slope; // sd function
  array[TPB_missing] real<lower=0> TPB_y_missing; // missing data parameter

  ///TMSB
  real<lower=0> TMSB_constant; // mean function
  real<lower=0> TMSB_exponent; // mean function
  real<lower=0> TMSB_offset; // mean function
  real<lower=0> TMSB_sigma; // sd function 
  array[TMSB_missing] real<lower=0> TMSB_y_missing; // missing data parameter

  ///TDB
  real<lower=0> TDB_constant; // mean function
  real<lower=0> TDB_exponent; // mean function
  real<lower=0> TDB_offset; // mean function 
  real<lower=0> TDB_noise_intercept; // sd function
  real<lower=0> TDB_noise_slope; // sd function
  array[TDB_missing] real<lower=0> TDB_y_missing; // missing data parameter

  ///FBDL
  real<lower=0> FBDL_constant; // mean function
  real<lower=0> FBDL_exponent; // mean function
  real<lower=0> FBDL_offset; // mean function 
  real<lower=0> FBDL_noise_intercept; // sd function
  real<lower=0> FBDL_noise_slope; // sd function
  array[FBDL_missing] real<lower=0> FBDL_y_missing; // missing data parameter

  ///HDL
  real<lower=0> HDL_constant; // mean function
  real<lower=0> HDL_exponent; // mean function
  real<lower=0> HDL_offset; // mean function
  real<lower=0> HDL_noise_intercept; // sd function 
  real<lower=0> HDL_noise_slope; // sd function
  array[HDL_missing] real<lower=0> HDL_y_missing; // missing data parameter

  ///HPB
  real<lower=0> HPB_constant; // mean function
  real<lower=0> HPB_exponent; // mean function
  real<lower=0> HPB_offset; // mean function 
  real<lower=0> HPB_noise_intercept; // sd function
  real<lower=0> HPB_noise_slope; // sd function
  array[HPB_missing] real<lower=0> HPB_y_missing; // missing data parameter

  ///HMSB
  real<lower=0> HMSB_constant; // mean function
  real<lower=0> HMSB_exponent; // mean function
  real<lower=0> HMSB_offset; // mean function
  real<lower=0> HMSB_sigma; // sd function 
  array[HMSB_missing] real<lower=0> HMSB_y_missing; // missing data parameter

  ///HDB
  real<lower=0> HDB_constant; // mean function
  real<lower=0> HDB_exponent; // mean function
  real<lower=0> HDB_offset; // mean function 
  real<lower=0> HDB_noise_intercept; // sd function
  real<lower=0> HDB_noise_slope; // sd function
  array[HDB_missing] real<lower=0> HDB_y_missing; // missing data parameter

  ///RDL
  real<lower=0> RDL_constant; // mean function
  real<lower=0> RDL_exponent; // mean function
  real<lower=0> RDL_offset; // mean function
  real<lower=0> RDL_sigma; // sd function 
  array[RDL_missing] real<lower=0> RDL_y_missing; // missing data parameter

  ///RPB
  real<lower=0> RPB_constant; // mean function
  real<lower=0> RPB_exponent; // mean function
  real<lower=0> RPB_offset; // mean function 
  real<lower=0> RPB_noise_intercept; // sd function
  real<lower=0> RPB_noise_slope; // sd function
  array[RPB_missing] real<lower=0> RPB_y_missing; // missing data parameter
  
  ///RMSB
  real<lower=0> RMSB_constant; // mean function
  real<lower=0> RMSB_exponent; // mean function
  real<lower=0> RMSB_offset; // mean function 
  real<lower=0> RMSB_sigma; // sd function
  array[RMSB_missing] real<lower=0> RMSB_y_missing; // missing data parameter
  
  ///RDB
  real<lower=0> RDB_constant; // mean function
  real<lower=0> RDB_exponent; // mean function
  real<lower=0> RDB_offset; // mean function 
  real<lower=0> RDB_noise_intercept; // sd function
  real<lower=0> RDB_noise_slope; // sd function
  array[RDB_missing] real<lower=0> RDB_y_missing; // missing data parameter
  
  ///UDL
  real<lower=0> UDL_constant; // mean function
  real<lower=0> UDL_exponent; // mean function
  real<lower=0> UDL_offset; // mean function 
  real<lower=0> UDL_noise_intercept; // sd function
  real<lower=0> UDL_noise_slope; // sd function
  array[UDL_missing] real<lower=0> UDL_y_missing; // missing data parameter
  
  ///UMSB
  real<lower=0> UMSB_constant; // mean function
  real<lower=0> UMSB_exponent; // mean function
  real<lower=0> UMSB_offset; // mean function
  real<lower=0> UMSB_sigma; //sd function 
  array[UMSB_missing] real<lower=0> UMSB_y_missing; // missing data parameter

  // Ordinal Parameters
  
  ///Dentition
  
  ///max_m1
  real<lower=0> max_m1_constant;
  array[N] real<lower=0, upper=1> max_m1_u; // nuisance pars
  ordered[t_dental] max_m1_t_pars; // threshold pars

  ///max_m2
  real<lower=0> max_m2_constant;
  array[N] real<lower=0, upper=1> max_m2_u; // nuisance pars
  ordered[t_dental] max_m2_t_pars; // threshold pars

  ///max_m3
  real<lower=0> max_m3_constant;
  array[N] real<lower=0, upper=1> max_m3_u; // nuisance pars
  ordered[t_dental] max_m3_t_pars; // threshold pars

  ///max_pm1
  real<lower=0> max_pm1_constant;
  array[N] real<lower=0, upper=1> max_pm1_u; // nuisance pars
  ordered[t_dental] max_pm1_t_pars; // threshold pars

  ///max_pm2
  real<lower=0> max_pm2_constant;
  array[N] real<lower=0, upper=1> max_pm2_u; // nuisance pars
  ordered[t_dental] max_pm2_t_pars; // threshold pars

  ///max_c
  real<lower=0> max_c_constant;
  array[N] real<lower=0, upper=1> max_c_u; // nuisance pars
  ordered[t_dental] max_c_t_pars; // threshold pars

  ///max_i1
  real<lower=0> max_i1_constant;
  array[N] real<lower=0, upper=1> max_i1_u; // nuisance pars
  ordered[t_dental] max_i1_t_pars; // threshold pars

  ///max_i2
  real<lower=0> max_i2_constant;
  array[N] real<lower=0, upper=1> max_i2_u; // nuisance pars
  ordered[t_dental] max_i2_t_pars; // threshold pars

  ///man_m1
  real<lower=0> man_m1_constant;
  array[N] real<lower=0, upper=1> man_m1_u; // nuisance pars
  ordered[t_dental] man_m1_t_pars; // threshold pars

  ///man_m2
  real<lower=0> man_m2_constant;
  array[N] real<lower=0, upper=1> man_m2_u; // nuisance pars
  ordered[t_dental] man_m2_t_pars; // threshold pars
 
  ///man_m3
  real<lower=0> man_m3_constant;
  array[N] real<lower=0, upper=1> man_m3_u; // nuisance pars
  ordered[t_dental] man_m3_t_pars; // threshold pars 

  ///man_pm1
  real<lower=0> man_pm1_constant;
  array[N] real<lower=0, upper=1> man_pm1_u; // nuisance pars
  ordered[t_dental] man_pm1_t_pars; // threshold pars

  ///man_pm2
  real<lower=0> man_pm2_constant;
  array[N] real<lower=0, upper=1> man_pm2_u; // nuisance pars
  ordered[t_dental] man_pm2_t_pars; // threshold pars

  ///man_c
  real<lower=0> man_c_constant;
  array[N] real<lower=0, upper=1> man_c_u; // nuisance pars
  ordered[t_dental] man_c_t_pars; // threshold pars

  ///man_i1
  real<lower=0> man_i1_constant;
  array[N] real<lower=0, upper=1> man_i1_u; // nuisance pars
  ordered[t_dental] man_i1_t_pars; // threshold pars
  
  ///man_i2
  real<lower=0> man_i2_constant;
  array[N] real<lower=0, upper=1> man_i2_u; // nuisance pars
  ordered[t_dental] man_i2_t_pars; // threshold pars

  ///FH_EF
  real<lower=0> fh_ef_constant;
  array[N] real<lower=0, upper=1> fh_ef_u; // nuisance pars
  ordered[t_lb_ef] fh_ef_t_pars; // threshold pars

  ///FGT_EF
  real<lower=0> fgt_ef_constant;
  array[N] real<lower=0, upper=1> fgt_ef_u; // nuisance pars
  ordered[t_lb_ef] fgt_ef_t_pars; // threshold pars

  ///FLT_EF
  real<lower=0> flt_ef_constant;
  array[N] real<lower=0, upper=1> flt_ef_u; // nuisance pars
  ordered[t_lb_ef] flt_ef_t_pars; // threshold pars

  ///FDE_EF
  real<lower=0> fde_ef_constant;
  array[N] real<lower=0, upper=1> fde_ef_u; // nuisance pars
  ordered[t_lb_ef] fde_ef_t_pars; // threshold pars

  ///TPE_EF
  real<lower=0> tpe_ef_constant;
  array[N] real<lower=0, upper=1> tpe_ef_u; // nuisance pars
  ordered[t_lb_ef] tpe_ef_t_pars; // threshold pars

  ///TDE_EF
  real<lower=0> tde_ef_constant;
  array[N] real<lower=0, upper=1> tde_ef_u; // nuisance pars
  ordered[t_lb_ef] tde_ef_t_pars; // threshold pars

  ///FBPE_EF
  real<lower=0> fbpe_ef_constant;
  array[N] real<lower=0, upper=1> fbpe_ef_u; // nuisance pars
  ordered[t_lb_ef] fbpe_ef_t_pars; // threshold pars
  
  ///FBDE_EF
  real<lower=0> fbde_ef_constant;
  array[N] real<lower=0, upper=1> fbde_ef_u; // nuisance pars
  ordered[t_lb_ef] fbde_ef_t_pars; // threshold pars
  
  ///HPE_EF
  real<lower=0> hpe_ef_constant;
  array[N] real<lower=0, upper=1> hpe_ef_u; // nuisance pars
  ordered[t_lb_ef] hpe_ef_t_pars; // threshold pars
  
  ///HDE_EF
  real<lower=0> hde_ef_constant;
  array[N] real<lower=0, upper=1> hde_ef_u; // nuisance pars
  ordered[t_lb_ef] hde_ef_t_pars; // threshold pars
  
  ///HME_EF
  real<lower=0> hme_ef_constant;
  array[N] real<lower=0, upper=1> hme_ef_u; // nuisance pars
  ordered[t_lb_ef] hme_ef_t_pars; // threshold pars
  
  ///RPE_EF
  real<lower=0> rpe_ef_constant;
  array[N] real<lower=0, upper=1> rpe_ef_u; // nuisance pars
  ordered[t_lb_ef] rpe_ef_t_pars; // threshold pars
  
  ///RDE_EF
  real<lower=0> rde_ef_constant;
  array[N] real<lower=0, upper=1> rde_ef_u; // nuisance pars
  ordered[t_lb_ef] rde_ef_t_pars; // threshold pars
  
  ///UPE_EF
  real<lower=0> upe_ef_constant;
  array[N] real<lower=0, upper=1> upe_ef_u; // nuisance pars
  ordered[t_lb_ef] upe_ef_t_pars; // threshold pars
  
  ///UDE_EF
  real<lower=0> ude_ef_constant;
  array[N] real<lower=0, upper=1> ude_ef_u; // nuisance pars
  ordered[t_lb_ef] ude_ef_t_pars; // threshold pars
  
  ///CT_EF
  real<lower=0> ct_ef_constant;
  array[N] real<lower=0, upper=1> ct_ef_u; // nuisance pars
  ordered[t_lb_ef] ct_ef_t_pars; // threshold pars

  ///ISPR_EF
  real<lower=0> ispr_ef_constant;
  array[N] real<lower=0, upper=1> ispr_ef_u; // nuisance pars
  ordered[t_pelvis_ef] ispr_ef_t_pars; // threshold pars
  
  ///ILIS_EF
  real<lower=0> ilis_ef_constant;
  array[N] real<lower=0, upper=1> ilis_ef_u; // nuisance pars
  ordered[t_pelvis_ef] ilis_ef_t_pars; // threshold pars

  /// Carpal and Tarsal Ossification
  
  ///CC_Oss
  real<lower=0> cc_oss_constant;
  array[N] real<lower=0, upper=1> cc_oss_u; // nuisance pars
  ordered[t_carpal] cc_oss_t_pars; // threshold pars
  
  ///TC_Oss
  real<lower=0> tc_oss_constant;
  array[N] real<lower=0, upper=1> tc_oss_u; // nuisance pars
  ordered[t_tarsal] tc_oss_t_pars; // threshold pars  
  
}
transformed parameters{

  // Continuous Parameters
  array[N] real FDL_mean;
  array[N] real FDL_noise;
  array[N] real FMSB_mean;
  array[N] real FMSB_noise;
  array[N] real FDB_mean;
  array[N] real FDB_noise;
  array[N] real TDL_mean;
  array[N] real TDL_noise;
  array[N] real TPB_mean;
  array[N] real TPB_noise;
  array[N] real TMSB_mean;
  array[N] real TMSB_noise;
  array[N] real TDB_mean;
  array[N] real TDB_noise;
  array[N] real FBDL_mean;
  array[N] real FBDL_noise;
  array[N] real HDL_mean;
  array[N] real HDL_noise;
  array[N] real HPB_mean;
  array[N] real HPB_noise;
  array[N] real HMSB_mean;
  array[N] real HMSB_noise;
  array[N] real HDB_mean;
  array[N] real HDB_noise;
  array[N] real RDL_mean;
  array[N] real RDL_noise;
  array[N] real RPB_mean;
  array[N] real RPB_noise;
  array[N] real RMSB_mean;
  array[N] real RMSB_noise;
  array[N] real RDB_mean;
  array[N] real RDB_noise;
  array[N] real UDL_mean;
  array[N] real UDL_noise;
  array[N] real UMSB_mean;
  array[N] real UMSB_noise;

  // Dentition Parameters
  array[N] real max_m1_mean;
  array[N] real max_m2_mean;
  array[N] real max_m3_mean;
  array[N] real max_pm1_mean;
  array[N] real max_pm2_mean;
  array[N] real max_c_mean;
  array[N] real max_i1_mean;
  array[N] real max_i2_mean;
  array[N] real man_m1_mean;
  array[N] real man_m2_mean;
  array[N] real man_m3_mean;
  array[N] real man_pm1_mean;
  array[N] real man_pm2_mean;
  array[N] real man_c_mean;
  array[N] real man_i1_mean;
  array[N] real man_i2_mean;

  // LB ef parameters
  array[N] real fh_ef_mean;
  array[N] real fgt_ef_mean;
  array[N] real flt_ef_mean;
  array[N] real fde_ef_mean;
  array[N] real tpe_ef_mean;
  array[N] real tde_ef_mean;
  array[N] real fbpe_ef_mean;
  array[N] real fbde_ef_mean;
  array[N] real hpe_ef_mean;
  array[N] real hde_ef_mean;
  array[N] real hme_ef_mean;
  array[N] real rpe_ef_mean;
  array[N] real rde_ef_mean;
  array[N] real upe_ef_mean;
  array[N] real ude_ef_mean;
  array[N] real ct_ef_mean;

  // Pelvis ef parameters
  array[N] real ispr_ef_mean;
  array[N] real ilis_ef_mean;

  // Ossification parameters
  array[N] real cc_oss_mean;
  array[N] real tc_oss_mean;
  
  // Mean and noise estimation
  
  for(n in 1:N){
    
    // Continuous parameters
    FDL_mean[n] = FDL_constant*pow(x[n],FDL_exponent) + FDL_offset;
    FDL_noise[n] = FDL_noise_intercept + x[n]*FDL_noise_slope;
    FMSB_mean[n] = FMSB_constant*pow(x[n],FMSB_exponent) + FMSB_offset;
    FMSB_noise[n] = FMSB_noise_intercept + x[n]*FMSB_noise_slope;
    FDB_mean[n] = FDB_constant*pow(x[n],FDB_exponent) + FDB_offset;
    FDB_noise[n] = FDB_noise_intercept + x[n]*FDB_noise_slope;
    TDL_mean[n] = TDL_constant*x[n]^TDL_exponent + TDL_offset;
    TDL_noise[n] = TDL_noise_intercept + x[n]*TDL_noise_slope;
    TPB_mean[n] = TPB_constant*pow(x[n],TPB_exponent) + TPB_offset;
    TPB_noise[n] = TPB_noise_intercept + x[n]*TPB_noise_slope;
    TMSB_mean[n] = TMSB_constant*pow(x[n],TMSB_exponent) + TMSB_offset;
    TMSB_noise[n] = TMSB_sigma;
    TDB_mean[n] = TDB_constant*pow(x[n],TDB_exponent) + TDB_offset;
    TDB_noise[n] = TDB_noise_intercept + x[n]*TDB_noise_slope;
    FBDL_mean[n] = FBDL_constant*pow(x[n],FBDL_exponent) + FBDL_offset;
    FBDL_noise[n] = FBDL_noise_intercept + x[n]*FBDL_noise_slope;
    HDL_mean[n] = HDL_constant*pow(x[n],HDL_exponent) + HDL_offset;
    HDL_noise[n] = HDL_noise_intercept + x[n]*HDL_noise_slope;
    HPB_mean[n] = HPB_constant*pow(x[n],HPB_exponent) + HPB_offset;
    HPB_noise[n] = HPB_noise_intercept + x[n]*HPB_noise_slope;
    HMSB_mean[n] = HMSB_constant*pow(x[n],HMSB_exponent) + HMSB_offset;
    HMSB_noise[n] = HMSB_sigma;
    HDB_mean[n] = HDB_constant*pow(x[n],HDB_exponent) + HDB_offset;
    HDB_noise[n] = HDB_noise_intercept + x[n]*HDB_noise_slope;
    RDL_mean[n] = RDL_constant*pow(x[n],RDL_exponent) + RDL_offset;
    RDL_noise[n] = RDL_sigma;
    RPB_mean[n] = RPB_constant*pow(x[n],RPB_exponent) + RPB_offset;
    RPB_noise[n] = RPB_noise_intercept + x[n]*RPB_noise_slope;
    RMSB_mean[n] = RMSB_constant*pow(x[n],RMSB_exponent) + RMSB_offset;
    RMSB_noise[n] = RMSB_sigma;
    RDB_mean[n] = RDB_constant*pow(x[n],RDB_exponent) + RDB_offset;
    RDB_noise[n] = RDB_noise_intercept + x[n]*RDB_noise_slope;
    UDL_mean[n] = UDL_constant*pow(x[n],UDL_exponent) + UDL_offset;
    UDL_noise[n] = UDL_noise_intercept + x[n]*UDL_noise_slope;
    UMSB_mean[n] = UMSB_constant*pow(x[n],UMSB_exponent) + UMSB_offset;
    UMSB_noise[n] = UMSB_sigma;

    // Ordinal parameters
    max_m1_mean[n] = max_m1_constant*x[n];
    max_m2_mean[n] = max_m2_constant*x[n];
    max_m3_mean[n] = max_m3_constant*x[n];
    max_pm1_mean[n] = max_pm1_constant*x[n];
    max_pm2_mean[n] = max_pm2_constant*x[n];
    max_c_mean[n] = max_c_constant*x[n];
    max_i1_mean[n] = max_i1_constant*x[n];
    max_i2_mean[n] = max_i2_constant*x[n];
    man_m1_mean[n] = man_m1_constant*x[n];
    man_m2_mean[n] = man_m2_constant*x[n];
    man_m3_mean[n] = man_m3_constant*x[n];
    man_pm1_mean[n] = man_pm1_constant*x[n];
    man_pm2_mean[n] = man_pm2_constant*x[n];
    man_c_mean[n] = man_c_constant*x[n];
    man_i1_mean[n] = man_i1_constant*x[n];
    man_i2_mean[n] = man_i2_constant*x[n];
    fh_ef_mean[n] = fh_ef_constant*x[n];
    fgt_ef_mean[n] = fgt_ef_constant*x[n];
    flt_ef_mean[n] = flt_ef_constant*x[n];
    fde_ef_mean[n] = fde_ef_constant*x[n];
    tpe_ef_mean[n] = tpe_ef_constant*x[n];
    tde_ef_mean[n] = tde_ef_constant*x[n];
    fbpe_ef_mean[n] = fbpe_ef_constant*x[n];
    fbde_ef_mean[n] = fbde_ef_constant*x[n];
    hpe_ef_mean[n] = hpe_ef_constant*x[n];
    hde_ef_mean[n] = hde_ef_constant*x[n];
    hme_ef_mean[n] = hme_ef_constant*x[n];
    rpe_ef_mean[n] = rpe_ef_constant*x[n];
    rde_ef_mean[n] = rde_ef_constant*x[n];
    upe_ef_mean[n] = upe_ef_constant*x[n];
    ude_ef_mean[n] = ude_ef_constant*x[n];
    ct_ef_mean[n] = ct_ef_constant*x[n];
    ispr_ef_mean[n] = ispr_ef_constant*x[n];
    ilis_ef_mean[n] = ilis_ef_constant*x[n];
    cc_oss_mean[n] = cc_oss_constant*x[n];
    tc_oss_mean[n] = tc_oss_constant*x[n];
    
  }
  
}
model{
  
  // Complete continuous vectors with missing parameters
  array[N] real y_all_FDL;
  array[N] real y_all_FMSB;
  array[N] real y_all_FDB;
  array[N] real y_all_TDL;
  array[N] real y_all_TPB;
  array[N] real y_all_TMSB;
  array[N] real y_all_TDB;
  array[N] real y_all_FBDL;
  array[N] real y_all_HDL;
  array[N] real y_all_HPB;
  array[N] real y_all_HMSB;
  array[N] real y_all_HDB;
  array[N] real y_all_RDL;
  array[N] real y_all_RPB;
  array[N] real y_all_RMSB;
  array[N] real y_all_RDB;
  array[N] real y_all_UDL;
  array[N] real y_all_UMSB;

  
  ///FDL
  for(n in 1:FDL_complete){
    y_all_FDL[FDL_index_present[n]] = y_FDL[n];
  }
  for(n in 1:FDL_missing){
    y_all_FDL[FDL_index_missing[n]] = FDL_y_missing[n];
  }
  
  ///FMSB
  for(n in 1:FMSB_complete){
    y_all_FMSB[FMSB_index_present[n]] = y_FMSB[n];
  }
  for(n in 1:FMSB_missing){
    y_all_FMSB[FMSB_index_missing[n]] = FMSB_y_missing[n];
  }

  ///FDB
  for(n in 1:FDB_complete){
    y_all_FDB[FDB_index_present[n]] = y_FDB[n];
  }
  for(n in 1:FDB_missing){
    y_all_FDB[FDB_index_missing[n]] = FDB_y_missing[n];
  }

  ///TDL
  for(n in 1:TDL_complete){
    y_all_TDL[TDL_index_present[n]] = y_TDL[n];
  }
  for(n in 1:TDL_missing){
    y_all_TDL[TDL_index_missing[n]] = TDL_y_missing[n];
  }

  ///TPB
  for(n in 1:TPB_complete){
    y_all_TPB[TPB_index_present[n]] = y_TPB[n];
  }
  for(n in 1:TPB_missing){
    y_all_TPB[TPB_index_missing[n]] = TPB_y_missing[n];
  }

  ///TMSB
  for(n in 1:TMSB_complete){
    y_all_TMSB[TMSB_index_present[n]] = y_TMSB[n];
  }
  for(n in 1:TMSB_missing){
    y_all_TMSB[TMSB_index_missing[n]] = TMSB_y_missing[n];
  }

  ///TDB
  for(n in 1:TDB_complete){
    y_all_TDB[TDB_index_present[n]] = y_TDB[n];
  }
  for(n in 1:TDB_missing){
    y_all_TDB[TDB_index_missing[n]] = TDB_y_missing[n];
  }

  ///FBDL
  for(n in 1:FBDL_complete){
    y_all_FBDL[FBDL_index_present[n]] = y_FBDL[n];
  }
  for(n in 1:FBDL_missing){
    y_all_FBDL[FBDL_index_missing[n]] = FBDL_y_missing[n];
  }

  ///HDL
  for(n in 1:HDL_complete){
    y_all_HDL[HDL_index_present[n]] = y_HDL[n];
  }
  for(n in 1:HDL_missing){
    y_all_HDL[HDL_index_missing[n]] = HDL_y_missing[n];
  }

  ///HPB
  for(n in 1:HPB_complete){
    y_all_HPB[HPB_index_present[n]] = y_HPB[n];
  }
  for(n in 1:HPB_missing){
    y_all_HPB[HPB_index_missing[n]] = HPB_y_missing[n];
  }

  ///HMSB
  for(n in 1:HMSB_complete){
    y_all_HMSB[HMSB_index_present[n]] = y_HMSB[n];
  }
  for(n in 1:HMSB_missing){
    y_all_HMSB[HMSB_index_missing[n]] = HMSB_y_missing[n];
  }

  ///HDB
  for(n in 1:HDB_complete){
    y_all_HDB[HDB_index_present[n]] = y_HDB[n];
  }
  for(n in 1:HDB_missing){
    y_all_HDB[HDB_index_missing[n]] = HDB_y_missing[n];
  }

  ///RDL
  for(n in 1:RDL_complete){
    y_all_RDL[RDL_index_present[n]] = y_RDL[n];
  }
  for(n in 1:RDL_missing){
    y_all_RDL[RDL_index_missing[n]] = RDL_y_missing[n];
  }

  ///RPB
  for(n in 1:RPB_complete){
    y_all_RPB[RPB_index_present[n]] = y_RPB[n];
  }
  for(n in 1:RPB_missing){
    y_all_RPB[RPB_index_missing[n]] = RPB_y_missing[n];
  }
  
  ///RMSB
  for(n in 1:RMSB_complete){
    y_all_RMSB[RMSB_index_present[n]] = y_RMSB[n];
  }
  for(n in 1:RMSB_missing){
    y_all_RMSB[RMSB_index_missing[n]] = RMSB_y_missing[n];
  }
  
  ///RDB
  for(n in 1:RDB_complete){
    y_all_RDB[RDB_index_present[n]] = y_RDB[n];
  }
  for(n in 1:RDB_missing){
    y_all_RDB[RDB_index_missing[n]] = RDB_y_missing[n];
  }
  
  ///UDL
  for(n in 1:UDL_complete){
    y_all_UDL[UDL_index_present[n]] = y_UDL[n];
  }
  for(n in 1:UDL_missing){
    y_all_UDL[UDL_index_missing[n]] = UDL_y_missing[n];
  }
  
  ///UMSB
  for(n in 1:UMSB_complete){
    y_all_UMSB[UMSB_index_present[n]] = y_UMSB[n];
  }
  for(n in 1:UMSB_missing){
    y_all_UMSB[UMSB_index_missing[n]] = UMSB_y_missing[n];
  }
  
  // Priors
  
  /// Global Priors
  L ~ lkj_corr_cholesky(50.0);
  
  /// Continuous Priors
  
  ///FDL
  FDL_constant ~ normal(65,1); 
  FDL_exponent ~ std_normal(); 
  FDL_offset ~ normal(65,1);  
  FDL_noise_intercept ~ normal(6,1); 
  FDL_noise_slope ~ std_normal();
  
  
  ///FMSB
  FMSB_constant ~ normal(5,1); 
  FMSB_exponent ~ std_normal(); 
  FMSB_offset ~ normal(6,1);  
  FMSB_noise_intercept ~ std_normal(); 
  FMSB_noise_slope ~ std_normal();

  ///FDB
  FDB_constant ~ normal(18,1); 
  FDB_exponent ~ std_normal(); 
  FDB_offset ~ normal(16,1);  
  FDB_noise_intercept ~ normal(3,1); 
  FDB_noise_slope ~ std_normal();  

  ///TDL
  TDL_constant ~ normal(50,1); 
  TDL_exponent ~ std_normal(); 
  TDL_offset ~ normal(58,1);  
  TDL_noise_intercept ~ normal(6,1); 
  TDL_noise_slope ~ normal(1,1); 

  ///TPB
  TPB_constant ~ normal(14,1); 
  TPB_exponent ~ std_normal(); 
  TPB_offset ~ normal(14,1);  
  TPB_noise_intercept ~ normal(2,1); 
  TPB_noise_slope ~ std_normal();


  ///TMSB
  TMSB_constant ~ normal(4,1); 
  TMSB_exponent ~ std_normal(); 
  TMSB_offset ~ normal(6,1); 
  TMSB_sigma ~ std_normal();

  ///TDB
  TDB_constant ~ normal(9,1); 
  TDB_exponent ~ std_normal(); 
  TDB_offset ~ normal(9,1);  
  TDB_noise_intercept ~ std_normal(); 
  TDB_noise_slope ~ std_normal();

  ///FBDL
  FBDL_constant ~ normal(52,1); 
  FBDL_exponent ~ normal(0,1); 
  FBDL_offset ~ normal(52,1);  
  FBDL_noise_intercept ~ normal(6,1); 
  FBDL_noise_slope ~ normal(2,1); 

  ///HDL
  HDL_constant ~ normal(46,1); 
  HDL_exponent ~ std_normal(); 
  HDL_noise_slope ~ normal(5,1); 
  HDL_offset ~ normal(57,1);
  HDL_noise_intercept ~ normal(5,1);

  ///HPB
  HPB_constant ~ normal(9, 1); 
  HPB_exponent ~ std_normal(); 
  HPB_offset ~ normal(11,1);  
  HPB_noise_intercept ~ std_normal(); 
  HPB_noise_slope ~ std_normal();

  ///HMSB
  HMSB_constant ~ normal(6,1); 
  HMSB_exponent ~ std_normal(); 
  HMSB_sigma ~ std_normal(); 
  HMSB_offset ~ normal(4,1);

  ///HDB
  HDB_constant ~ normal(12,1); 
  HDB_exponent ~ std_normal(); 
  HDB_offset ~ normal(14,1);  
  HDB_noise_intercept ~ std_normal(); 
  HDB_noise_slope ~ std_normal(); 

  ///RDL
  RDL_constant ~ normal(33,2); 
  RDL_exponent ~ normal(0.5,0.5);  
  RDL_offset ~ normal(46,1);
  RDL_sigma ~ student_t(3,0,2.5);


  ///RPB
  RPB_constant ~ normal(3,1); 
  RPB_exponent ~ normal(0,1); 
  RPB_offset ~ normal(5,1);  
  RPB_noise_intercept ~ std_normal(); 
  RPB_noise_slope ~ std_normal();
  
  ///RMSB
  RMSB_constant ~ normal(4,1); 
  RMSB_exponent ~ std_normal(); 
  RMSB_offset ~ normal(4,1);  
  RMSB_sigma ~ std_normal(); 
  
  ///RDB
  RDB_constant ~ normal(5,1); 
  RDB_exponent ~ std_normal(); 
  RDB_offset ~ normal(10,1);  
  RDB_noise_intercept ~ std_normal(); 
  RDB_noise_slope ~ std_normal(); 
  
  ///UDL
  UDL_constant ~ normal(37,1); 
  UDL_exponent ~ normal(0,1); 
  UDL_offset ~ normal(54,1);  
  UDL_noise_intercept ~ std_normal(); 
  UDL_noise_slope ~ normal(5,1); 
  
  ///UMSB
  UMSB_constant ~ normal(4,1); 
  UMSB_exponent ~ std_normal(); 
  UMSB_offset ~ normal(4,1);
  UMSB_sigma ~ std_normal();

  // Ordinal Priors
	
  ///Dentition
  
  ///max_m1
  max_m1_constant ~ std_normal();
  
  for(i in 1:size(max_m1_t_pars)){
    max_m1_t_pars[i] ~ normal(i+1,1);
  }

  ///max_m2
  max_m2_constant ~ std_normal();
  
  for(i in 1:size(max_m2_t_pars)){
    max_m2_t_pars[i] ~ normal(i+1,1);
  }

  ///max_m3
  max_m3_constant ~ std_normal();
  
  for(i in 1:size(max_m3_t_pars)){
    max_m3_t_pars[i] ~ normal(i+1,1);
  }

  ///max_pm1
  max_pm1_constant ~ std_normal();
  
  for(i in 1:size(max_pm1_t_pars)){
    max_pm1_t_pars[i] ~ normal(i+1,1);
  }

  ///max_pm2
  max_pm2_constant ~ std_normal();
  
  for(i in 1:size(max_pm2_t_pars)){
    max_pm2_t_pars[i] ~ normal(i+1,1);
  }

  ///max_c
  max_c_constant ~ std_normal();
  
  for(i in 1:size(max_c_t_pars)){
    max_c_t_pars[i] ~ normal(i+1,1);
  }

  ///max_i1
  max_i1_constant ~ std_normal();
  
  for(i in 1:size(max_i1_t_pars)){
    max_i1_t_pars[i] ~ normal(i+1,1);
  }

  ///max_i2
  max_i2_constant ~ std_normal();
  
  for(i in 1:size(max_i2_t_pars)){
    max_i2_t_pars[i] ~ normal(i+1,1);
  }

  ///man_m1
  man_m1_constant ~ std_normal();
  
  for(i in 1:size(man_m1_t_pars)){
    man_m1_t_pars[i] ~ normal(i+1,1);
  }

  ///man_m2
  man_m2_constant ~ std_normal();
  
  for(i in 1:size(man_m2_t_pars)){
    man_m2_t_pars[i] ~ normal(i+1,1);
  }

  ///man_m3
  man_m3_constant ~ std_normal();

  for(i in 1:size(man_m3_t_pars)){
    man_m3_t_pars[i] ~ normal(i+1,1);
  }

  ///man_pm1
  man_pm1_constant ~ std_normal();
  
  for(i in 1:size(man_pm1_t_pars)){
    man_pm1_t_pars[i] ~ normal(i+1,1);
  }

  ///man_pm2
  man_pm2_constant ~ std_normal();
  
  for(i in 1:size(man_pm2_t_pars)){
    man_pm2_t_pars[i] ~ normal(i+1,1);
  }

  ///man_c
  man_c_constant ~ std_normal();
  
  for(i in 1:size(man_c_t_pars)){
    man_c_t_pars[i] ~ normal(i+1,1);
  }

  ///man_i1
  man_i1_constant ~ std_normal();
  
  for(i in 1:size(man_i1_t_pars)){
    man_i1_t_pars[i] ~ normal(i+1,1);
  }
  
  ///man_i2
  man_i2_constant ~ std_normal();
  
  for(i in 1:size(man_i2_t_pars)){
    man_i2_t_pars[i] ~ normal(i+1,1);
  }

  ///lb ef
  
  ///FH_EF
  fh_ef_constant ~ std_normal();
  
  for(i in 1:size(fh_ef_t_pars)){
    fh_ef_t_pars[i] ~ normal(i+1,1);
  }

  ///FGT_EF
  fgt_ef_constant ~ std_normal();
  
  for(i in 1:size(fgt_ef_t_pars)){
    fgt_ef_t_pars[i] ~ normal(i+1,1);
  }

  ///FLT_EF
  flt_ef_constant ~ std_normal();
  
  for(i in 1:size(flt_ef_t_pars)){
    flt_ef_t_pars[i] ~ normal(i+1,1);
  }

  ///FDE_EF
  fde_ef_constant ~ std_normal();
  
  for(i in 1:size(fde_ef_t_pars)){
    fde_ef_t_pars[i] ~ normal(i+1,1);
  }

  ///TPE_EF
  tpe_ef_constant ~ std_normal();
  
  for(i in 1:size(tpe_ef_t_pars)){
    tpe_ef_t_pars[i] ~ normal(i+1,1);
  }

 ///TDE_EF
  tde_ef_constant ~ std_normal();
  
  for(i in 1:size(tde_ef_t_pars)){
    tde_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///FBPE_EF
  fbpe_ef_constant ~ std_normal();
  
  for(i in 1:size(fbpe_ef_t_pars)){
    fbpe_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///FBDE_EF
  fbde_ef_constant ~ std_normal();
  
  for(i in 1:size(fh_ef_t_pars)){
    fbde_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///HPE_EF
  hpe_ef_constant ~ std_normal();
  
  for(i in 1:size(hpe_ef_t_pars)){
    hpe_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///HDE_EF
  hde_ef_constant ~ std_normal();
  
  for(i in 1:size(hde_ef_t_pars)){
    hde_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///HME_EF
  hme_ef_constant ~ std_normal();
  
  for(i in 1:size(hme_ef_t_pars)){
    hme_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///RPE_EF
  rpe_ef_constant ~ std_normal();
  
  for(i in 1:size(rpe_ef_t_pars)){
    rpe_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///RDE_EF
  rde_ef_constant ~ std_normal();
  
  for(i in 1:size(rde_ef_t_pars)){
    rde_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///UPE_EF
  upe_ef_constant ~ std_normal();
  
  for(i in 1:size(upe_ef_t_pars)){
    upe_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///UDE_EF
  ude_ef_constant ~ std_normal();
  
  for(i in 1:size(ude_ef_t_pars)){
    ude_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///CT_EF
  ct_ef_constant ~ std_normal();
  
  for(i in 1:size(ct_ef_t_pars)){
    ct_ef_t_pars[i] ~ normal(i+1,1);
  }

  /// Pelvis Fusion
  
  ///ISPR_EF
  ispr_ef_constant ~ std_normal();
  
  for(i in 1:size(ispr_ef_t_pars)){
    ispr_ef_t_pars[i] ~ normal(i+1,1);
  }
  
  ///ILIS_EF
  ilis_ef_constant ~ std_normal();
  
  for(i in 1:size(ilis_ef_t_pars)){
    ilis_ef_t_pars[i] ~ normal(i+1,1);
  }

  /// Carpal and Tarsal Ossification
  
  ///CC_Oss
  cc_oss_constant ~ std_normal();
  
  for(i in 1:size(cc_oss_t_pars)){
    cc_oss_t_pars[i] ~ normal(i+1,1);
  }
  
  ///TC_Oss
  tc_oss_constant ~ std_normal();
  
  for(i in 1:size(tc_oss_t_pars)){
    tc_oss_t_pars[i] ~ normal(i+1,1);
  }

  
  // LIKELIHOOD
  
  target += centered_gaussian_copula_cholesky_(
    {normal_marginal(y_all_FDL, FDL_mean, FDL_noise, N),
    normal_marginal(y_all_FMSB, FMSB_mean, FMSB_noise, N),
    normal_marginal(y_all_FDB, FDB_mean, FDB_noise, N),
    normal_marginal(y_all_TDL, TDL_mean, TDL_noise, N),
    normal_marginal(y_all_TPB, TPB_mean, TPB_noise, N),
    normal_marginal(y_all_TMSB, TMSB_mean, TMSB_noise, N),
    normal_marginal(y_all_TDB, TDB_mean, TDB_noise, N),
    normal_marginal(y_all_FBDL, FBDL_mean, FBDL_noise, N),
    normal_marginal(y_all_HDL, HDL_mean, HDL_noise, N),
    normal_marginal(y_all_HPB, HPB_mean, HPB_noise, N),
    normal_marginal(y_all_HMSB, HMSB_mean, HMSB_noise, N),
    normal_marginal(y_all_HDB, HDB_mean, HDB_noise, N),
    normal_marginal(y_all_RDL, RDL_mean, RDL_noise, N),
    normal_marginal(y_all_RPB, RPB_mean, RPB_noise, N),
    normal_marginal(y_all_RMSB, RMSB_mean, RMSB_noise, N),
    normal_marginal(y_all_RDB, RDB_mean, RDB_noise, N),
    normal_marginal(y_all_UDL, UDL_mean, UDL_noise, N),
    normal_marginal(y_all_UMSB, UMSB_mean, UMSB_noise, N),
    probit_marginal(y_max_m1, max_m1_mean, max_m1_u, max_m1_t_pars),
    probit_marginal(y_max_m2, max_m2_mean, max_m2_u, max_m2_t_pars),
    probit_marginal(y_max_m3, max_m3_mean, max_m3_u, max_m3_t_pars),
    probit_marginal(y_max_pm1, max_pm1_mean, max_pm1_u, max_pm1_t_pars),
    probit_marginal(y_max_pm2, max_pm2_mean, max_pm2_u, max_pm2_t_pars),
    probit_marginal(y_max_c, max_c_mean, max_c_u, max_c_t_pars),
    probit_marginal(y_max_i1, max_i1_mean, max_i1_u, max_i1_t_pars),
    probit_marginal(y_max_i2, max_i2_mean, max_i2_u, max_i2_t_pars),
    probit_marginal(y_man_m1, man_m1_mean, man_m1_u, man_m1_t_pars),
    probit_marginal(y_man_m2, man_m2_mean, man_m2_u, man_m2_t_pars),
    probit_marginal(y_man_m3, man_m3_mean, man_m3_u, man_m3_t_pars),
    probit_marginal(y_man_pm1, man_pm1_mean, man_pm1_u, man_pm1_t_pars),
    probit_marginal(y_man_pm2, man_pm2_mean, man_pm2_u, man_pm2_t_pars),
    probit_marginal(y_man_c, man_c_mean, man_c_u, man_c_t_pars),
    probit_marginal(y_man_i1, man_i1_mean, man_i1_u, man_i1_t_pars),
    probit_marginal(y_man_i2, man_i2_mean, man_i2_u, man_i2_t_pars),
    probit_marginal(y_fh_ef, fh_ef_mean, fh_ef_u, fh_ef_t_pars),
    probit_marginal(y_fgt_ef, fgt_ef_mean, fgt_ef_u, fgt_ef_t_pars),
    probit_marginal(y_flt_ef, flt_ef_mean, flt_ef_u, flt_ef_t_pars),
    probit_marginal(y_fde_ef, fde_ef_mean, fde_ef_u, fde_ef_t_pars),
    probit_marginal(y_tpe_ef, tpe_ef_mean, tpe_ef_u, tpe_ef_t_pars),
    probit_marginal(y_tde_ef, tde_ef_mean, tde_ef_u, tde_ef_t_pars),
    probit_marginal(y_fbpe_ef, fbpe_ef_mean, fbpe_ef_u, fbpe_ef_t_pars),
    probit_marginal(y_fbde_ef, fbde_ef_mean, fbde_ef_u, fbde_ef_t_pars),
    probit_marginal(y_hpe_ef, hpe_ef_mean, hpe_ef_u, hpe_ef_t_pars),
    probit_marginal(y_hde_ef, hde_ef_mean, hde_ef_u, hde_ef_t_pars),
    probit_marginal(y_hme_ef, hme_ef_mean, hme_ef_u, hme_ef_t_pars),
    probit_marginal(y_rpe_ef, rpe_ef_mean, rpe_ef_u, rpe_ef_t_pars),
    probit_marginal(y_rde_ef, rde_ef_mean, rde_ef_u, rde_ef_t_pars),
    probit_marginal(y_upe_ef, upe_ef_mean, upe_ef_u, upe_ef_t_pars),
    probit_marginal(y_ude_ef, ude_ef_mean, ude_ef_u, ude_ef_t_pars),
    probit_marginal(y_ct_ef, ct_ef_mean, ct_ef_u, ct_ef_t_pars),
    probit_marginal(y_ispr_ef, ispr_ef_mean, ispr_ef_u, ispr_ef_t_pars),
    probit_marginal(y_ilis_ef, ilis_ef_mean, ilis_ef_u, ilis_ef_t_pars),
    probit_marginal(y_cc_oss, cc_oss_mean, cc_oss_u, cc_oss_t_pars),
    probit_marginal(y_tc_oss, tc_oss_mean, tc_oss_u, tc_oss_t_pars)}, L);

  
}
generated quantities{
  
  // Correlation Matrix
  
  corr_matrix[M] corr_mat = multiply_lower_tri_self_transpose(L);
  
}
