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

  real weibull_mixture_lpdf(array[] real x, vector lambda, vector alpha, vector beta) {
    int N = num_elements(x);
    int K = num_elements(lambda);
    real lp = 0.0;
    vector[K] log_lambda = log(lambda);
    
    for (n in 1:N) {
      vector[K] lps = log_lambda;
      for (k in 1:K) {
        lps[k] += weibull_lpdf(x[n]+0.002 | alpha[k], beta[k]);
      }
      lp += log_sum_exp(lps);
    }
    return lp;
  }

}
data{

  // Data Information

  int N;
  int M; // total number of response variables (columns)
  int LB_C; // long bone categories
  int DD_C; // dental categories

  // Response Variables

  array[N] real FDL;
  array[N] int HPE;
  array[N] int max_M1;

  // Model Parameters
  
  matrix[M,M] corr_mat;

  real FDL_constant; // mean function
  real FDL_exponent; // mean function
  real FDL_offset; // mean function 
  real FDL_noise_intercept; // sd function
  real FDL_noise_slope; // sd function

  real HPE_ef_constant;
  vector[LB_C - 1] HPE_ef_t_pars;

  real max_M1_constant;
  vector[DD_C - 1] max_M1_t_pars;

}
transformed data{

  array[N] real<lower=0, upper=1> hpe_u;
  array[N] real<lower=0, upper=1> m1_u;

  for(n in 1:N){

    hpe_u[n] = uniform_rng(0,1);
    m1_u[n] = uniform_rng(0,1);
  }

}
parameters{

  array[N] real<lower=0, upper=23> x;

}
transformed parameters{

  array[N] real FDL_mean;
  array[N] real FDL_noise;

  array[N] real HPE_mean;
  array[N] real M1_mean;

  for(n in 1:N){

    FDL_mean[n] = FDL_constant*x[n]^FDL_exponent + FDL_offset;
    FDL_noise[n] = FDL_noise_intercept + x[n]*FDL_noise_slope;

    HPE_mean[n] = HPE_ef_constant*x[n];
    M1_mean[n] = max_M1_constant*x[n];

  }

}
model{

  target += weibull_mixture_lpdf(x | [0.269, 0.497, 0.234]', [0.651,8.623,0.736]', [2.30,18.38,3.19]');

  target += centered_gaussian_copula_cholesky_(
    {normal_marginal(FDL, FDL_mean, FDL_noise, 1),
     probit_marginal(HPE, HPE_mean, hpe_u, HPE_ef_t_pars),
     probit_marginal(max_M1, M1_mean, m1_u, max_M1_t_pars)}, cholesky_decompose(corr_mat));

}
