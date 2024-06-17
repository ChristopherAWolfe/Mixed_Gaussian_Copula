data{

  int N; // size of input
  int M; // size of full correlation matrix
  int K; // size of sample correlation matrix

  array[N] real y1;
  array[N] int y2;
  array[N] real x;

  int y2_c;
  int y3_c;

}
transformed data{

  array[3] int v = {9,30,50}; //indices of matrix to subset

}
parameters{

  matrix[M,M] corr_mat;
  real HDL_constant; // mean function
  real HDL_exponent; // mean function
  real HDL_offset; // mean function 
  real HDL_noise_intercept; // sd function
  real HDL_noise_slope; // sd function
  real man_pm1_constant; // mean function
  vector[y2_c-1] man_pm1_t_pars; // mean function
  real ct_ef_constant; // mean function 
  vector[y3_c-1] ct_ef_t_pars; // sd function

}
transformed parameters{

  array[N] real mu_y1;
  array[N] real mu_y2;
  array[N] real mu_y3;
  array[N] real sd_y1;

  for(n in 1:N){
  
    mu_y1[n] = HDL_constant*x[n]^HDL_exponent + HDL_offset;
    mu_y2[n] = man_pm1_constant*x[n];
    mu_y3[n] = ct_ef_constant*x[n];
    sd_y1[n] = HDL_noise_intercept + x[n]*HDL_noise_slope; 
  
  }


  matrix [K,K] sub_cor = corr_mat[v,v];

}
generated quantities{

  array[N] int ypred;
  array[N] int expectation;

  {

    for(n in 1:N){

      real latent;
 
      if (y2[n] == 1) {
        latent = man_pm1_t_pars[1] - 1; // Below the first threshold
      } else if (y2[n] == y2_c) {
        latent = man_pm1_t_pars[y2_c-1] + 1; // Above the last threshold
      } else {
        latent = (man_pm1_t_pars[y2[n] - 1] + man_pm1_t_pars[y2[n]]) / 2;
      }
      
    real mu_cond = mu_y3[n] + sub_cor[1,3]*(1/sd_y1[n])*(y1[n] - mu_y1[n]) + sub_cor[2,3]*(1/1)*(latent - mu_y2[n]);
    real sd_cond = 1*sqrt(1 - ((sub_cor[1,3]^2 + sub_cor[2,3]^2 - 2*sub_cor[1,3]*sub_cor[2,3]*sub_cor[1,2]) / (1-sub_cor[1,2]^2)));

    real latent_pred = mu_cond + sd_cond*std_normal_rng();

    int category = 0;
    if (mu_cond <= ct_ef_t_pars[1]) {
      category = 1;
    } else if (mu_cond <= ct_ef_t_pars[2]) {
      category = 2;
    } else if (mu_cond <= ct_ef_t_pars[3]) {
      category = 3;
    } else if (mu_cond <= ct_ef_t_pars[4]) {
      category = 4;
    } else if (mu_cond <= ct_ef_t_pars[5]) {
      category = 5;
    } else if (mu_cond <= ct_ef_t_pars[6]) {
      category = 6;
    } else {
      category = 7;
    } 
    // Assign predicted value
    ypred[n] = category;

    int category1 = 0;
    if (latent_pred <= ct_ef_t_pars[1]) {
      category1 = 1;
    } else if (latent_pred <= ct_ef_t_pars[2]) {
      category1 = 2;
    } else if (latent_pred <= ct_ef_t_pars[3]) {
      category1 = 3;
    } else if (latent_pred <= ct_ef_t_pars[4]) {
      category1 = 4;
    } else if (latent_pred <= ct_ef_t_pars[5]) {
      category1 = 5;
    } else if (latent_pred <= ct_ef_t_pars[6]) {
      category1 = 6;
    } else {
      category1 = 7;
    } 
    // Assign predicted value
    expectation[n] = category1;

    }

    }

  

}
