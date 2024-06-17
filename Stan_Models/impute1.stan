data{

  int N; // size of input
  int M; // size of full correlation matrix
  int K; // size of sample correlation matrix

  array[N] real y1;
  array[N] real y2;
  array[N] real x;

}
parameters{

  matrix[M,M] corr_mat;
  real FMSB_constant; // mean function
  real FMSB_exponent; // mean function
  real FMSB_offset; // mean function 
  real FMSB_noise_intercept; // sd function
  real FMSB_noise_slope; // sd function
  real FDB_constant; // mean function
  real FDB_exponent; // mean function
  real FDB_offset; // mean function 
  real FDB_noise_intercept; // sd function
  real FDB_noise_slope; // sd function
  real FDL_constant; // mean function
  real FDL_exponent; // mean function
  real FDL_offset; // mean function 
  real FDL_noise_intercept; // sd function
  real FDL_noise_slope; // sd function

}
transformed parameters{

  array[N] real mu_y1;
  array[N] real mu_y2;
  array[N] real mu_y3;
  array[N] real sd_y1;
  array[N] real sd_y2;
  array[N] real sd_y3;

  for(n in 1:N){
  
    mu_y1[n] = FMSB_constant*x[n]^FMSB_exponent + FMSB_offset;
    mu_y2[n] = FDB_constant*x[n]^FDB_exponent + FDB_offset;
    mu_y3[n] = FDL_constant*x[n]^FDL_exponent + FDL_offset;
    sd_y1[n] = FMSB_noise_intercept + x[n]*FMSB_noise_slope;
    sd_y2[n] = FDB_noise_intercept + x[n]*FDB_noise_slope;
    sd_y3[n] = FDL_noise_intercept + x[n]*FDL_noise_slope;
 
  
  }


  matrix [K,K] sub_cor;
  sub_cor = add_diag(sub_cor, 1);
  sub_cor[1,2] = corr_mat[3,2]; 
  sub_cor[1,3] = corr_mat[1,2] ;
  sub_cor[2,1] = corr_mat[3,2];
  sub_cor[3,1] = corr_mat[1,2]; 
  sub_cor[3,2] = corr_mat[1,3];
  sub_cor[2,3] = corr_mat[1,3];
}
generated quantities{

  array[N] real expectation;
  array[N] real ypred;

  {

    for(n in 1:N){

    real mu_cond = mu_y3[n] + sub_cor[1,3]*(sd_y3[n]/sd_y1[n])*(y1[n] - mu_y1[n]) + sub_cor[2,3]*(sd_y3[n]/sd_y2[n])*(y2[n] - mu_y2[n]);
    real sd_cond = sd_y3[n]*sqrt(1 - ((sub_cor[1,3]^2 + sub_cor[2,3]^2 - 2*sub_cor[1,3]*sub_cor[2,3]*sub_cor[1,2]) / (1-sub_cor[1,2]^2)));

    expectation[n] = mu_cond;
    ypred[n] = mu_cond + sd_cond*std_normal_rng();  

    }  

  }

}
