
data {
  int<lower=1> N;  // total number of observations
  int<lower=2> ncolY;  // number of categories
  int<lower=2> ncolX; // number of predictor levels
  matrix[N,ncolX] X; // predictor design matrix
  int<lower=0,upper=3600> Y[N,ncolY];  // response variable
  
  int num_data; // length of Hour vector
  int num_basis; // nrow of bsplines
  
  int<lower=1> nd; //number of unique mice
  int<lower=1> ndays; //number of unique days
  int<lower=1> nbatch; //number of unique batches

  int<lower=1, upper=nd> mice_id[N]; // mice id full vector
  int<lower=1, upper=ndays> days[N]; // mice id full vector
  int<lower=1, upper=nbatch> batch[N]; // mice id full vector
  
  matrix[num_basis, num_data] Bsplines; // cyclic-splines computed in R
  
}

parameters {
  matrix[ncolY-1,ncolX] params_mu; // coefficients for mean mu
  matrix[ncolY-1,ncolX] params_sig; // coefficients for dispersion sigma
  matrix[ncolY-1,ncolX] params_pi; // coefficients for pi
  vector<lower=0,upper=1>[ncolY-1] Zed[N]; // independent beta random variables
  real<lower=0,upper=10> tau_id_mu;
  real<lower=0,upper=10> tau_days_mu;
  real<lower=0,upper=100> tau_batch_mu;
  real<lower=0,upper=10> tau_id_pi;
  real<lower=0,upper=10> tau_days_pi;
  real<lower=0,upper=200> tau_batch_pi;
  
  
  // matrix[ncolY-1,ncolX] beta1; // coefficients (raw)
  vector[ncolY-1] A1; //parameters for the Genotype HOM and Spline interactions
  
  // matrix[ncolY-1,ncolX] beta2; // coefficients (raw)
  vector[ncolY-1] A2; //parameters for the Genotype HET and Spline interactions
  
  // matrix[ncolY-1,ncolX] beta3; // coefficients (raw)
  vector[ncolY-1] A3; //parameters for the Genotype HOM and Spline interactions
  
  vector[ncolY-1] A4; //parameters for the Genotype HET and Spline interactions
  
  vector[nd] id; //mouse id matrix
  vector[ndays] days_id; //days id matrix
  vector[nbatch] batch_id; //batch id matrix
  
  // vector[nd] id2; //mouse id matrix
  // vector[ndays] days_id2; //days id matrix
  // vector[nbatch] batch_id2; //batch id matrix
  
  vector[nd] id3; //mouse id matrix
  vector[ndays] days_id3; //days id matrix
  vector[nbatch] batch_id3; //batch id matrix
  
  row_vector[num_basis] a1_1; // coefficients for bsplines
  row_vector[num_basis] a2_1; // coefficients for bsplines
  row_vector[num_basis] a3_1; // coefficients for bsplines
  row_vector[num_basis] a4_1; // coefficients for bsplines
  row_vector[num_basis] a5_1; // coefficients for bsplines
  row_vector[num_basis] a6_1; // coefficients for bsplines
  row_vector[num_basis] a7_1; // coefficients for bsplines
  row_vector[num_basis] a8_1; // coefficients for bsplines
  
  // row_vector[num_basis] a1_2; // coefficients for bsplines
  // row_vector[num_basis] a2_2; // coefficients for bsplines
  // row_vector[num_basis] a3_2; // coefficients for bsplines
  // row_vector[num_basis] a4_2; // coefficients for bsplines
  // row_vector[num_basis] a5_2; // coefficients for bsplines
  // row_vector[num_basis] a6_2; // coefficients for bsplines
  // row_vector[num_basis] a7_2; // coefficients for bsplines
  // row_vector[num_basis] a8_2; // coefficients for bsplines
  
  row_vector[num_basis] a1_3; // coefficients for bsplines
  row_vector[num_basis] a2_3; // coefficients for bsplines
  row_vector[num_basis] a3_3; // coefficients for bsplines
  row_vector[num_basis] a4_3; // coefficients for bsplines
  row_vector[num_basis] a5_3; // coefficients for bsplines
  row_vector[num_basis] a6_3; // coefficients for bsplines
  row_vector[num_basis] a7_3; // coefficients for bsplines
  row_vector[num_basis] a8_3; // coefficients for bsplines
}

transformed parameters{
  simplex[ncolY] Prob[N]; // unit-simplex of probabilities for the multinomial distribution
  vector<lower=0,upper=1>[ncolY-1] mu_logits[N]; // mean mu following logit transformation
  vector<lower=0,upper=1>[ncolY-1] sig_logits[N]; // dispersion sigma following logit transformation
  vector<lower=0,upper=1>[ncolY-1] pi_logits[N]; // bernoulli parameter following logit transformation
  vector<lower=0>[ncolY-1] alpha[N]; // alpha parameter for beta distribution
  vector<lower=0>[ncolY-1] beta[N]; // beta parameter for beta distribution
  
  vector[num_data] drink_splines; // vector of bsplines times coefficients
  vector[num_data] eat_splines; // vector of bsplines times coefficients
  vector[num_data] ebh_splines; // vector of bsplines times coefficients
  vector[num_data] groom_splines; // vector of bsplines times coefficients
  vector[num_data] hang_splines; // vector of bsplines times coefficients
  vector[num_data] rear_splines; // vector of bsplines times coefficients
  vector[num_data] rest_splines; // vector of bsplines times coefficients
  vector[num_data] sniff_splines; // vector of bsplines times coefficients
  // vector[num_data] walk_splines; // vector of bsplines times coefficients
  vector[ncolY-1] beta_splines[N]; // beta_splines matrix for each behaviour
  
  // vector[num_data] drink_splines2; // vector of bsplines times coefficients
  // vector[num_data] eat_splines2; // vector of bsplines times coefficients
  // vector[num_data] ebh_splines2; // vector of bsplines times coefficients
  // vector[num_data] groom_splines2; // vector of bsplines times coefficients
  // vector[num_data] hang_splines2; // vector of bsplines times coefficients
  // vector[num_data] rear_splines2; // vector of bsplines times coefficients
  // vector[num_data] rest_splines2; // vector of bsplines times coefficients
  // vector[num_data] sniff_splines2; // vector of bsplines times coefficients
  // // vector[num_data] walk_splines2; // vector of bsplines times coefficients
  // vector[ncolY-1] beta_splines2[N]; // beta_splines matrix for each behaviour
  
  vector[num_data] drink_splines3; // vector of bsplines times coefficients
  vector[num_data] eat_splines3; // vector of bsplines times coefficients
  vector[num_data] ebh_splines3; // vector of bsplines times coefficients
  vector[num_data] groom_splines3; // vector of bsplines times coefficients
  vector[num_data] hang_splines3; // vector of bsplines times coefficients
  vector[num_data] rear_splines3; // vector of bsplines times coefficients
  vector[num_data] rest_splines3; // vector of bsplines times coefficients
  vector[num_data] sniff_splines3; // vector of bsplines times coefficients
  // vector[num_data] walk_splines3; // vector of bsplines times coefficients
  vector[ncolY-1] beta_splines3[N]; // beta_splines matrix for each behaviour
  
  for (i in 1:N){ // modelling stick breaking process
    real sumP = 0;
    Prob[i,1] = Zed[i,1];
    sumP = Prob[i,1];
    for (n in 2:(ncolY - 1)) {
      Prob[i,n] = (1 - sumP) * Zed[i,n];
      sumP += Prob[i,n];
    }
    Prob[i,ncolY] = 1 - sumP;
  }
  
  
  drink_splines = to_vector(a1_1*Bsplines);
  eat_splines = to_vector(a2_1*Bsplines);
  ebh_splines = to_vector(a3_1*Bsplines);
  groom_splines = to_vector(a4_1*Bsplines);
  hang_splines = to_vector(a5_1*Bsplines);
  rear_splines = to_vector(a6_1*Bsplines);
  rest_splines = to_vector(a7_1*Bsplines);
  sniff_splines = to_vector(a8_1*Bsplines);
  
  // drink_splines2 = to_vector(a1_2*Bsplines);
  // eat_splines2 = to_vector(a2_2*Bsplines);
  // ebh_splines2 = to_vector(a3_2*Bsplines);
  // groom_splines2 = to_vector(a4_2*Bsplines);
  // hang_splines2 = to_vector(a5_2*Bsplines);
  // rear_splines2 = to_vector(a6_2*Bsplines);
  // rest_splines2 = to_vector(a7_2*Bsplines);
  // sniff_splines2 = to_vector(a8_2*Bsplines);
  
  drink_splines3 = to_vector(a1_3*Bsplines);
  eat_splines3 = to_vector(a2_3*Bsplines);
  ebh_splines3 = to_vector(a3_3*Bsplines);
  groom_splines3 = to_vector(a4_3*Bsplines);
  hang_splines3 = to_vector(a5_3*Bsplines);
  rear_splines3 = to_vector(a6_3*Bsplines);
  rest_splines3 = to_vector(a7_3*Bsplines);
  sniff_splines3 = to_vector(a8_3*Bsplines);

  for (n in 1:N){
    beta_splines[n] = [drink_splines[n], eat_splines[n], ebh_splines[n], groom_splines[n], hang_splines[n], rear_splines[n], rest_splines[n], sniff_splines[n]]';
    // beta_splines2[n] = [drink_splines2[n], eat_splines2[n], ebh_splines2[n], groom_splines2[n], hang_splines2[n], rear_splines2[n], rest_splines2[n], sniff_splines2[n]]';
    beta_splines3[n] = [drink_splines3[n], eat_splines3[n], ebh_splines3[n], groom_splines3[n], hang_splines3[n], rear_splines3[n], rest_splines3[n], sniff_splines3[n]]';
  }
  
  for (n in 1:N) { // linking the mean and dispersion to parameters & design matrix through logit transformation
    for (m in 1:(ncolY-1)){
      mu_logits[n,m] = X[n,] * transpose(params_mu[m,]) + beta_splines[n,m] + A1[m] * X[n,2] * beta_splines[n,m]  + A2[m] * X[n,3] * beta_splines[n,m] + batch_id[batch[n]] + id[mice_id[n]] + days_id[days[n]];
      sig_logits[n,m] = X[n,] * transpose(params_sig[m,]);
      mu_logits[n,m] = exp(mu_logits[n,m]) / (1 + exp(mu_logits[n,m]));
      sig_logits[n,m] = exp(sig_logits[n,m]) / (1 + exp(sig_logits[n,m]));
      pi_logits[n,m] = X[n,] * transpose(params_pi[m,]) + beta_splines3[n,m] + A3[m] * X[n,2] * beta_splines3[n,m] + A4[m] * X[n,3] * beta_splines3[n,m] + batch_id3[batch[n]] + id3[mice_id[n]] + days_id3[days[n]];
      pi_logits[n,m] = exp(pi_logits[n,m]) / (1 + exp(pi_logits[n,m]));
      alpha[n,m] = mu_logits[n,m] * ((1 / sig_logits[n,m]) - 1); // reparameterizing mu and sigma in terms of alpha and beta
      beta[n,m] = (1 - mu_logits[n,m]) * ((1 / sig_logits[n,m]) - 1);
    }
  }
  
}

model {
  // priors
  for (j in 1:ncolX) {
    for (k in 1:(ncolY-1)) {
      params_mu[k,j] ~ normal(0,5);
      params_sig[k,j] ~ normal(0,5);
      params_pi[k,j] ~ normal(0,5);
    }
  }
  // Interaction priors
  for (i in 1:(ncolY-1)){
    A1[i] ~ normal(0, 5);
    A2[i] ~ normal(0, 5);
    A3[i] ~ normal(0, 5);
    A4[i] ~ normal(0, 5);
  }
  // random effect priors
  to_vector(id) ~ normal(0, tau_id_mu);
  to_vector(days_id) ~ normal(0,tau_days_mu);
  to_vector(batch_id) ~ normal(0,tau_batch_mu);
  // random effect priors
  // to_vector(id2) ~ normal(0, 1);
  // to_vector(days_id2) ~ normal(0,1);
  // to_vector(batch_id2) ~ normal(0,1);
  // random effect priors
  to_vector(id3) ~ normal(0, tau_id_pi);
  to_vector(days_id3) ~ normal(0,tau_days_pi);
  to_vector(batch_id3) ~ normal(0,tau_batch_pi);
  // spline priors
  a1_1 ~ normal(0,1);
  a2_1 ~ normal(0,1);
  a3_1 ~ normal(0,1);
  a4_1 ~ normal(0,1);
  a5_1 ~ normal(0,1);
  a6_1 ~ normal(0,1);
  a7_1 ~ normal(0,1);
  a8_1 ~ normal(0,1);
  // spline priors
  // a1_2 ~ normal(0,1);
  // a2_2 ~ normal(0,1);
  // a3_2 ~ normal(0,1);
  // a4_2 ~ normal(0,1);
  // a5_2 ~ normal(0,1);
  // a6_2 ~ normal(0,1);
  // a7_2 ~ normal(0,1);
  // a8_2 ~ normal(0,1);
  // spline priors
  a1_3 ~ normal(0,1);
  a2_3 ~ normal(0,1);
  a3_3 ~ normal(0,1);
  a4_3 ~ normal(0,1);
  a5_3 ~ normal(0,1);
  a6_3 ~ normal(0,1);
  a7_3 ~ normal(0,1);
  a8_3 ~ normal(0,1);
  
  // likelihood
  for (n in 1:N) {
    for (m in 1:(ncolY-1)){
      if (Y[n,m] == 0){
        target += log(pi_logits[n,m]);
      } else {
        target += log(1 - pi_logits[n,m]) + beta_lpdf(Zed[n,m] | alpha[n,m], beta[n,m]);
      }
    }
    target += multinomial_lpmf(Y[n,] | Prob[n,]);
  }

}


generated quantities {
  int<lower=0,upper=3600> yrep[N, ncolY]; // y replicates

  for (n in 1:N) {
    yrep[n,] = multinomial_rng(Prob[n,], 3600);
  }

  
}



