/// GAUSSIAN PROCESS AGE INDEPENDENT BACK-CALCULATION

// Unlike reference model, vaguely informative prior
// on inverse length-scale

functions{
  // Function to create transition matrices
  matrix tmat_fct(int t, vector q, matrix d, vector ones) {
    matrix[9, 9] Lambda;
    
    Lambda = diag_matrix(ones); // Setting diagonal entries to one
    
    Lambda[1,1] = (1-q[1]) * (1-d[1,t]);
    Lambda[1,2] = q[1] * (1 - d[1,t]);
    Lambda[1,6] = d[1,t];
    
    Lambda[2,2] = (1-q[2]) * (1 - d[2,t]);
    Lambda[2,3] = q[2] * (1 - d[2,t]);
    Lambda[2,7] = d[2,t];
    
    Lambda[3,3] = (1-q[3]) * (1 - d[3,t]);
    Lambda[3,4] = q[3] * (1 - d[3,t]);
    Lambda[3,8] = d[3,t];
    
    Lambda[4,4] = (1-q[4]) * (1 - d[4,t]);
    Lambda[4,5] = q[4] * (1 - d[4,t]);
    Lambda[4,9] = d[4,t];
    
    return Lambda;
  }
}

data{
  int<lower=1> nquar; // number of quarters considered
  int<lower=0> HIV[nquar] ; // data with new HIV diag by year and age
  int<lower=0> AIDS[nquar] ; // data with new AIDS diag by year and age
  int<lower=0> CD4[nquar,4] ; // CD4 data by number of years and ages and CD4 states
  vector<lower=0, upper=1>[4] q; // Fixed progression probabilities
  vector<lower=0>[4] init_prev; // Initial prevalence (as for age-indept model)
  real x[nquar]; // inputs of GP - 1:nquar
  real<lower=0> sigma_sq; // Arbitrary noise of GP
}

transformed data{
  vector[9] ones = rep_vector(1, 9);
  vector[8] zeroes = rep_vector(0, 8);
  vector[nquar] mu = rep_vector(0,nquar);// Prior mean of GP
}

parameters{
  vector[nquar] gamma; // log-infections -- GP
  real<lower=0> eta;
  real<lower=0> rho;
  vector[nquar] delta_raw1; 
  vector[nquar] delta_raw2; 
  vector[nquar] delta_raw3; 
  vector[nquar] delta_raw4; 
  vector<lower=0>[4] vardelta; // variances of diagnoses
  real<lower=0, upper=1> under_rep;
}

transformed parameters{
  matrix[nquar,nquar] Sigma;
  matrix[nquar,nquar] L;
  vector<lower=0>[nquar] h; // number of infs at time t (exp(gamma))
  matrix<lower=0, upper=1>[4, nquar] d; // matrix of diag probs
  matrix<lower=0>[9,nquar] cum_diags ; // expected cum number of individuals in state s at time t 
  row_vector<lower=0>[nquar] incidence_state5; // expected incidence in state 5
  row_vector<lower=0>[nquar] incidence_state6; // expected incidence in state 6
  row_vector<lower=0>[nquar] incidence_state7; // expected incidence in state 7
  row_vector<lower=0>[nquar] incidence_state8; // expected incidence in state 8
  row_vector<lower=0>[nquar] incidence_state9; // expected incidence in state 9
  simplex[4] exp_p[nquar]; // proportion of HIV diagnoses in state k (at time t)
  row_vector<lower=0>[nquar] exp_HIV_dx; // expected number of HIV dx
  row_vector<lower=0>[nquar] exp_AIDS_dx; // expected number of HIV dx
  matrix<lower=0,upper=1>[9,9] PA[nquar]; // collection of model transition matrices
  vector[9] tmpvec; 
  vector[9] tmpvec1; 
  vector[9] init_vec; 
  vector[nquar] delta1;  // random walk diagnoses
  vector[nquar] delta2; 
  vector[nquar] delta3; 
  vector[nquar] delta4; 
  matrix<lower=0>[4,nquar] prev_mat; // expected prevalence
  real<lower=0> inv_rho=inv(rho);

  // INCIDENCE
  // GP model for infections
  
  // COVARIANCE MATRIX OF GP  
  // using cov_exp_quar
  // Note: inv_rho multiplied by two,
  // hence difference (in priors) wrt non cov_exp_quad version
  Sigma = cov_exp_quad(x, eta, inv_rho);
  
  for (i in 1:nquar) {
    Sigma[i,i] = Sigma[i,i] + sigma_sq; // ensures numerical PD
  }
  
  L = cholesky_decompose(Sigma);
  
  // Incidence is exp of log inc-rw
  h = exp(gamma);
  
  // Initializing delta
  delta1[1] = -3.2 + 0.2*delta_raw1[1]; // Diag in 1995 - st1
  delta2[1] = -3.2 + 0.2*delta_raw2[1]; // Diag in 1995 - st2
  delta3[1] = -3.0 + 0.2*delta_raw3[1]; // Diag in 1995 - st3
  delta4[1] = -2.5 + 0.3*delta_raw4[1]; // Diag in 1995 - st4
  
  // Non-centered reparametrization
  // Can probably optimize that 
  // using tail and head
  for(t in 2:nquar){
    delta1[t] = delta1[t-1] + sqrt(vardelta[1])*delta_raw1[t];
    delta2[t] = delta2[t-1] + sqrt(vardelta[2])*delta_raw2[t];
    delta3[t] = delta3[t-1] + sqrt(vardelta[3])*delta_raw3[t];
    delta4[t] = delta4[t-1] + sqrt(vardelta[4])*delta_raw4[t];
  }

  // DIAGNOSES 
  // Post 1984 (as model starts in 1995)
  for(t in 1:nquar){
    d[1,t] = inv_logit(delta1[t]);
    d[2,t] = inv_logit(delta2[t]);
    d[3,t] = inv_logit(delta3[t]);
    d[4,t] = inv_logit(delta4[t]);
  }
  
  // EPIDEMIC EVOLUTION
  
  // Creating transition matrices
  for(t in 1:nquar){
    PA[t,,] = tmat_fct(t, q, d, ones);
  }
  
  // Initialization at time 1
  init_vec = append_row(h[1],zeroes); // New infections at time 1
  tmpvec1 = append_row(init_prev,zeroes[1:5]); // 0 initial prevalence in dx states... augmenting vec
  cum_diags[,1] = PA[1,,]' * tmpvec1 + init_vec; // Initial prevalence advancing
  
  // Times 2, ..., nquar
  for(t in 2:nquar){
    tmpvec = append_row(h[t],zeroes);
    tmpvec1 = cum_diags[,t-1];
    cum_diags[,t] = PA[t,,]' * tmpvec1 + tmpvec;
  }
  
  // Incidence
  // Difference in cumulative diagnoses
  
  // No incidence at t=1
  incidence_state5[1] = cum_diags[5,1];
  incidence_state6[1] = cum_diags[6,1];
  incidence_state7[1] = cum_diags[7,1];
  incidence_state8[1] = cum_diags[8,1];
  incidence_state9[1] = cum_diags[9,1];
  
  // t = 2, ..., nquar
  incidence_state5[2:nquar] = tail(cum_diags[5,], nquar - 1) - head(cum_diags[5,], nquar - 1);
  incidence_state6[2:nquar] = tail(cum_diags[6,], nquar - 1) - head(cum_diags[6,], nquar - 1);
  incidence_state7[2:nquar] = tail(cum_diags[7,], nquar - 1) - head(cum_diags[7,], nquar - 1);
  incidence_state8[2:nquar] = tail(cum_diags[8,], nquar - 1) - head(cum_diags[8,], nquar - 1);
  incidence_state9[2:nquar] = tail(cum_diags[9,], nquar - 1) - head(cum_diags[9,], nquar - 1);
  
  // Prevalence matrix
  prev_mat = cum_diags[1:4,];
  
  // Expected number of HIV diagnoses
  // Sum of incidences in states 5 to 8
  exp_HIV_dx = incidence_state6 + incidence_state7 + incidence_state8 + incidence_state9;
   
  // Expected number of AIDS diagnoses
  // Include under-reporting from year 2000 (quarter 20)
  // Avoiding if-else loop
  exp_AIDS_dx[1:20] = incidence_state5[1:20];
  exp_AIDS_dx[21:nquar] = incidence_state5[21:nquar]*under_rep; 
  
  // Including CD4 counts
  // Reliable CD4 data only available from 1991 (quarter=53) onwards
  for(t in 1:nquar){
    exp_p[t,1] = incidence_state6[t] / exp_HIV_dx[t];
    exp_p[t,2] = incidence_state7[t] / exp_HIV_dx[t];
    exp_p[t,3] = incidence_state8[t] / exp_HIV_dx[t];
    exp_p[t,4] = incidence_state9[t] / exp_HIV_dx[t];
  }
  
}

model{
  
  // PRIORS
  //Priors on hyperparameters GP
  rho ~ student_t(4, 0, 1);
  eta ~ normal(4, 1);
  
  // LIKELIHOOD 
  
  // GP for log-infections
  gamma ~ multi_normal_cholesky(mu, L);
  
  // Non centered parametrization for diagnoses
  delta_raw1 ~ normal(0, 1);
  delta_raw2 ~ normal(0, 1);
  delta_raw3 ~ normal(0, 1);
  delta_raw4 ~ normal(0, 1);
  vardelta ~ gamma(1,32); // prior for logit random walk
  
  // Under reporting
  under_rep ~ beta(236, 118);
  
  // Poisson likelihood HIV data
  HIV ~ poisson(exp_HIV_dx);
  
  // Poisson likelihood AIDS data
  AIDS ~ poisson(exp_AIDS_dx);
  
  // Multinomial term CD4
  for(t in 1:nquar){
    CD4[t,] ~ multinomial(exp_p[t,]);
  }
}


generated quantities{
// Log-lik calculations to get WAIC / LOO/ DIC
  vector[nquar] log_lik; 
  for(t in 1:nquar){
    log_lik[t] = poisson_lpmf(HIV[t] | exp_HIV_dx[t]) + poisson_lpmf(AIDS[t] | exp_AIDS_dx[t]) +  multinomial_lpmf(CD4[t,] | exp_p[t,]) ; 
  }
}
  
