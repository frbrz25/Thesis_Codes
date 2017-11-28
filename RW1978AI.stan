/// FRAN'S ATTEMPT TO CONVERT AGE INDEPENDENT BACK-CALCULATION USING STAN 

// NOTE FOR EQUIVALENCE WITH JAGS MODEL
// JAGS parametrizes normal in terms of precision = 1 / var
// And in Paul's paper priors were on variances
// STAN parametrizes normal in terms of sd ( = sqrt(var) ) 

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
  int<lower=0> CD4[nquar-52,4] ; // CD4 data by number of years and ages and CD4 states
  vector<lower=0, upper=1>[4] q; // Fixed progression probabilities
}

transformed data{
  vector[9] ones;
  vector[8] zeroes;
  int<lower=1> ndiagpars; // number of years with dx (note dx start in 1985, q28)
  int<lower=1> cd4_l; // number of years with CD4 (note CD4 start in 1991, q52)


  ones = rep_vector(1, 9);
  zeroes = rep_vector(0, 8);
  ndiagpars = nquar - 28; 
  cd4_l = nquar - 52;
}

parameters{
  vector[nquar] gamma; // log-infections (random walk)
  real<lower=0> var1; // log infs rw's var for times t = 1, ... , 36
  real<lower=0> var2; // log infs rw's var for times t = 37, ... , nquar
  real<lower=0, upper=1> c_78; // level of dx from 1978 to 1984
  real<lower=0, upper=1> c_84; // level of dx in 1984
  real<lower=0, upper=1> under_rep;
  vector[ndiagpars] delta_raw1;  // random walk diagnoses
  vector[ndiagpars] delta_raw2; 
  vector[ndiagpars] delta_raw3; 
  vector[ndiagpars] delta_raw4; 
  vector<lower=0>[4] vardelta; // variances of diagnoses
}

transformed parameters{
  vector<lower=0>[nquar] sd_vec; // maps the right sd for log-infs rw to right time
  vector<lower=0>[nquar] h; // number of infs at time t (exp(gamma))
  matrix<lower=0, upper=1>[4, nquar] d; // matrix of diag probs
  matrix<lower=0>[9,nquar] cum_diags ; // expected cum number of individuals in state s at time t 
  row_vector<lower=0>[nquar] incidence_state5; // expected incidence in state 5
  row_vector<lower=0>[nquar] incidence_state6; // expected incidence in state 6
  row_vector<lower=0>[nquar] incidence_state7; // expected incidence in state 7
  row_vector<lower=0>[nquar] incidence_state8; // expected incidence in state 8
  row_vector<lower=0>[nquar] incidence_state9; // expected incidence in state 9
  simplex[4] exp_p[cd4_l]; // proportion of HIV diagnoses in state k (at time t)
  row_vector<lower=0>[nquar] exp_HIV_dx; // expected number of HIV dx
  row_vector<lower=0>[nquar] exp_AIDS_dx; // expected number of HIV dx
  matrix<lower=0,upper=1>[9,9] PA[nquar]; // collection of model transition matrices
  vector[9] tmpvec; 
  vector[9] tmpvec1; 
  matrix<lower=0>[4,nquar] prev_mat; // expected incidence in state 9
  vector[ndiagpars] delta1;  // random walk diagnoses
  vector[ndiagpars] delta2; 
  vector[ndiagpars] delta3; 
  vector[ndiagpars] delta4; 
  
  // INCIDENCE
  
  // Creating vector of s.d.
  sd_vec[1:36] = rep_vector(sqrt(var1),36);
  sd_vec[37:nquar] = rep_vector(sqrt(var2),nquar-36);
  
  // Incidence is exp of log inc-rw
  h = exp(gamma);
  
  // DIAGNOSES 
  
  // Initializing delta
  delta1[1] = -3.086429 + 0.5635968*delta_raw1[1]; 
  delta2[1] = -3.086429 + 0.5635968*delta_raw2[1]; 
  delta3[1] = -3.086429 + 0.5635968*delta_raw3[1]; 
  delta4[1] = -3.086429 + 0.5635968*delta_raw4[1]; 

  // Non-centered reparametrization
  // Can probably optimize that 
  // using tail and head
  for(t in 2:ndiagpars){
    delta1[t] = delta1[t-1] + vardelta[1]*delta_raw1[t];
    delta2[t] = delta2[t-1] + vardelta[2]*delta_raw2[t];
    delta3[t] = delta3[t-1] + vardelta[3]*delta_raw3[t];
    delta4[t] = delta4[t-1] + vardelta[4]*delta_raw4[t];
  }
  
  // Diag model made of 3 cases
  // Pre 1984 diag
  // 1984 diag
  // Post 1984 diag (intro of test)
  // See Birrell's Paper for further details
  
  // Pre 1984
  for(t in 1:24){
    d[1,t] = c_78 * c_84 * inv_logit(delta1[1]);
    d[2,t] = c_78 * c_84 * inv_logit(delta2[1]);
    d[3,t] = c_78 * c_84 * inv_logit(delta3[1]);
    d[4,t] = c_78 * c_84 * inv_logit(delta4[1]);
  }
  
  // 1984
  for(t in 25:28){
    d[1,t] = c_84 * inv_logit(delta1[1]);
    d[2,t] = c_84 * inv_logit(delta2[1]);
    d[3,t] = c_84 * inv_logit(delta3[1]);
    d[4,t] = c_84 * inv_logit(delta4[1]);
  }
  
  // Post 1984
  for(t in 29:nquar){
    d[1,t] = inv_logit(delta1[t-28]);
    d[2,t] = inv_logit(delta2[t-28]);
    d[3,t] = inv_logit(delta3[t-28]);
    d[4,t] = inv_logit(delta4[t-28]);
  }
  
  // EPIDEMIC EVOLUTION
  
  // Creating transition matrices
  for(t in 1:nquar){
    PA[t,,] = tmat_fct(t, q, d, ones);
  }
  
  // Initialization at time 1
  cum_diags[1,1] = h[1];
  cum_diags[2:9,1] = zeroes;
  
  // Times 2, ..., nquar
  for(t in 2:nquar){
    tmpvec[1] = h[t];
    tmpvec[2:9] = zeroes;
    tmpvec1 = cum_diags[,t-1];
    cum_diags[,t] = PA[t,,]' * tmpvec1 + tmpvec;
  }
  
  // Incidence
  // Difference in cumulative diagnoses
  
  // No incidence at t=1
  incidence_state5[1] = 0;
  incidence_state6[1] = 0;
  incidence_state7[1] = 0;
  incidence_state8[1] = 0;
  incidence_state9[1] = 0;
  
  // t = 2, ..., nquar
  incidence_state5[2:nquar] = tail(cum_diags[5,], nquar - 1) - head(cum_diags[5,], nquar - 1);
  incidence_state6[2:nquar] = tail(cum_diags[6,], nquar - 1) - head(cum_diags[6,], nquar - 1);
  incidence_state7[2:nquar] = tail(cum_diags[7,], nquar - 1) - head(cum_diags[7,], nquar - 1);
  incidence_state8[2:nquar] = tail(cum_diags[8,], nquar - 1) - head(cum_diags[8,], nquar - 1);
  incidence_state9[2:nquar] = tail(cum_diags[9,], nquar - 1) - head(cum_diags[9,], nquar - 1);
  
  prev_mat = cum_diags[1:4,];
  
  // Expected number of HIV diagnoses
  // Sum of incidences in states 5 to 8
  exp_HIV_dx = incidence_state6 + incidence_state7 + incidence_state8 + incidence_state9;
   
  // Expected number of AIDS diagnoses
  // Include under-reporting from year 2000 (quarter 89)
  exp_AIDS_dx[1:88] = incidence_state5[1:88];
  exp_AIDS_dx[89:nquar] = incidence_state5[89:nquar] * under_rep;
  
  // Including CD4 counts
  // Reliable CD4 data only available from 1991 (quarter=53) onwards
  for(t in 1:cd4_l){
    exp_p[t,1] = incidence_state6[t+52] / exp_HIV_dx[t+52];
    exp_p[t,2] = incidence_state7[t+52] / exp_HIV_dx[t+52];
    exp_p[t,3] = incidence_state8[t+52] / exp_HIV_dx[t+52];
    exp_p[t,4] = incidence_state9[t+52] / exp_HIV_dx[t+52];
  }
  
}

model{
  
  // PRIORS
  
  // Infection
  gamma[1] ~ normal(0.69315, 0.308642);
  
  var1 ~ gamma(1,32);
  var2 ~ gamma(1,32);
  
  vardelta ~ gamma(1,32);
  
  // Rest
  under_rep ~ beta(236, 118);
  c_78 ~ uniform(0, 1);
  c_84 ~ beta(3.413376, 25.679993);
  
  // LIKELIHOOD 
  
  // Random walk model for infections
  tail(gamma, nquar - 1) ~ normal( head(gamma, nquar - 1), head(sd_vec, nquar - 1) );

  // Random walk for diagnoses
  // NON CENTERED PARAMETRIZATION
  delta_raw1 ~ normal(0, 1);
  delta_raw2 ~ normal(0, 1);
  delta_raw3 ~ normal(0, 1);
  delta_raw4 ~ normal(0, 1);
  
  // Poisson likelihood HIV data
  // HIV dx can only occur from quarter 2
  HIV[2:nquar] ~ poisson(exp_HIV_dx[2:nquar]);
  
  // Poisson likelihood AIDS data
  // AIDS dx can only occur from quarter 5
  AIDS[5:nquar] ~ poisson(exp_AIDS_dx[5:nquar]);
  
  // Multinomial term CD4
  for(t in 1:cd4_l){
    CD4[t,] ~ multinomial(exp_p[t,]);
  }
}
