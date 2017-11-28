// TENSOR PRODUCT FOR REAL DATA

// SAME AS FOR SIMULATION STUDY (but under-reporting)
// YEARLY SCALE
// QUARTER OF YEARLY INFECTIONS ASSUMED TO OCCUR AT THE BEGINNING OF EACH QUARTER
// AGE-DEPENDENT DIAGNOSES
// UNDER-REPORTING

// AGE DEPENDENT DIAGNOSES 
// In usual age-dependent setting there is an
// age specific intercept for each age group
// So the age specific intercept affects all diag rw
// (state 1, state2, state 3 and 4)

functions{
  vector row_sums(matrix X) {
    vector[rows(X)] s ;  
    s =X * rep_vector(1, cols(X)) ;
    return s ;
  }
  
  row_vector col_sums(matrix X) {
    row_vector[cols(X)] s ;
    s =rep_row_vector(1, rows(X)) * X ;
    return s ;
  }
  
  // Q_fct creates quarterly transition (only, no diag) matrices 
  // Assume diagnosis occurs first and if not then progression
  matrix Q_fct(matrix d, matrix q, int t, int a){
    matrix[4,4] Q;
    Q[1,1] = (1 - d[1,t]) * (1 - q[1,a]) ;
    Q[1,2] = (1 - d[1,t]) * q[1,a] ;
    Q[1,3] = 0 ;
    Q[1,4] = 0 ;
    Q[2,1] = 0 ;
    Q[2,2] = (1 - d[2,t]) * (1 - q[2,a]) ;
    Q[2,3] = (1 - d[2,t]) * q[2,a] ;
    Q[2,4] = 0 ;
    Q[3,1] = 0 ;
    Q[3,2] = 0 ;
    Q[3,3] = (1 - d[3,t]) * (1 - q[3,a]) ;
    Q[3,4] = (1 - d[3,t]) * q[3,a] ;
    Q[4,1] = 0 ;
    Q[4,2] = 0 ;
    Q[4,3] = 0 ;
    Q[4,4] = (1 - d[4,t]) * (1 - q[4,a]) ;
    return Q;
  }
 
  
  // P_fct creates the progression part of a transition matrix
  // Takes as input a Q matrix defined by fct above
  // and diagnoses and progression matrices d and q
  // and times and ages t and a
  // P_fct creates the progression part of a transition matrix
  // Takes as input a Q matrix defined by fct above
  // and diagnoses and progression matrices d and q
  // and times and ages t and a
  matrix P_fct(matrix Q, matrix d, matrix q, int t, int a, int l, matrix I4, matrix zero_mat){
    matrix[4,4] Pl;
    matrix[4,5] PA;
    
    if(l==0){
      Pl = I4 + Q + (Q * Q) + (Q * Q * Q);
    }
    else if(l==1){
      Pl = I4 + Q + (Q * Q);
    }
    else if(l==2){
      Pl = I4 + Q ;
    }
    else if(l==3){
      Pl = I4 ;
    }
    else{
      Pl = zero_mat ;
    }
    
    for (i in 1:4) {
      for (k in 1:4) {
        // HIV diagnosis (k=1:4)
        PA[i,k] = Pl[i,k] * d[k,t] ;
      }
      // Progression to AIDS (k=5)
      PA[i,5] = Pl[i,4] * (1 - d[4,t]) * q[4,a] ;
    }
    return PA;
  }

  
  // Links age to age group
  // age 1-10 -> ag1
  // age 11-10 -> ag2
  // age 21-30 -> ag3
  // age 31-52 -> ag4
  int age_to_agroup(int a){
    int ag;
    if(a <= 10){
      ag = 1;
    }
    else if(a > 10 && a <= 20){
      ag = 2;
    }
    else if(a > 20 && a <= 30){
      ag = 3;
    }
    else {
      ag = 4;
    }
    return ag;
  }
  
}


data{
  int<lower=1> nyr; // number of years 
  int<lower=1> nage; // number of ages 
  int<lower=1> ninfpars; // number of infection parameters 
  int<lower=0> HIV[nyr,nage] ; // data with new HIV diag by year and age
  int<lower=0> AIDS[nyr,nage] ; // data with new AIDS diag by year and age
  int<lower=0> CD4[nyr,nage,4] ; // CD4 data by number of years and ages and CD4 states(3)
  matrix[4,nage] init_prev ; // Initial prevalence by age at inf and CD4 states(3) 
  matrix[4,nage] q;// matrix of progressions for different ages.. fixed
  matrix[nyr*nage,ninfpars] X; // design matrix ... 50 is number of inf params
  matrix[ninfpars-1,2*ninfpars-2] S1; 
}

transformed data{
  vector[3] empty_vec = rep_vector(0,3);
  vector[5] empty_vec1 = rep_vector(0,5);
  vector[4] empty_vec2 = rep_vector(0,4);
  matrix[nyr,5] empty_mat = rep_matrix(0,nyr,5);
  matrix[4,4] zero_mat  = rep_matrix(0, 4, 4);
  matrix[4,4] I4 = diag_matrix(rep_vector(1.0, 4));
  int nagroup = 4; // number of age groups
  vector[ninfpars-1] zero = rep_vector(0,ninfpars-1);
}

parameters{
  vector[ninfpars] beta; //b are the infection parameters
  vector<lower=0>[2] sigma; // 1/lambda ... i.e. sd not precision (stan likes that more) 
  vector<lower=0>[4] vardelta;
  vector[nyr] delta_raw1; 
  vector[nyr] delta_raw2; 
  vector[nyr] delta_raw3; 
  vector[nyr] delta_raw4; 
  vector[nagroup] alpha; // age group specific (ag) and state (k) specific intercept for dx probs --- rule to map to vector [(ag-1)*4 + k - 1]
  real<lower=0, upper=1> under_rep;
}


transformed parameters{
  // Here define all parameters which are used in model block
  // All other parameters are defined locally (efficiency)
  matrix[nyr,nage] exp_HIV; // expected number of HIV diag at time t age a
  matrix[nyr,nage] exp_AIDS; // expected number of HIV diag at time t age a
  simplex[4] p_tak[nyr,nage]; // proportion of HIV diagnoses in state k (time t and age a)
  matrix[ninfpars-1,ninfpars-1] K1_inv; // precision matrix prior for spline coefficient
  matrix[4,nyr] d[4];// matrix of diagnoses for different ages.. now fixed, eventually estimated
  vector[nyr] delta1; 
  vector[nyr] delta2; 
  vector[nyr] delta3; 
  vector[nyr] delta4; 
  vector<lower=0>[2] lambda; // smoothing parameters

  lambda = inv(square(sigma));
  // PRIORS for the tensor product
  K1_inv = S1[1:(ninfpars-1),1:(ninfpars-1)] * lambda[1]  + S1[1:(ninfpars-1),ninfpars:(2*ninfpars-2)] * lambda[2];

  // local parameters block 
  {
    vector[nyr*nage] H; //vector of expected log responses
    matrix[4,4] QA[nyr, nage, nagroup]; // probabilities of transition within latent states given time t and age of inf a0... for already latent individuals
    matrix[4,4] QA1[nyr, nage, nagroup]; // given infections in Q1
    matrix[4,4] QA2[nyr, nage, nagroup]; // given infections in Q2
    matrix[4,4] QA3[nyr, nage, nagroup]; // given infections in Q3
    matrix[4,5] PA[nyr, nage, nagroup]; // probabilities of transition from latent to observed states given time t and age of inf a0... for already latent individuals
    matrix[4,5] PA1[nyr, nage, nagroup]; // given infections in Q1
    matrix[4,5] PA2[nyr, nage, nagroup]; // given infections in Q2
    matrix[4,5] PA3[nyr, nage, nagroup]; // given infections in Q3
    matrix[nyr, 5] diag_arr[nyr,nage]; // expected number of diag at time t age a given infection at time t0 state k
    row_vector[5] exp_diag[nyr,nage]; // expected number of diag at time t age a state k

    // INCIDENCE
    H = exp(X*beta); // These are log infection...
    
   // Initializing delta
   delta1[1] = -3.2 + 0.2*delta_raw1[1]; // Diag in 1995 - st1
   delta2[1] = -3.2 + 0.2*delta_raw2[1]; // Diag in 1995 - st2
   delta3[1] = -3.0 + 0.2*delta_raw3[1]; // Diag in 1995 - st3
   delta4[1] = -2.5 + 0.3*delta_raw4[1]; // Diag in 1995 - st4
    
   // Non-centered reparametrization
   // Can probably optimize that 
   // using tail and head
   for(t in 2:nyr){
     delta1[t] = delta1[t-1] + vardelta[1]*delta_raw1[t];
     delta2[t] = delta2[t-1] + vardelta[2]*delta_raw2[t];
     delta3[t] = delta3[t-1] + vardelta[3]*delta_raw3[t];
     delta4[t] = delta4[t-1] + vardelta[4]*delta_raw4[t];
   }
      
    // diag probs are inv logit of delta
    // so that are normalised btwn 0 and 1
    // inv logit not vectorized -> use loop
    // diag probs are inv logit of delta
    // so that are normalised btwn 0 and 1
    // inv logit not vectorized -> use loop
    for(ag in 1:nagroup){
      for(i in 1:nyr){
        d[ag,1,i] = inv_logit(delta1[i]+alpha[ag]);
        d[ag,2,i] = inv_logit(delta2[i]+alpha[ag]);
        d[ag,3,i] = inv_logit(delta3[i]+alpha[ag]);
        d[ag,4,i] = inv_logit(delta4[i]+alpha[ag]);
      }
    }
    
    // Looping over t and a to get trans mat
    // for all ages and times
    for (t in 1:nyr) {
      for (a in 1:nage) {
        for(ag in 1:nagroup){
          matrix[4,4] Q = Q_fct(d[ag,,], q, t, a);
          // Quarterly probabilities of progression between latent HIV states
          // Q is the transition matrix for progressions (age and time specific)

          //Annual progression probabilities between latent HIV states
          QA[t,a,ag,,] = Q * Q * Q * Q;
          QA1[t,a,ag,,] = Q * Q * Q; // If infection in Q1, no moves in first quarter and then three quarters of transition 
          QA2[t,a,ag,,] = Q * Q; // If infection in Q2, no moves in quarter two and then two quarters of transition 
          QA3[t,a,ag,,] = Q; // If infection in Q3, no moves in quarter three and then one quarter of transition ... no prog in Q4 new infs only enter the model
        
          //Annual progression probabilities to observed states
          PA[t,a,ag,,] = P_fct(Q, d[ag,,], q, t, a, 0, I4, zero_mat); // If diagnosis occur at some (unspecified time in year), if inf occured in previous year
          PA1[t,a,ag,,] = P_fct(Q, d[ag,,], q, t, a, 1, I4, zero_mat); // If inf occurs in Q1, and diagnosed within three quarters
          PA2[t,a,ag,,] = P_fct(Q, d[ag,,], q, t, a, 2, I4, zero_mat); // If inf occurs in Q2, and diagnosed within two quarters
          PA3[t,a,ag,,] = P_fct(Q, d[ag,,], q, t, a, 3, I4, zero_mat); // If inf occurs in Q3, and diagnosed within one quarter... no dx in Q4 as can only enter model
        }
        // Also Setting all entries of diag_arr to 0
        // (default: no new arrival)
        // Needed: if not get nan
        diag_arr[t,a] = empty_mat;
      }
    }

    for(t0 in 1:nyr){
      for(a0 in 1:nage){
        // Local variables
        int ag0;
        vector[4] start;
        vector[4] lat_p; // Initial prevalence ... still latent after time 1
        vector[5] diag_p; // Initial prevalence ... diagnosed after time 1
        vector[4] lat_arr[nyr]; // expected number of latent individuals at certain time
  
        // At time = time of inf
        // New infs enter model in state 1
        
        // At time 1 we need to consider initial prevalence
        // Prevalence is latent in 1994Q4 (t = 0)
        // So consider diagnosis and progressions
        lat_p = empty_vec2;
        diag_p = empty_vec1;
        ag0 = age_to_agroup(a0);

        if(t0 == 1){
            lat_p = (QA[1,a0,ag0,,]' * init_prev[,a0]);
            diag_p =  (PA[1,a0,ag0,,]' * init_prev[,a0]);
        }
        
        // Infections occur in quarters
        // So a quarter of infections enters model
        // At the beginning of each quarter
        start = append_row( H[(a0-1)*nyr+t0]/4, empty_vec);
  
        lat_arr[t0,] = (QA1[t0,a0,ag0,,]' * start + QA2[t0,a0,ag0,,]' * start + QA3[t0,a0,ag0,,]' * start + start) + lat_p;
        diag_arr[t0,a0,t0,] = (PA1[t0,a0,ag0,,]' * start + PA2[t0,a0,ag0,,]' * start + PA3[t0,a0,ag0,,]' * start)' + diag_p'; 
        
        for(t in (t0+1):nyr){ // might worry if t = nyr, but in STAN loops in (nyr+1):nyr are not evaluated (unlike R, like JAGS) --- saving an if loop
          int a = min(a0+t-t0,nage); // age a given age a0 and time t0 of in and current time t
          int ag =  age_to_agroup(a);
          vector[4] temp = lat_arr[t-1,]; // expected number of diag at time t age a state k

          // Here bringing forward in time epidemics
          // of infected individuals... so move on yearly steps (QA, PA)
          // No new infections enter the model
          lat_arr[t,] = (QA[t,a0,ag,,]' * temp);
          if(a == nage){
            row_vector[5] x;
            x = diag_arr[t,a,t0,];
            diag_arr[t,a,t0,] = x + (PA[t,a0,ag,,]' * temp)';
          }else{
            diag_arr[t,a,t0,] = (PA[t,a0,ag,,]' * temp)';
          }
        }
      }
    }
    

    
    for(t in 1:nyr){
      for(a in 1:nage){
        exp_diag[t,a,] = col_sums(diag_arr[t,a,,]);
        // expected HIV diagnoses
        // expected AIDS diagnoses
        exp_HIV[t,a] = exp_diag[t,a,1] + exp_diag[t,a,2] + exp_diag[t,a,3] + exp_diag[t,a,4]; 
        
        if(t <= 5){ ## Under-reporting from 2000 onwards (t=6)
          exp_AIDS[t,a] = exp_diag[t,a,5]; 
        }else{
          exp_AIDS[t,a] = exp_diag[t,a,5]*under_rep; 
        }
      
        // Now we want the proportion of CD4 that
        // is diangosed with HIV in state k at time t and age a (p_kta)
        p_tak[t,a,1] = exp_diag[t,a,1] / exp_HIV[t,a];
        p_tak[t,a,2] = exp_diag[t,a,2] / exp_HIV[t,a];
        p_tak[t,a,3] = exp_diag[t,a,3] / exp_HIV[t,a];
        p_tak[t,a,4] = exp_diag[t,a,4] / exp_HIV[t,a];
      }
    }
  }
  
}


model{
  // Priors
  // Splines
  sigma ~ student_t(4,0,200); // Weakly informative priors.. favours simple models 

  beta[1] ~ normal(0,30) ;// 30 s.d corresponding to precision given in JAGS
  beta[2:ninfpars] ~ multi_normal_prec(zero,K1_inv) ;
  
  vardelta ~ gamma(1,32); // prior for logit random walk

  // Under reporting
  under_rep ~ beta(236, 118);

  alpha ~ normal(0,1); // Priors for age specific intercepts

  // Non centered parametrization for diagnoses
  delta_raw1 ~ normal(0, 1);
  delta_raw2 ~ normal(0, 1);
  delta_raw3 ~ normal(0, 1);
  delta_raw4 ~ normal(0, 1);

  // Likelihood of model
  for(t in 1:nyr){
     HIV[t,1] ~ poisson(exp_HIV[t,1]); 
     for(a in 2:nage){ 
        // a is in 2 ages, because for age = 1 -> at some t get Po(0)
        // Eg if a=1, t=5 -> t0=5
        // then cannot get to AIDS cos need 4 transitions (nor HIV < 200)
        // Above only applies to CD4 and AIDS
        AIDS[t,a] ~ poisson(exp_AIDS[t,a]);
        CD4[t,a,] ~ multinomial(p_tak[t,a,]); 
        HIV[t,a] ~ poisson(exp_HIV[t,a]); 
        
      }
  }
}

generated quantities{
  // Log-lik calculations to get WAIC / LOO/ DIC
  vector[nyr*nage] log_lik; // has to be vector to be consistent 
    for(t in 1:nyr){
      log_lik[(1-1)*nyr+t] = poisson_lpmf(HIV[t,1] | exp_HIV[t,1]); 
      for(a in 2:nage){
          log_lik[(a-1)*nyr+t] = poisson_lpmf(HIV[t,a] | exp_HIV[t,a]) + poisson_lpmf(AIDS[t,a] | exp_AIDS[t,a]) +  multinomial_lpmf(CD4[t,a,] | p_tak[t,a,]) ; 
      }
    }
}
