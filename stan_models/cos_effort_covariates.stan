data {
  int<lower=0> N;   // number of species -- sites (time series)
  int<lower=2> K;   // number of clusters
  int<lower=0> T;   // half length of the time series
  int<lower=0> P;   // number of covariates
  vector[2*T] x;    // a vector from 1:2*T 
  int  y[N,2*T];    //  observations at site / species 
  vector[T] freq;   // frequency vector
  matrix[N,P] cov;  // covariate value for each site /species
  int  S[N,2*T];    // sampling effort
}

parameters {
      real<lower=-pi(),upper=pi()> phi[K,T]; 
      real<lower=0> ph;
      matrix<lower=0> [K,T] Mu;
      vector[N] Mu0;
      matrix[K,P] B;          // Coefficients for the covariates (for each cluster)
      real alpha[K];        // Intercept for each cluster (for the covariate model)
      matrix<lower=0> [K,N] eps; //random effect
      real<lower=0> sigma_k; //random effect
      real<lower=0> sigma_tau;
      real<lower=0> tau;
      vector<lower=0>[K] gam;
      matrix<lower=0> [K,T] lam;
      real<lower=0> sigma_cov;
      real<lower=0> tau_cov;
      vector<lower=0>[K] gam_cov;

}
transformed parameters{
  simplex[K] theta[N];

  real log_lambda_k[K,2*T];
  for(t in 1:2*T){
    for (k in 1:K){
       log_lambda_k[k,] =  to_array_1d(Mu[k,]*cos((2*pi()*freq*to_row_vector(x)) + transpose(rep_matrix((to_row_vector(phi[k,])),2*T))));
    }
  }

  // Prior for site wise cluster probabilities

  for (i in 1:N) {
    vector[K] logits;
    
    // Compute logits (alpha[k] + dot_product(cov[i,], B[,k]))
    for (k in 1:K) {
      logits[k] = alpha[k] + dot_product(cov[i,], B[k,]) +eps[k,i];
    }
    
    // Apply the softmax to the logits
    theta[i, ] = exp(logits) / sum(exp(logits));
  }
}


model {
  for (k in 1:K){
    alpha[k]~ normal(0,0.5);
  }
 
  sigma_cov ~ student_t(4,0,1);
  tau_cov ~ cauchy(0,sigma_cov);
  for (k in 1:K){
     gam_cov[k] ~ gamma(0.5,1);

     B[k,]~ normal(0, tau_cov * gam_cov[k]);
     
  }
 

  //shrinkage prior for mu
  sigma_tau ~ student_t(4,0,1);
  tau ~ cauchy(0,sigma_tau);
  for (k in 1:K){
    gam[k] ~ gamma(0.5,1);
    for (t in 1:T){
      lam[k,t] ~ inv_gamma(0.5,1);
      Mu[k,t]~ normal(0,lam[k,t]*gam[k]*tau);;
     
    }
  }
  // Prior for site wise intercept
  Mu0 ~ normal (0,2);

  for (k in 1:K){
    sigma_k ~ normal(0, 0.5);
    for (i in 1:N){
      eps[k,i]~ normal(0, sigma_k);
    }
   }
  ph ~ lognormal(0, 1);
  for (k in 1:K){
    for (t in 1:T){
        phi[k,t] ~ normal(0,pi()/2);
    }  
  }

  vector[K] lps = rep_vector(0, K);
  for (i in 1:N) {
    for (k in 1:K){
     lps[k] =log(theta[i,k])+ neg_binomial_2_lpmf(y[i,] | to_row_vector(S[i,]).*(exp(Mu0[i] + to_row_vector(log_lambda_k[k,]))),ph);
     //   lps[k] =log(theta[i,k])+ poisson_lpmf (y[i,] | (exp(Mu0[i] + to_row_vector(log_lambda_k[k,]))));
     }
    target += log_sum_exp(lps);
  }
}

