functions{
  
  int num_matches(vector p, real y) {
  int n = 0;
  for (i in 1:rows(p))
    if (p[i] == y)
      n += 1;
  return n;
}
  
int[] which_equal(vector p, real y) {
  int match_positions[num_matches(p, y)];
  int pos = 1;
  for (i in 1:size(p)) {
    if (p[i] == y) {
      match_positions[pos] = i;
      pos += 1;
    }
  }
  return match_positions;
}
}

data {
  int<lower=0> N;
  int<lower=2> K;
  int<lower=0> T;
  int<lower=0> P;
  vector[2*T] x;
  vector[T] freq;
  int<lower=1> M;
  int<lower=1> J;
  int<lower=1> KK;
  array[M] int<lower=1, upper=J> jj;
  array[M] int<lower=1, upper=KK> kk;
  array[M] int yy;
  array[M] int ss;
  matrix[N,P] cov;

  
}

parameters {
     
      real<lower=0> ph;
      matrix<lower=0> [K,T] Mu;
      matrix<lower=0,upper=1> [K,T] a_;
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
  matrix[K,2*T] log_lambda_k;
  matrix [K,T] a_square;
  matrix[K,T] b_square;
  a_square = a_.* Mu.*Mu;
  b_square = Mu.*Mu - a_square;
  // loop over the cluster
  for (k in 1:K){
    log_lambda_k[k,] = sqrt(a_square[k,])*sin(2*pi()*freq*to_row_vector(x))+ sqrt(b_square[k,])*cos(2*pi()*freq*to_row_vector(x));
    //log_lambda_k[k,] =  Mu[k,]*cos((2*pi()*freq*to_row_vector(x)) + transpose(rep_matrix((phi[k,]),2*T)));
  
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
    for (t in 1:T){
      a_[k,t]~ beta(1, 1);
    }
   }
 
  for (k in 1:K){
    sigma_k ~ normal(0, 0.5);
    for (i in 1:N){
      eps[k,i]~ normal(0, sigma_k);
    }
   }
  ph ~ lognormal(0, 1);

  vector[K] lps = rep_vector(0, K);
  for (i in 1:N) {
    for (k in 1:K){
       //  lps[k] =log(theta[i,k])+ poisson_lpmf(yy[(which_equal(to_vector(jj), i))] | to_vector(ss[(which_equal(to_vector(jj), i))] ).* to_vector(exp(Mu0[i] + to_row_vector(log_lambda_k[k, kk[which_equal(to_vector(jj), i)]]))));
      lps[k] =log(theta[i,k])+ neg_binomial_2_lpmf(yy[(which_equal(to_vector(jj), i))] | to_vector(ss[(which_equal(to_vector(jj), i))] ).* to_vector(exp(Mu0[i] + to_row_vector(log_lambda_k[k, kk[which_equal(to_vector(jj), i)]]))),ph);
     }
    target += log_sum_exp(lps);
  }
}

