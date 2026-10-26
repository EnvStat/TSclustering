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
  vector[2*T] x;
  vector[T] freq;
  int<lower=1> M;
  int<lower=1> J;
  int<lower=1> KK;
  array[M] int<lower=1, upper=J> jj;
  array[M] int<lower=1, upper=KK> kk;
  array[M] int yy;
  vector<lower=0.01>[K] alpha;

}

parameters {
      simplex[K] theta[N];

      real<lower=0> ph;
      matrix<lower=0> [K,T] Mu;
      matrix<lower=0,upper=1> [K,T] a_;
      vector[N] Mu0;
      real<lower=0> sigma_tau;
      real<lower=0> tau;
      vector<lower=0>[K] gam;
      matrix<lower=0> [K,T] lam;
}
transformed parameters{
  matrix[K,2*T] log_lambda_k;
  matrix [K,T-1] a_square;
  matrix[K,T-1] b_square;
  a_square = a_.* Mu.*Mu;
  b_square = Mu.*Mu - a_square;
  // loop over the cluster
  for (k in 1:K){
    log_lambda_k[k,] = sqrt(a_square[k,])*sin(2*pi()*freq*to_row_vector(x))+ sqrt(b_square[k,])*cos(2*pi()*freq*to_row_vector(x));
   
  
  }
  // Prior for site wise cluster probabilities

  //  theta[i,] = to_vector(to_row_vector(c) + exp(cov[i,]*transpose(B)))/sum(M + exp(cov[i,]*transpose(B))));
  
  }


model {
 
 // for (k in 1:K){
   // B[k,] ~ normal(0.5,2);
 // }


  // Prior for site wise cluster probabilities
  for (i in 1:N){
    theta[i,] ~ dirichlet(alpha);
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

  vector[K] lps = rep_vector(0, K);
  for (i in 1:N) {
    for (k in 1:K){
     lps[k] = log(theta[i,k])+ poisson_lpmf (yy[(which_equal(to_vector(jj), i))]|to_vector(exp(Mu0[i] + to_row_vector(log_lambda_k[k, kk[which_equal(to_vector(jj), i)]]))));
     // lps[k] =log(theta[i,k])+ neg_binomial_2_lpmf(yy[(which_equal(to_vector(jj), i))] | to_vector(exp(Mu0[i] + to_row_vector(log_lambda_k[k, kk[which_equal(to_vector(jj), i)]]))),ph);
     }
    target += log_sum_exp(lps);
  }
}

