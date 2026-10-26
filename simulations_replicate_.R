#install.packages("rstan")
library(rstan)
library(dplyr)
library(tidyr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(rstan)


n=c(100,200,300)

load("output_cov_bal200.RData")


for (d in 1:3){
  N= n[d]
    
    # Initialize the vector of site wise intercepts (~log mean abundance over time)


  C0.vec <- matrix(nrow = n[d],ncol = 1)
  C0.vec <- rep(1.7,n[d])

# Initialize the matrix of site-wise cluster probabilities
# and cluster identities
  K=3
 
#install.packages("LaplacesDemon")
  library(LaplacesDemon)

# Vector of sampling times 
  sampling.times = 1:20

### sample frequency 

  freq_clustered = matrix(nrow=K,ncol=10) 
  for (i in 1:K){
  freq_clustered[i,] = sample(1:10,10,replace=F)
  freq_clustered[i,]=sort(freq_clustered[i,])/20
}

  print(freq_clustered)
# Sample the site-wise cluster identities

  N<-n[d]
  for(h in 1:100){
    
  for (i in 1:N){
    theta[i,] = rdirichlet(1,rep(1,K))
  # sample cluster identity per site
    zz[i,] = rmultinom(1,1,theta[i,])  
  }
  print(colSums(zz))

# Sample the log(lambda) for each site
  log.lambda = matrix(nrow=n,ncol=20,0)

 
  
  ###### Getting log lambda for n sites and 20 time points
  
  for (i in 1:N){
    # pick the frequency vector for a site
    freeq = zz[i,] %*% freq_clustered 
    # Initialize log.lambda to intercept
    log.lambda[i,] = C0.vec[i]
    amp_sin= zz[i,] %*% amp_sin_new 
    amp_cos= zz[i,] %*% amp_cos_new 
    # add the frequency component contributions to log.lambda 
    #amp_cos=
    #amp_sin=
    for (j in 1:dim(freq_clustered)[2]) {
      
      log.lambda[i,] = log.lambda[i,] + amp_sin[j] * sin(2*pi*sampling.times*freeq[j]) +amp_cos[j] * cos(2*pi*sampling.times*freeq[j]) 
      
      
    }
  }
  
  T<-10
  ynewz<-matrix(0,N,20)
  for (i in 1:N){
    # The actual log.lambda and the expected log.lambda
    
    ynewz[i,] = rpois(2*T, exp(log.lambda[i,]))
    
  }


    
  fit_gmm <- function(K) {
    stan(
      file = "sin_cos_sparse_simulations.stan",
      data = list(
        N = dim(ynewz)[1],
        y=ynewz,
        K=K,
        T=10,
        x=1:20,
        freq=1:10/20,
        T=20,
        x=1:20,
        K = K, alpha=rep(0.5,K)),
      chains = 1, iter = 4000, warmup = 2000)
  }
  Kmax <- 3
  fits <- lapply(3, function(K) fit_gmm(K))
  names(fits) <- paste0("K", 3)

  fit_ss <- rstan::extract(fits$K3,permuted=TRUE) 
  filename <- paste0("output_", N + h, ".RData")
  print(filename)
  
  
  save(log.lambda,fit_ss,zz,ynewz,theta, amp_cos_new,amp_sin_new, file =filename)

}



}