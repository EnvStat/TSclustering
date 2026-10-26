#install.packages("rstan")
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(rstan)

#install.packages("V8")
#install.packages("V8", repos = "https://cloud.r-project.org/")
library(LaplacesDemon)
set.seed(1)

load("output_cov_bal200.RData")
K=3
n=c(100,200,300)
for (d in 1:3){
    N=n[d]
    # Initialize the vector of site wise intercepts (~log mean abundance over time)

    C0.vec <- matrix(nrow = n[d],ncol = 1)
    C0.vec <- rep(1.7,n[d])
  ###B matrix for fixed coefficients generated from normal distribution
    B <- matrix(rnorm(K * P), nrow = K, ncol = P)
    alpha<-rnorm(3, mean = 0, sd = 1)
    
### an example of a fixed B matrix
#B <- matrix(c(
#   2.0,  1.2, -1.2,
#  -1.2,  2.0,  1.2,
#   1.2, -1.2,  2.0
#), nrow = 3, byrow = TRUE)
#B_centered <- sweep(B, 2, colMeans(B), "-")
#B<-B_centered


# Initialize the matrix of site-wise cluster probabilities
# and cluster identities
 

    theta = matrix(nrow=n[d], ncol=K)
    zz = matrix(nrow=n[d], ncol=K)



# Vector of sampling times 
    sampling.times = 1:20

### sample frequency 

    freq_clustered = matrix(nrow=K,ncol=10) 
    for (i in 1:K){
      freq_clustered[i,] = sample(1:10,10,replace=F)
      freq_clustered[i,]=sort(freq_clustered[i,])/20
    }



    library(nnet)  # for multinom
    set.seed(1)
    P <- 3
    K <- 3
    for(h in 1:100){

# --- Covariates: no cluster structure by construction ---
    cov <- matrix(rnorm(N * P), nrow = N, ncol = P)
    cov <- scale(cov)
    colnames(cov) <- paste0("X", 1:P)



# --- Logits and class probabilities ---
    logits <- matrix(0, nrow = N, ncol = K)
    for (i in 1:N) {
        for (k in 1:K) {
            logits[i, k] <- alpha[k] + sum(cov[i, ] * B[k, ])
        }
    }          

softmax <- function(eta) {
  exp_eta <- exp(eta)  # subtract max for numerical stability
  exp_eta / sum(exp_eta)
}

theta <- t(apply(logits, 1, softmax))

# --- Two label-generation schemes ---

    # (a) Stochastic draw — what you had, via zz
    zz <- t(apply(theta, 1, function(p) {
    z <- numeric(K)
    z[sample(1:K, size = 1, prob = p)] <- 1
    z
    }))
    true_labels_sampled <- apply(zz, 1, which.max)


# (b) Deterministic argmax — the "noise-free" ground truth
#true_labels_argmax <- apply(theta, 1, which.max)

# Deterministic one-hot encoding
    #labels <- max.col(theta, ties.method = "first")


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
      file = "cos_simulations_cov_poisson2.stan",
      data = list(
        N = dim(ynewz)[1],
        y=ynewz,
        K=K,
        P=3,
        T=10,
        x=1:20,
        cov= cov,
        freq=1:10/20,
        K = K),
      chains = 1, iter = 4000, warmup = 2000)
  }
  Kmax <- 3
  fits <- lapply(3, function(K) fit_gmm(K))
  names(fits) <- paste0("K", 3)
  fit_ss <- rstan::extract(fits$K3,permuted=TRUE) 
  filename <- paste0("output_cov_", N + h, ".RData")

  save(log.lambda,fit_ss,zz,ynewz,theta, B, cov, alpha, amp_cos_new,amp_sin_new, file =filename)

}


}