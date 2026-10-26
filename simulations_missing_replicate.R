#install.packages("rstan")
library(rstan)
library(dplyr)
library(tidyr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(rstan)

#install.packages("V8")
#install.packages("V8", repos = "https://cloud.r-project.org/")


x = 1:20

n=c(100,200,300)
 K=3
###loading the amplitude values##########
load("output_cov_bal200.RData")

for (d in 1:3){
  N=n[d]  #number of time series in each simulation
# Initialize the vector of site wise intercepts (~log mean abundance over time)
  C0.vec <- matrix(nrow = n[d],ncol = 1)
  C0.vec <- rep(1.7,n[d])

# Initialize the matrix of site-wise cluster probabilities
# and cluster identities
 
  theta = matrix(nrow=n[d], ncol=K)
  zz = matrix(nrow=n[d], ncol=K)
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

  for(h in 1:100){
    
  for (i in 1:N){
    theta[i,] = rdirichlet(1,rep(1,K))
  # sample cluster identity per site
    zz[i,] = rmultinom(1,1,theta[i,])  
  }
  print(colSums(zz))

# Sample the log(lambda) for each site
  log.lambda = matrix(nrow=N,ncol=20,0)


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


  y_missing<-matrix(0,nrow=100,ncol=20)
  y_missing<-ynewz
  # Select 25 time series
   series_selected <- sample(1:300, 75)

  # For each selected time series
  for (i in series_selected) {
  
      # Random number of missing values: 1 to 10
      n_missing <- sample(1:10, 1)
  
  # Random time points
      time_points <- sample(1:20, n_missing)
  
  # Replace with NA
  y_missing[i, time_points] <- NA
  }
  # print(y_missing)
  df <- as.data.frame(y_missing) %>%
  dplyr::mutate(jj = dplyr::row_number()) %>%
  tidyr::pivot_longer(
    cols = -jj,
    names_to = "kk",
    values_to = "yy"
  ) %>%
  dplyr::mutate(
    kk = readr::parse_number(kk)  
  )
  
  df <- df %>% 
  drop_na()     # removes all rows containing NA
  # Extract vector of all yy values
  yy <- df$yy

  M  <- nrow(df)
  jj <- df$jj
  kk <- df$kk
  yy <- df$yy

  fit_gmm <- function(K) {
    stan(
      file = "sin_cos_sparse_simulations.stan",
      data = list(
        J = nrow(y_missing),
        KK = ncol(y_missing),
        N = dim(ynewz)[1],
        y=y_missing,
        K=K,
        T=10,
        x=1:20,
        freq=1:9/20,
         N=N,
         M  = length(jj),
        jj = as.integer(jj),   # must be integer
        kk = as.integer(kk),
        yy = as.integer(yy),
        T=20,
        x=1:20,
        freq=1:10/20,
        K = K, alpha=rep(0.5,K)),
      chains = 1, iter = 2000, warmup = 1000)
  }
  Kmax <- 3
  fits <- lapply(3, function(K) fit_gmm(K))
  names(fits) <- paste0("K", 3)

  fit_ss <- rstan::extract(fits$K3,permuted=TRUE) 
  filename <- paste0("output_bal_missing_", N + h, ".RData")
  print(filename)
  
  
  save(log.lambda,fit_ss,zz,ynewz,theta, amp_cos_new,amp_sin_new, file =filename)

}

}