##############packages#################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(rstan)
SEED <- 48927 # set random seed for reproducability
set.seed(SEED)
library(dplyr)

#############data######################################
##################abundance and sampling effort########

load(file = "abundance_species_trait_30.RData")

###################trait data for each species#########

load(file = "traits_162_after_scaling.RData")
#df_cov<-as.data.frame(X_traits[-86,])


######standardize covariates
#df_cov<-X_traits
#df_scaled <- df_cov %>%
 # mutate(across(where(is.numeric), scales::rescale))
#df_scaled
#x_cov<-df_scaled
#df_scaled
#x_cov<-df_scaled[,c(3:14)] 

###getting y and sampling matrix
y<-matrix(0,162,30)  ##times series of each species
sampling<-matrix(0,162,30)  ##sampling effort 

y<-Vector_abun
sampling<-Vector_effort


cov<-x_cov
K <-3
fit_gmm <- function(K) {
  stan(
    file = "cos_effort_covariates.stan",
    data = list(
    N =dim(y)[1],
    P=dim(cov)[2],
    T=15,
    x=1:30,
    freq=1:15/30,
    y=y,
    cov=cov,
    S=sampling, K = K),
    chains = 1, iter = 4000, warmup =2000  )
}

Kmax <- 10
fits <- lapply(3:5, function(K) fit_gmm(K))
names(fits) <- paste0("K", 3:5)



file_name <- paste("fits_new_sin_cos_species_neg_3_cov_family3_5_new",".RData")

save(fits, file = file_name)




