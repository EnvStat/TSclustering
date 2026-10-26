#############Getting the accuracy and coverage for parameters###################################
###################################packages#############################################
library(label.switching)
library(clue)
library(dplyr)

#####################################initialization############################

results <- data.frame()
d<-100  ###or 200 or 300 depending on number of time series in the simulations
levels <- c(0.85,0.9,0.95)

#################################coverage for amplitude########################
p_j <- list(
  `1` = c(8),
  `2` = c( 8, 9),
  `3` = c(2, 4)
)

#############p_j represents the places where amplitudes not equal to zero###### 
for (p in 1:3) {
  
  for (j in p_j[[as.character(p)]]) {
    
    covered <- rep(NA,100)
    accuracy100<-numeric(100)
    for (level in levels) {
      
      for (s in 101:200) { ##for d=200, 201 to 300. for d=300, 301 to 400
        load(sprintf("output_cov%d.RData", s))
        n<-dim(fit_ss$theta)[1]
        
        true_labels <- apply(zz, 1, which.max)
        run<-stephens(fit_ss$theta,1e-06,200)
        fit_ss_theta_new<-array(c(fit_ss$theta),dim = c(n,3,d))
        reordered_theta<-permute.mcmc(fit_ss_theta_new, run$permutations)
        result_theta <- array(c(reordered_theta$output), dim = c(n,d,3))
        avgData_theta <- apply(result_theta, c(2,3), mean)

        est_labels <- max.col(avgData_theta)
        
        tab <- table(true_labels, est_labels)
        tab
        
        if (nrow(tab) > ncol(tab)) {
          message("Skipping simulation ", s,
                  ": tab has ", nrow(tab), " rows and ",
                  ncol(tab), " columns.")
          next
        }
        
        
        perm <- solve_LSAP(tab, maximum = TRUE)
        
        est_aligned <- rep(NA, length(est_labels))
        
        for (k in 1:ncol(tab)) {
          est_aligned[est_labels == k] <- which(perm == k)
        }
        
        accuracy100[s-d] <- mean(true_labels == est_aligned)
        
        print(mean(true_labels == est_aligned))
        
        # perm[i] = which estimated cluster matches true cluster i
        
        print(perm)
        
        posterior <- fit_ss$Mu[, perm[p], j]
        alpha_level <- 1 - level
        
        ci <- quantile(
          posterior,
          probs = c(alpha_level / 2,
                    1 - alpha_level / 2)
        )
        
        
        true_alpha <-sqrt(amp_cos_new^2 + amp_sin_new^2)
        true_value <- true_alpha[p, j]
        
       
        if (is.finite(true_value)) {
          covered[s-d] <- 
            true_value >= ci[1] &&
            true_value <= ci[2]
          
        }
      }
      results <- rbind(
        results,
        data.frame(
          p = p,
          j = j,
          credible_level = level,
          coverage = mean(covered, na.rm = TRUE),
          n = sum(!is.na(covered))
        )
      )
      
      print(results$coverage)
      
    }
  }
}

save(accuracy100,results,file="results_cov_accuracy_amp_coverage.RData")


##############################################################################
###############coverage for covariates#########################################
##############################################################################

results2 <- data.frame()

levels <- c(0.85,0.9,0.95)



results2<-data.frame()
for (i in 1:3) {
  
  for (j in 1:3) {
    
    covered <- numeric(100)
    accuracy100<-numeric(100)
    for (level in levels) {
      
      for (s in 101:200) {
        load(sprintf("output_cov%d.RData", s))
        
        n<-dim(fit_ss$theta)[1]
        
        true_labels <- apply(zz, 1, which.max)
        run<-stephens(fit_ss$theta,1e-06,200)
        fit_ss_theta_new<-array(c(fit_ss$theta),dim = c(n,3,d))
        reordered_theta<-permute.mcmc(fit_ss_theta_new, run$permutations)
        result_theta <- array(c(reordered_theta$output), dim = c(n,d,3))
        avgData_theta <- apply(result_theta, c(2,3), mean)

        est_labels <- max.col(avgData_theta)
        
        
        tab <- table(true_labels, est_labels)
        tab
        
        if (nrow(tab) > ncol(tab)) {
          message("Skipping simulation ", s,
                  ": tab has ", nrow(tab), " rows and ",
                  ncol(tab), " columns.")
          next
        }
        
        perm <- solve_LSAP(tab, maximum = TRUE)
        est_aligned <- rep(NA, length(est_labels))
        
        for (k in 1:ncol(tab)) {
          est_aligned[est_labels == k] <- which(perm == k)
        }
        
        accuracy100[s-d] <- mean(true_labels == est_aligned)
        
        print(mean(true_labels == est_aligned))
        
        # perm[i] = which estimated cluster matches true cluster i
        
        print(perm)
        
        
        b<-fit_ss$B
        
        reordered_b<-permute.mcmc(b, run$permutations)
        
        
        ncov<-dim(fit_ss$B)[3]
        ncov
        result_b <- array(c(reordered_b$output), dim = c(2000,K,ncov))
        
        posterior <- result_b[, perm[i], j]
        alpha_level <- 1 - level
        
        ci <- quantile(
          posterior,
          probs = c(alpha_level / 2,
                    1 - alpha_level / 2)
        )
        
        print(ci)
        
        true_value <- B[i, j]
        print(true_value)
        
      
        if (is.finite(true_value)) {
          covered[s-d] <- 
            true_value >= ci[1] &&
            true_value <= ci[2]
          
        }
      }
      results2 <- rbind(
        results2,
        data.frame(
          i = i,
          j = j,
          credible_level = level,
          coverage = mean(covered, na.rm = TRUE),
          n = sum(!is.na(covered))
        )
      )
      
      print(results2$coverage)
      
    }
  }
}

save(results2,results,accuracy100,file="results_cov_accuracy_covergage_amp_cov.RData")

##############another way to get the coverage of covariates#############
K<-3

level   <- 0.95
n_reps  <- 100
covered <- matrix(NA, nrow = n_reps, ncol = (K - 1) * 3)

################once, before the loop#############
ref_to <- function(M, r = 1) M - matrix(M[r, ], nrow(M), ncol(M), byrow = TRUE)

coverage_one <- function(draws, B, level = 0.95, ref = 1) {
  # draws: S x K x ncov array, already aligned to the true classes
  K <- dim(draws)[2]; ncov <- dim(draws)[3]
  a <- 1 - level
  B_true_ref <- ref_to(B, ref)
  
  draws_ref <- array(NA_real_, dim(draws))
  for (m in seq_len(dim(draws)[1])) draws_ref[m, , ] <- ref_to(draws[m, , ], ref)
  
  out <- numeric(0)
  for (i in setdiff(1:K, ref)) {
    for (j in seq_len(ncov)) {
      ci  <- quantile(draws_ref[, i, j], c(a / 2, 1 - a / 2))
      out <- c(out, B_true_ref[i, j] >= ci[1] && B_true_ref[i, j] <= ci[2])
    }
  }
  out                                   # length (K - 1) * ncov
}


# ---- after all replications ----
colMeans(covered, na.rm = TRUE)
mean(covered, na.rm = TRUE)

