######################################################################
#######################################Packages#######################
######################################################################
install.packages("Rmisc")
install.packages("label.switching")
install.packages("data.table")
install.packages("dplyr")
install.packages("MASS") 
install.packages("reshape2") 
install.packages("reshape") 
install.packages("ggthemes") 
library(label.switching)
library(Rmisc)
library("data.table")
library(dplyr)
library(ggplot2)

library(ggthemes) # 
library(rstan)
library(MASS) 
library(reshape2) 
library(reshape)
library(gridExtra)

####################Extract the posterior############################
fit_ss <- rstan::extract(fits$K3,permuted=TRUE) 

K<-dim(fit_ss$Mu)[2]## number of clusters 
T<-dim(fit_ss$Mu)[3] ##number of frequencies
num<-dim(fit_ss$Mu)[1]  ##number of iterations after warmup
d<-dim(fit_ss$theta)[2] ##number of curves - sites - species...
dim(fit_ss$log_lambda_k)
###implementing the label switching
run<-stephens(fit_ss$theta,1e-06,200)
run
#####################################################################
###After that, we do the permutation based on the previous results 
###########################we do this for Mu#########################
mcmc<-fit_ss$Mu
dim(mcmc)
reordered_mu<-permute.mcmc(mcmc, run$permutations)
result_mu <- array(c(reordered_mu$output), dim = c(num,K,T))
################for  b effect size of the covariate
a<-fit_ss$alpha
b<-fit_ss$B
reordered_b<-permute.mcmc(b, run$permutations)
ncov<-dim(fit_ss$B)[3]
ncov
result_b <- array(c(reordered_b$output), dim = c(num,K,ncov))

B_matrix_mean<-apply(result_b,c(2,3),mean)
alpha_mean<-apply(fit_ss$alpha,c(2),mean)
####permuted theta
fit_ss_theta_new<-array(c(fit_ss$theta),dim = c(num,K,d))
reordered_theta<-permute.mcmc(fit_ss_theta_new, run$permutations)
result_theta <- array(c(reordered_theta$output), dim = c(num,d,K))
fit_ss_mu0_new<-fit_ss$Mu0

###for phase###############
if (length(dim(fit_ss$lambda))==4){
  mcmc_phi<-fit_ss$phi
  reordered_phi<-mcmc_phi
  for (i in 1:d){
    result<-permute.mcmc(mcmc_phi[,,i,], run$permutations)
    reordered_phi[,,i,]<-result$output
    
  }
  result_phi <- reordered_phi
  
  
  #####for lambda (when there is phase)
  
  mcmc_lambda<-fit_ss$log_lambda_k
  reordered_lambda<-mcmc_lambda
  for (i in 1:d){
    result<-permute.mcmc(mcmc_lambda[,,i,], run$permutations)
    reordered_lambda[,,i,]<-result$output
    
  }
  result_lambda<-reordered_lambda
  ##### Getting mean over iterations of phi
  phi_mean<-apply(result_phi, c(2,3,4), mean)
  phi_mean
  dim(phi_mean)
  freq_comp_cluster_phi<-result_phi
  
  phi<-reshape2::melt(freq_comp_cluster_phi)
  
  
  colnames(phi)<-c("id","K","T","value")
  
  }else{
  
    #####permuted lambda (when there are covariates) #####
    
    fit_ss_lambda_new<-fit_ss$log_lambda_k
    reordered_lambda<-permute.mcmc(fit_ss_lambda_new, run$permutations)
    dim(reordered_lambda)
    result_lambda <- array(c(reordered_lambda$output), dim = c(num,K,30))
    
    dim(result_lambda)
    
    
    
  }
#####################################################################################
#############################prepare mu for plotting#################################
#####################################################################################
freq_comp_cluster_Mu<-result_mu

mu<-reshape2::melt(freq_comp_cluster_Mu)

colnames(mu)<-c("id","K","T","value")

#####getting the probability of belonging to a cluster theta
avgData_theta <- apply(result_theta, c(2,3), mean)
cluster_index<-max.col(avgData_theta)  ###getting the cluster index
cluster_index_with_theta<-cbind(1:d,cluster_index,avgData_theta)
newdata <- cluster_index_with_theta[order(cluster_index),]
thetax <- paste0("theta",c(1:(K)))# 
colnames(newdata)<-c("site_id","cluster_index",thetax)
newdata<-data.frame(newdata)
newdata
colMeans(avgData_theta)
df<-newdata
df %>%
  group_by(cluster_index) %>%
  summarise(
    mean = mean(value),
    sd = sd(value)
  )


library(dplyr)

result <- df %>%
  group_by(cluster_index) %>%
  summarise(
    across(
      c(theta1, theta2, theta3),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

###getting the indices for each cluster
indices_cluster_1= newdata$site_id[which(newdata$cluster_index==1)]

indices_cluster_1
indices_cluster_2= newdata$site_id[which(newdata$cluster_index==2)]
indices_cluster_2

indices_cluster_3= newdata$site_id[which(newdata$cluster_index==3)]
indices_cluster_3

indices_cluster_4= newdata$site_id[which(newdata$cluster_index==4)]

indices_cluster_5= newdata$site_id[which(newdata$cluster_index==5)]

indices_cluster_5

#############################Plotting the time series####################
#########################################################################
#########################################################################

library(dygraphs)
library(xts) # To make the convertion data-frame / xts format
# Format 3: Several variables for each date

plot_abundance<-function(indices,ynew,cluster_number,s){
  data <- data.frame(
    time=seq(from=1, to=2*T, by=1 ), 
    site1=ynew[indices[1],]/s[indices[1],]
  )
  
  
  
  index<-c()
  for (i in 2:length(indices))
  {
    
    data<-cbind(data,ynew[indices[i],]/s[indices[i],])
    index<-cbind(index,paste0("s",indices[i])) 
  }
  index
  data
  colnames(data)<-c("time",paste0("s",indices[1]),index)
  data
  p <- dygraph(data, main=paste0("Cluster",cluster_number), xlab = "Time (years)", ylab = "Raw density of abundance")%>%# dyAxis("y", label="Abundance", valueRange = c(0, 60)) %>%
    dyLegend(show = "never", hideOnMouseOut = FALSE)   
  p
  
  
}



plot_abundance_relative<-function(indices,ynew,cluster_number,s){
  data <- data.frame(
    time=seq(from=1, to=30, by=1 ), 
    site1=ynew[indices[1],]/s[indices[1],]/mean(ynew[indices[1],]/s[indices[1],])
  )
  
  
  index<-c()
  for (i in 2:length(indices))
  {
    
    data<-cbind(data,ynew[indices[i],]/s[indices[i],]/mean(ynew[indices[i],]/s[indices[i],]))
    index<-cbind(index,paste0("s",indices[i])) 
  }
  index
  data
  colnames(data)<-c("time",paste0("s",indices[1]),index)
  data
  p <- dygraph(data, main=paste0("Cluster",cluster_number), xlab = "Time (years)", ylab = "raw relative density in abundance") %>%# dyAxis("y", label="Abundance", valueRange = c(0, 60)) %>%
    dyLegend(show = "never", hideOnMouseOut = FALSE)   
  p
  
  
}

####Plotting the clustered abundance curve in each cluster


plot_abundance(indices_cluster_1,ynewzz,1,s) 
plot_abundance(c(indices_cluster_2),ynewzz,2,s)
plot_abundance(indices_cluster_3,ynewzz,3,s)
plot_abundance(indices_cluster_4,ynewzz,4,s)
plot_abundance(indices_cluster_5,ynewzz,5,s)


plot_abundance_relative(indices_cluster_1,ynewzz,1,s)
plot_abundance_relative(indices_cluster_2,ynewzz,2,s)
plot_abundance_relative(indices_cluster_3,ynewzz,3,s)
plot_abundance_relative(indices_cluster_4,ynewzz,4,s)
plot_abundance_relative(indices_cluster_5,ynewzz,5,s)
plot_abundance_relative(indices_cluster_6,ynewzz,6,s)

######################################################################
##################Different plotting for frequency components#########
######################################################################
library(dplyr)
library(Rmisc)
library(ggplot2)
data<-mu
K
colnames(data)<-c("id","K","freq","value")
data
unique(data$K)

data<-data[,1:4] %>% 
  group_by(K, freq) %>% 
  dplyr::summarise(mean = mean(value),
                   lower =  quantile(value,c(0.025)),
                   upper =  quantile(value,c(0.975)))%>%
  mutate(freq = freq %>% as.factor())

###########frequncies will be##############################
############################################################
freq_new<-round(30/c(1:15),1)

data$freq_new<-rep(freq_new, times =K)


data

####################Amplitude plotting####################
###################Frequnecy profile######################

colnames(data)<-c("K","freq","Mean_posterior","lower","upper","freq_new")
data$K<-factor(data$K)
p1<-ggplot(data, aes(x = factor(freq_new), y= Mean_posterior, fill= K)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Corresponding frequencies for each cluster",
       x= "years per cycle",
       y="Frequency amplitude",
       fill= "Cluster") +
  
  theme_classic() +
  
  ##adding a new layer and specifying the criterion for the facets - or how to "split" the data
  facet_grid(~ K, labeller = as_labeller(function(x) paste("Cluster", x)))+
  scale_fill_manual(values = c(
    "1" = "#F28E2B",  # orange
    "2" = "#72C2B5",  # green → light blue
    "3" = "#9B7BC1"   # purple
  ))+theme(legend.position = "none")+theme(
    strip.text = element_text(
      size = 14,
      face = "bold"),
    axis.title.y = element_text(face = "bold",size=14),axis.title.x = element_text(face = "bold",size=14),
    axis.text.y = element_text(face = "bold",size=12),axis.text.x = element_text(face = "bold",size=12)
  )+scale_x_discrete(
    breaks = c( "2.3", "3","3.8","5","6","7.5","10","15","30"))


p1
########################################################
######################LAMBDA plotting#####################
##########################################################
lambda_log<-result_lambda
loglambda<-reshape2::melt(lambda_log)
loglambda$value<-exp(loglambda$value)
dim(loglambda)
colnames(loglambda)<-c("id","K","time","value")
unique(loglambda$K)
clusters = unique(loglambda$K)

n_per_year <- 2000*K
years <- 1993:2022

# repeat each year 1500 times
loglambda$year <- rep(years, each = n_per_year)

loglambda
##############################################################################
####plotting frequency components with confidence interval and true points####
##############################################################################
mean_cl_quantile <- function(x, q = c(0.05, 0.95), na.rm = TRUE){
  dat <- data.frame(y = mean(x, na.rm = na.rm),
                    ymin = quantile(x, probs = q[1], na.rm = na.rm),
                    ymax = quantile(x, probs = q[2], na.rm = na.rm))
  return(dat)
}
###################################################################
##################################################################

cluster_plots = list()
for(cluster_ in 1:3) {
  cluster_plots[[cluster_]] = ggplot(loglambda %>% filter(K == cluster_), aes(x = year, y = value)) + 
    ylim(-4, 6) + labs(y=expression(bold("Log relative intensity, log "~ lambda))) +
    xlab("year") +
    stat_summary(
      geom = "line",
      fun = mean,
      linewidth = 1.2,
      colour = "black"
    
    ) +
    stat_summary(
      geom = "ribbon",
      fun.data = mean_cl_quantile,
      alpha = 0.3
    ) +
    geom_rangeframe() +
    scale_x_continuous(
      breaks = c(2000, 2010, 2020)
    ) +
    theme_tufte() +
    theme(
      # Remove grid lines and borders
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Make the axes darker/thicker
      axis.line = element_line(colour = "black", linewidth = 0.8),
  
      # Axis titles
      axis.title.x = element_text(
        size = 20,
        face = "bold",
        colour = "black"
      ),
      axis.title.y = element_text(
        size = 16,
        face = "bold",
        colour = "black"
      ),
      
      # Axis values
      axis.text.x = element_text(
        size = 20,
        face="bold",
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 20,
        face="bold",
        colour = "black"
     
      ),
      
      legend.position = "none"
    )
  print(cluster_plots[[cluster_]])
  
  # ggsave(cluster_plots[[cluster_]], file=paste0("lambda_162_5", cluster_,".png"), width = 44.45, height = 27.78, units = "cm", dpi=300)
}

cluster_plots[[2]]

p<-grid.arrange(grobs = cluster_plots, ncol=K)

##########################################################################################
############Plotting the frequency profiles and representative time series together#######
##########################################################################################
library(cowplot)
library(ggplot2)

###########Figure 2 and 4 in the manuscript#########
ggdraw() +
  draw_plot(p1, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(p, x = 0, y = 0, width = 1, height = 0.5)



#############################################################################################
########################Plotting the probability of belonging of each cluster in rectangle#############
##############################################################################################

#########Appendix figures###########
bar_height <- 0.4 
no_highlight_col <- 'black'
average_highlight_col <- 'grey40'
avgData_theta<-as.data.frame(avgData_theta)
###if site use this
avgData_theta<-cbind(newdata,site_49[newdata$site_id])
###if species use this
avgData_theta<-cbind(newdata,species2[newdata$site_id])
colnames(avgData_theta)<-c("site_id","cluster_index","V1","V2","V3","SiteID")
sorted_dat <- avgData_theta %>% 
  mutate(num = 1:d)%>%
  arrange(num)


####for theta plotting 

library(ggplot2)
library(dplyr)

bar_height <- 0.4

no_highlight_col <- "black"
average_highlight_col <- "grey40"

y_axis_text_size <- 3
axis_label_size <- 8
text_size <- 18


# Plot
p <- ggplot(sorted_dat) +
  
  # Cluster 1
  geom_rect(
    aes(
      xmin = 0,
      xmax = V1,
      ymin = num - bar_height,
      ymax = num + bar_height
    ),
    fill = "#F28E2B"
  ) +
  
  # Cluster 2
  geom_rect(
    aes(
      xmin = V1,
      xmax = V1 + V2,
      ymin = num - bar_height,
      ymax = num + bar_height
    ),
    fill = "#72C2B5"
  ) +
  
  # Cluster 3
  geom_rect(
    aes(
      xmin = V1 + V2,
      xmax = 1,
      ymin = num - bar_height,
      ymax = num + bar_height
    ),
    fill = "#9B7BC1"
  ) +
  
  # Site labels
  geom_text(
    aes(
      x = -0.02,
      y = num,
      label = SiteID
    ),
    hjust = 1,
    size = y_axis_text_size,
    color = no_highlight_col,
    fontface = "bold"
  ) +
  
  labs(
    x = "Probability",
    y = NULL
  ) +
  
  scale_x_continuous(
    limits = c(-0.3, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  
  coord_cartesian(
    ylim = c(0, 162),  ### you can use here (0,62) and then (62,124) and finally (123,162) to have three separated figures 
    clip = "off"
  ) +
  
  theme_minimal() +
  
  theme(
    text = element_text(
      size = text_size,
      color = average_highlight_col
    ),
    panel.grid = element_blank(),
    
    axis.title.x = element_text(
      size = axis_label_size,
      color = no_highlight_col
    ),
    
    axis.text.x = element_text(
      size = axis_label_size,
      color = no_highlight_col
    ),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    plot.margin = margin(10, 10, 10, 100)
  )

p


#########################################################################################
##############################Response to categorical variable and effect size##########
#########################################################################################




######standardize covariates

df_cov<-as.data.frame(X_traits[-86,])
######standardize covariates
df_scaled <- df_cov %>%
  mutate(across(where(is.numeric), scales::rescale))


x_cov<-df_scaled[,c(3:14)]  
print(x_cov)

#####################plotting overwintering stage vs cluster##################

df_category<- data.frame(
  cat_A = x_cov[,5],
  cat_B = x_cov[,6],
  cat_C = x_cov[,7],
  cat_D = x_cov[,8]
)


label_map2 <- c(
  A = "Adult",
  B = "Caterpillar",
  C = "Egg",
  D="Pupa"
)
#####################use this if you are plotting diet vs cluster##################

df_category<- data.frame(
  cat_A = x_cov[,9],
  cat_B = x_cov[,10],
  cat_C = x_cov[,11],
  cat_D = x_cov[,12]
)

label_map2 <- c(
  A = "Monophagous",
  B = "Mushroom",
  C = "Oligophagous",
  D="Polyphagous"
)


####################################################################################

cat_cols <- df_category %>% dplyr::select(starts_with("cat_"))
cat_cols

df_final <- df_category %>%
  mutate(
    category = colnames(cat_cols)[max.col(cat_cols, ties.method = "first")],
    category = ifelse(rowSums(cat_cols) == 0, "Other", category),
    category = gsub("cat_", "", category),
    category = recode(category, !!!label_map2))


df_final
cat_cols2 <- df_category %>%
  dplyr::select(starts_with("cat_"))

df_final <- df_category %>%
  mutate(
    another_category = colnames(cat_cols2)[
      max.col(as.matrix(cat_cols2), ties.method = "first")
    ],
    another_category = gsub("cat_", "", another_category),
    another_category = recode(another_category, !!!label_map2)
  )

df_final

##############################Response along the covariates################################
####################Response along the covariates###################################
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(tidyr)
library(tidyverse)
library(scales)
###########################For species case############################
######standardize covariates############################################

df_cov<-as.data.frame(X_traits[-86,])
######standardize covariates
df_scaled <- df_cov %>%
  mutate(across(where(is.numeric), scales::rescale))


x_cov<-df_scaled[,c(3:14)]  
print(x_cov)

#####################plotting overwintering stage vs cluster##################

df_category<- data.frame(
  cat_A = x_cov[,5],
  cat_B = x_cov[,6],
  cat_C = x_cov[,7],
  cat_D = x_cov[,8]
)

#####################use this if you are plotting diet vs cluster##################

df_category<- data.frame(
  cat_A = x_cov[,9],
  cat_B = x_cov[,10],
  cat_C = x_cov[,11],
  cat_D = x_cov[,12]
)

df_final


label_map2 <- c(
  A = "Adult",
  B = "Caterpillar",
  C = "Egg",
  D="Pupa"
)
label_map2 <- c(
  A = "Monophagous",
  B = "Mushroom",
  C = "Oligophagous",
  D = "Polyphagous"
)

cat_cols <- df_category %>%
  dplyr::select(starts_with("cat_"))

df_final <- df_category %>%
  mutate(
    category = colnames(cat_cols)[max.col(as.matrix(cat_cols), ties.method = "first")],
    category = ifelse(rowSums(cat_cols) == 0, "Other", category),
    category = gsub("cat_", "", category),
    category = recode(category, !!!label_map2)
  )
###for overwintering 

df<-df_final
head(df)
newdata$overwintering<-c(df$category[indices_cluster_1],df$category[indices_cluster_2],df$category[indices_cluster_3])
################# for diet
newdata$diet<-c(df$category[indices_cluster_1],df$category[indices_cluster_2],df$category[indices_cluster_3])

df<-newdata
newdata

# Reshape data to long format
df_long <- df %>%
  pivot_longer(cols = starts_with("theta"), 
               names_to = "theta_type", 
               values_to = "theta_value")


filtered_data<-df_long


label_map <- c(
  "theta1" = "Cluster 1",
  "theta2" = "Cluster 2",
  "theta3" = "Cluster 3")

filtered_data

p2<-ggplot(filtered_data, aes(x = overwintering, y = theta_value, fill = overwintering)) +
  
  # Violin with transparency
  geom_violin(trim = FALSE, alpha = 0.6, color = NA,drop=FALSE) +
  
  # Add jittered points for raw data
  geom_quasirandom(width = 0.15, size = 1, alpha = 0.7, color = "black") +
  
  # Overlay mean + sd
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, color = "black") +
  
  # Facets
  facet_wrap(~theta_type, scales = "free_y",
             labeller = labeller(theta_type = label_map)) +
  
  # Colors: palette or custom
  scale_fill_brewer(palette = "Set2") +
  
  # Theme enhancements
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1,size=12,face="bold"),
    axis.text.y = element_text(size = 10,face="bold"),
    axis.title = element_text(size =  12,face="bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) + scale_x_discrete(labels = label_map2)+
  
  # Labels
  labs(
    # title = "Distribution of Theta Values by overwintering stage",
    x = "Overwintering Stage",
    y = expression("Probability of belonging to a group"~ theta)
  )


p2

########################For diet###########################
newdata$diet<-c(df$category[indices_cluster_1],df$category[indices_cluster_2],df$category[indices_cluster_3])

df<-newdata
# Reshape data to long format
df_long <- df %>%
  pivot_longer(cols = starts_with("theta"), 
               names_to = "theta_type", 
               values_to = "theta_value")


filtered_data<-df_long


label_map <- c(
  "theta1" = "Cluster 1",
  "theta2" = "Cluster 2",
  "theta3" = "Cluster 3")

p3<-ggplot(df_long, aes(x = diet, y = theta_value, fill = diet)) +
  geom_boxplot(alpha = 1) +  # Violin plot with semi-transparency
  scale_fill_manual(values = c(
    "A"         = "#999999",  # Coral
    "B"       = "#6B5B95",  # Purple
    "C"      = "#88B04B",  # Olive green
    "D"  = "#45B8AC" , # Teal,
    "Other" = "#F2B134"  # Warm Gold
    
  )) +scale_x_discrete(labels = label_map2)+
  facet_wrap(~theta_type, scales = "free_y",
             labeller = labeller(theta_type = label_map)) +  # Faceting based on theta_type
  theme_minimal() +
  labs(
    # title = "Violin Plot of Theta Values by Overwintering Strategy",
    x = "Larval breadth diet within host type",
    y = expression("Probability of belonging to a group "~ theta)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1,size=12,face="bold"),
    axis.text.y = element_text(size = 9, face="bold"),
    axis.title = element_text(size = 12,face="bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) 
# Place the legend at the top

p3
############################### Response to numerical covariates ######################
n<-d
K <- K# number of clusters
ncov<-dim(fit_ss$B)[3]
covariates_to_plot <-1:ncov 

plot_list <- list()
curve_list <- list()

x_cov<- as.matrix(x_cov)
ncov<-length(covariates_to_plot)
for (covar in covariates_to_plot) {
  
  x_temp <- matrix(0, n, ncov)
  
  # Set all covariates at mean
  for (j in 1:ncov) {
    x_temp[, j] <- mean(x_cov[, j])
  }
  
  # Vary only this covariate
  x_temp[, covar] <- seq(min(x_cov[,covar]), max(x_cov[,covar]), length.out = n)
  
  # Compute theta
  for (i in 1:n) {
    theta_temp <- matrix(NA, n, K)
    
    for (i in 1:n) {
      eta <- as.numeric(alpha_mean + x_temp[i,] %*% t(B_matrix_mean))
      theta_temp[i, ] <- exp(eta) / sum(exp(eta))
    }
  }
  
  df <- as.data.frame(theta_temp)
  df$X <- x_temp[, covar]
  df$Covariate <- paste0("Cov", covar)
  
  curve_list[[length(curve_list) + 1]] <- df
}

data_wide <- bind_rows(curve_list)
data_wide[, c("V1", "V2", "V3")]

data_long <- bind_rows(curve_list) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Cluster",
    values_to = "Response"
  ) %>%
  mutate(Cluster = factor(Cluster))

# Combine all covariates
data_long <- bind_rows(curve_list) %>%
  pivot_longer(cols = starts_with("V"), names_to = "Cluster", values_to = "Response") %>%
  mutate(Cluster = factor(Cluster))

# Plot
ggplot(data_long, aes(x = X, y = Response, group = Cluster)) +
  geom_line(aes(group = interaction(Covariate, Cluster)), alpha = 0.7) +
  facet_grid(Cluster ~ Covariate, scales = "free_x") +
  theme_minimal(base_size = 13) +
  labs(y = "Δ Response", x = "Covariate Value") +
  theme(strip.text = element_text(size = 12, face = "bold"))


covariates_to_plot <-1:ncov 
plot_list <- list()
curve_list <- list()

ncov<-length(covariates_to_plot)

curve_list <- list()
n_draws<-dim(fit_ss$alpha)[1]
alpha_samples<-fit_ss$alpha
B_samples <-result_b
r<-dim(result_b)[1]
ncov<-dim(x_cov)[2]
for (covar in covariates_to_plot) {
  # print(covar)
  
  x_temp <- matrix(0, n, ncov)
  
  # Set all covariates at mean
  for (j in 1:ncov) {
    x_temp[, j] <- mean(x_cov[, j])
  }
  
  # Vary only this covariate
  x_temp[, covar] <- seq(min(x_cov[,covar]), max(x_cov[,covar]), length.out = n)
  
  # Placeholder for storing theta from multiple draws
  theta_all_draws <- array(NA, dim = c(n, K, r))  # e.g., n_draws = 1000 posterior draws
  
  for (draw in 1:n_draws) {
    alpha_draw <- alpha_samples[draw, ]       # K-dimensional
    B_draw <- B_samples[draw, , ]             # K x ncov
    
    for (i in 1:n) {
      eta <- alpha_draw + x_temp[i, ] %*% t(B_draw)   # 1 x K
      theta_all_draws[i, , draw] <- exp(eta) / sum(exp(eta))
    }
  }
  
  # Compute mean and SE
  theta_mean <- apply(theta_all_draws, c(1,2), mean)
  theta_se   <- apply(theta_all_draws, c(1,2), sd) / sqrt(n_draws)
  
  df <- as.data.frame(theta_mean)
  df$X <- x_temp[, covar]
  df$Covariate <- paste0("Cov", covar)
  
  # Add SE for ribbon
  df_se <- as.data.frame(theta_se)
  names(df_se) <- paste0("SE_", names(df_se))
  df <- cbind(df, df_se)
  
  curve_list[[length(curve_list) + 1]] <- df
}

# Combine all covariates
data_long <- bind_rows(curve_list) %>%
  pivot_longer(cols = starts_with("V"), names_to = "Cluster", values_to = "Response") %>%
  pivot_longer(cols = starts_with("SE_V"), names_to = "SE_Cluster", values_to = "SE") %>%
  mutate(Cluster = factor(Cluster))

data_long
ggplot(data_long, aes(x = X, y = Response, group = Cluster)) +
  geom_line(aes(group = interaction(Covariate, Cluster)), alpha = 0.7) +
  geom_ribbon(aes(ymin = Response - SE, ymax = Response + SE, fill = Cluster), alpha = 0.2) +
  facet_grid(Cluster ~ Covariate, scales = "free_x") +
  theme_classic()+
  labs(y = "Δ Response", x = "Covariate Value") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  scale_fill_grey(start = 0.8, end = 0.5)

data_long$X<-as.numeric(data_long$X)
ggplot(data_long, aes(x = X, y = Response, group = Cluster)) +
  geom_line(aes(group = interaction(Covariate, Cluster)), alpha = 0.7) +
  geom_ribbon(aes(ymin = Response - SE, ymax = Response + SE, fill = Cluster), alpha = 0.2) +
  facet_grid(Cluster ~ Covariate, scales = "free_x") +
  theme_minimal(base_size = 13) +
  labs(y = "Δ Response", x = NULL) +  # remove x-axis title
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank()     # remove x-axis tick labels
  ) +
  scale_fill_grey(start = 0.8, end = 0.5)

# Combine all covariates
data_long <- bind_rows(curve_list)
data_long
# Pivot longer for mean response
mean_long <- data_long %>%
  pivot_longer(cols = starts_with("V"), names_to = "Cluster", values_to = "Response")

# Pivot longer for SE (make sure names match)
se_long <- data_long %>%
  pivot_longer(cols = starts_with("SE_V"), names_to = "Cluster", values_to = "SE") %>%
  mutate(Cluster = sub("SE_", "", Cluster))  # remove SE_ prefix to match Cluster names

# Join mean and SE
plot_data <- left_join(mean_long, se_long, by = c("X", "Covariate", "Cluster"))
plot_data
# Ensure X is numeric
plot_data$X <- as.numeric(plot_data$X)
plot_data$Cluster <- factor(plot_data$Cluster)
plot_data$Cluster
# Suppose your plot_data has a Cluster column like "V1", "V2", etc.
plot_data$Cluster <- factor(plot_data$Cluster, 
                            levels = paste0("V", 1:K),  # current levels
                            labels = paste("G", 1:K))  # new labels

cov_names <- c(
  Cov1= "STI-range",
  Cov2= "SPI-mean",
  Cov3 = "SPI-range",
  Cov4 = "Female size",
  Cov5 = "Adult",      
  Cov6= "Caterpillar",
  Cov7= "Egg",
  Cov8= "Pupa",
  Cov9="Monophagous",
  Cov10= "Mushroom",
  Cov11= "Oligophagous",          
  Cov12= "Polyphagous"
)           

plot_data$Covariate <- cov_names[plot_data$Covariate]

plot_data
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[1:ncov/3], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[5:8], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[9:ncov], ]
plot_data_subset




p3<-ggplot(plot_data_subset, aes(x = X, y = Response, group = Cluster)) +
  geom_ribbon(
    aes(ymin = Response - 2*SE, ymax = Response + 2*SE, fill = Cluster),
    alpha = 0.2
  ) +
  geom_line(
    aes(color = Cluster),
    alpha = 0.7
  ) +
  facet_grid(
    ~ Covariate,
    scales = "free_x",
    labeller = as_labeller(function(x) gsub("_prop$", "", x))
  ) +
  theme_classic() +
  labs(
    y = expression(Delta~Response),
    x = NULL,
    color = "Cluster",
    fill  = "Cluster"
    
  ) +
  theme(
    strip.text  = element_text(size = 8, face="bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 13),
    legend.position = "none"
  ) 


p2
p3
p4
#cov_names <-colnames(X_traits)[2:14]

curve_list <- list()
theta_samples_long <- NULL

for (covar in covariates_to_plot) {
  cat("Processing covariate", covar, "\n")
  
  covar_label <- cov_names[paste0("Cov", covar)]
  
  # ✅ Rebuild x_temp for this covariate
  x_temp <- matrix(0, n, ncov)
  for (j in 1:ncov) {
    x_temp[, j] <- mean(x_cov[, j])
  }
  x_temp[, covar] <- seq(min(x_cov[,covar]), max(x_cov[,covar]), length.out = n)
  
  # ✅ Recompute theta MEAN for this covariate
  theta_temp <- matrix(NA, n, K)
  for (i in 1:n) {
    eta <- as.numeric(alpha_mean + x_temp[i,] %*% t(B_matrix_mean))
    theta_temp[i, ] <- exp(eta) / sum(exp(eta))
  }
  
  df_mean <- as.data.frame(theta_temp)
  df_mean$X <- x_temp[, covar]
  df_mean$Covariate <- rep(covar_label, n)
  
  curve_list[[length(curve_list) + 1]] <- df_mean
  
  # ✅ Recompute FULL DRAW SAMPLES for this covariate
  sample_indices <- sample(1:n_draws, 20)  # reduce noise
  
  for (draw in sample_indices) {
    theta_draw <- matrix(NA, n, K)
    
    for (i in 1:n) {
      eta <- alpha_samples[draw,] + x_temp[i,] %*% t(B_samples[draw,,])
      theta_draw[i, ] <- exp(eta) / sum(exp(eta))
    }
    
    temp <- as.data.frame(theta_draw)
    temp$X <- x_temp[, covar]
    temp$Covariate <- rep(covar_label, n)
    temp$Draw <- draw
    
    theta_samples_long <- rbind(theta_samples_long, temp)
  }
}


plot_data <- bind_rows(curve_list) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "Cluster",
               values_to = "Response")

plot_data$Cluster <- factor(plot_data$Cluster,
                            levels = paste0("V", 1:K),
                            labels = paste("G", 1:K))


theta_samples_long <- theta_samples_long %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "Cluster",
               values_to = "Response")

theta_samples_long$Cluster <- factor(theta_samples_long$Cluster,
                                     levels = paste0("V", 1:K),
                                     labels = paste("G", 1:K))

theta_samples_long_subset <- theta_samples_long[theta_samples_long$Covariate %in% cov_names[1:4], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[1:4], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[7:12], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[13:ncov], ]


ggplot() +
  geom_point(data = theta_samples_long_subset,
             aes(x = X, y = Response, color = Cluster),
             size = 0.2, alpha = 0.4) +
  geom_line(data = plot_data_subset,
            aes(x = X, y = Response, group = Cluster),
            linewidth = 0.7, color = "black") +   # line on top
  facet_grid(Cluster ~ Covariate, scales = "free_x") +
  theme_classic() +
  labs(y = "Δ Response", x = NULL) +
  theme(
    strip.text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_color_grey(start = 0.6, end = 0.7)

###################################################################################
#########################For site case ############################################
#####site_abundance_covariates######################################################

df1<- df_covariates[,c(2:7)] %>%
  mutate(across(where(is.numeric), scales::rescale))

df2<-covariate_df[,c(2,3,4,5,6,7,8:13)]



x_cov<- bind_cols(df1, df2)

print(dim(x_cov))
head(x_cov)


##############################################################



cov_names2<-c( "BLforest_mean_prop",   "Cforest_mean_prop" ,   "Mforest_mean_prop"  ,  "SNHherb_mean_prop"  ,"clumpiness_mean_prop" , "entropy_mean_prop",    "BLforest_sd_prop" ,    "Cforest_sd_prop"  ,   "Mforest_sd_prop"  ,    "SNHherb_sd_prop" ,     "clumpiness_sd_prop",   "entropy_sd_prop")  
cov_names3<-colnames(df_covariates[2:7])
cov_names<- c( Cov1= cov_names3[1],
               Cov2= cov_names3[2],
               Cov3 = cov_names3[3],
               Cov4 = cov_names3[4],
               Cov5= cov_names3[5],
               Cov6 =cov_names3[6],Cov7=cov_names2[1],Cov8=cov_names2[2],Cov9=cov_names2[3],Cov10=cov_names2[4],Cov11=cov_names2[5],Cov12=cov_names2[6],Cov13=cov_names2[7],Cov14=cov_names2[8],Cov15=cov_names2[9],Cov16=cov_names2[10],Cov17=cov_names2[11],Cov18=cov_names2[12])
cov_names
cov_names <- gsub("_prop$", "", cov_names)


n<-d
K <- K# number of clusters
ncov<-dim(fit_ss$B)[3]
covariates_to_plot <-1:ncov 

plot_list <- list()
curve_list <- list()

x_cov<- as.matrix(x_cov)
ncov<-length(covariates_to_plot)
for (covar in covariates_to_plot) {
  
  x_temp <- matrix(0, n, ncov)
  
  # Set all covariates at mean
  for (j in 1:ncov) {
    x_temp[, j] <- mean(x_cov[, j])
  }
  
  # Vary only this covariate
  x_temp[, covar] <- seq(min(x_cov[,covar]), max(x_cov[,covar]), length.out = n)
  
  # Compute theta
  for (i in 1:n) {
    theta_temp <- matrix(NA, n, K)
    
    for (i in 1:n) {
      eta <- as.numeric(alpha_mean + x_temp[i,] %*% t(B_matrix_mean))
      theta_temp[i, ] <- exp(eta) / sum(exp(eta))
    }
  }
  
  df <- as.data.frame(theta_temp)
  df$X <- x_temp[, covar]
  df$Covariate <- paste0("Cov", covar)
  
  curve_list[[length(curve_list) + 1]] <- df
}

data_wide <- bind_rows(curve_list)
data_wide[, c("V1", "V2", "V3")]

data_long <- bind_rows(curve_list) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Cluster",
    values_to = "Response"
  ) %>%
  mutate(Cluster = factor(Cluster))

# Combine all covariates
data_long <- bind_rows(curve_list) %>%
  pivot_longer(cols = starts_with("V"), names_to = "Cluster", values_to = "Response") %>%
  mutate(Cluster = factor(Cluster))

# Plot
ggplot(data_long, aes(x = X, y = Response, group = Cluster)) +
  geom_line(aes(group = interaction(Covariate, Cluster)), alpha = 0.7) +
  facet_grid(Cluster ~ Covariate, scales = "free_x") +
  theme_minimal(base_size = 13) +
  labs(y = "Δ Response", x = "Covariate Value") +
  theme(strip.text = element_text(size = 12, face = "bold"))


covariates_to_plot <-1:ncov 
plot_list <- list()
curve_list <- list()

ncov<-length(covariates_to_plot)

curve_list <- list()
n_draws<-dim(fit_ss$alpha)[1]
alpha_samples<-fit_ss$alpha
B_samples <-result_b
r<-dim(result_b)[1]
ncov<-dim(x_cov)[2]
for (covar in covariates_to_plot) {
  # print(covar)
  
  x_temp <- matrix(0, n, ncov)
  
  # Set all covariates at mean
  for (j in 1:ncov) {
    x_temp[, j] <- mean(x_cov[, j])
  }
  
  # Vary only this covariate
  x_temp[, covar] <- seq(min(x_cov[,covar]), max(x_cov[,covar]), length.out = n)
  
  # Placeholder for storing theta from multiple draws
  theta_all_draws <- array(NA, dim = c(n, K, r))  # e.g., n_draws = 1000 posterior draws
  
  for (draw in 1:n_draws) {
    alpha_draw <- alpha_samples[draw, ]       # K-dimensional
    B_draw <- B_samples[draw, , ]             # K x ncov
    
    for (i in 1:n) {
      eta <- alpha_draw + x_temp[i, ] %*% t(B_draw)   # 1 x K
      theta_all_draws[i, , draw] <- exp(eta) / sum(exp(eta))
    }
  }
  
  # Compute mean and SE
  theta_mean <- apply(theta_all_draws, c(1,2), mean)
  theta_se   <- apply(theta_all_draws, c(1,2), sd) / sqrt(n_draws)
  
  df <- as.data.frame(theta_mean)
  df$X <- x_temp[, covar]
  df$Covariate <- paste0("Cov", covar)
  
  # Add SE for ribbon
  df_se <- as.data.frame(theta_se)
  names(df_se) <- paste0("SE_", names(df_se))
  df <- cbind(df, df_se)
  
  curve_list[[length(curve_list) + 1]] <- df
}

# Combine all covariates
data_long <- bind_rows(curve_list) %>%
  pivot_longer(cols = starts_with("V"), names_to = "Cluster", values_to = "Response") %>%
  pivot_longer(cols = starts_with("SE_V"), names_to = "SE_Cluster", values_to = "SE") %>%
  mutate(Cluster = factor(Cluster))

data_long
ggplot(data_long, aes(x = X, y = Response, group = Cluster)) +
  geom_line(aes(group = interaction(Covariate, Cluster)), alpha = 0.7) +
  geom_ribbon(aes(ymin = Response - SE, ymax = Response + SE, fill = Cluster), alpha = 0.2) +
  facet_grid(Cluster ~ Covariate, scales = "free_x") +
  theme_classic()+
  labs(y = "Δ Response", x = "Covariate Value") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  scale_fill_grey(start = 0.8, end = 0.5)

data_long$X<-as.numeric(data_long$X)
ggplot(data_long, aes(x = X, y = Response, group = Cluster)) +
  geom_line(aes(group = interaction(Covariate, Cluster)), alpha = 0.7) +
  geom_ribbon(aes(ymin = Response - SE, ymax = Response + SE, fill = Cluster), alpha = 0.2) +
  facet_grid(Cluster ~ Covariate, scales = "free_x") +
  theme_minimal(base_size = 13) +
  labs(y = "Δ Response", x = NULL) +  # remove x-axis title
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank()     # remove x-axis tick labels
  ) +
  scale_fill_grey(start = 0.8, end = 0.5)

# Combine all covariates
data_long <- bind_rows(curve_list)
data_long
# Pivot longer for mean response
mean_long <- data_long %>%
  pivot_longer(cols = starts_with("V"), names_to = "Cluster", values_to = "Response")

# Pivot longer for SE (make sure names match)
se_long <- data_long %>%
  pivot_longer(cols = starts_with("SE_V"), names_to = "Cluster", values_to = "SE") %>%
  mutate(Cluster = sub("SE_", "", Cluster))  # remove SE_ prefix to match Cluster names

# Join mean and SE
plot_data <- left_join(mean_long, se_long, by = c("X", "Covariate", "Cluster"))
plot_data
# Ensure X is numeric
plot_data$X <- as.numeric(plot_data$X)
plot_data$Cluster <- factor(plot_data$Cluster)
plot_data$Cluster
# Suppose your plot_data has a Cluster column like "V1", "V2", etc.
plot_data$Cluster <- factor(plot_data$Cluster, 
                            levels = paste0("V", 1:K),  # current levels
                            labels = paste("G", 1:K))  # new labels

plot_data$Covariate <- cov_names[plot_data$Covariate]

plot_data
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[1:ncov/3], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[5:8], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[9:ncov], ]
plot_data_subset




p3<-ggplot(plot_data_subset, aes(x = X, y = Response, group = Cluster)) +
  geom_ribbon(
    aes(ymin = Response - 2*SE, ymax = Response + 2*SE, fill = Cluster),
    alpha = 0.2
  ) +
  geom_line(
    aes(color = Cluster),
    alpha = 0.7
  ) +
  facet_grid(
    ~ Covariate,
    scales = "free_x",
    labeller = as_labeller(function(x) gsub("_prop$", "", x))
  ) +
  theme_classic() +
  labs(
    y = expression(Delta~Response),
    x = NULL,
    color = "Cluster",
    fill  = "Cluster"
    
  ) +
  theme(
    strip.text  = element_text(size = 8, face="bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 13),
    legend.position = "none"
  ) 


curve_list <- list()
theta_samples_long <- NULL

for (covar in covariates_to_plot) {
  cat("Processing covariate", covar, "\n")
  
  covar_label <- cov_names[paste0("Cov", covar)]
  
  # ✅ Rebuild x_temp for this covariate
  x_temp <- matrix(0, n, ncov)
  for (j in 1:ncov) {
    x_temp[, j] <- mean(x_cov[, j])
  }
  x_temp[, covar] <- seq(min(x_cov[,covar]), max(x_cov[,covar]), length.out = n)
  
  # ✅ Recompute theta MEAN for this covariate
  theta_temp <- matrix(NA, n, K)
  for (i in 1:n) {
    eta <- as.numeric(alpha_mean + x_temp[i,] %*% t(B_matrix_mean))
    theta_temp[i, ] <- exp(eta) / sum(exp(eta))
  }
  
  df_mean <- as.data.frame(theta_temp)
  df_mean$X <- x_temp[, covar]
  df_mean$Covariate <- rep(covar_label, n)
  
  curve_list[[length(curve_list) + 1]] <- df_mean
  
  # ✅ Recompute FULL DRAW SAMPLES for this covariate
  sample_indices <- sample(1:n_draws, 20)  # reduce noise
  
  for (draw in sample_indices) {
    theta_draw <- matrix(NA, n, K)
    
    for (i in 1:n) {
      eta <- alpha_samples[draw,] + x_temp[i,] %*% t(B_samples[draw,,])
      theta_draw[i, ] <- exp(eta) / sum(exp(eta))
    }
    
    temp <- as.data.frame(theta_draw)
    temp$X <- x_temp[, covar]
    temp$Covariate <- rep(covar_label, n)
    temp$Draw <- draw
    
    theta_samples_long <- rbind(theta_samples_long, temp)
  }
}


plot_data <- bind_rows(curve_list) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "Cluster",
               values_to = "Response")

plot_data$Cluster <- factor(plot_data$Cluster,
                            levels = paste0("V", 1:K),
                            labels = paste("G", 1:K))


theta_samples_long <- theta_samples_long %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "Cluster",
               values_to = "Response")

theta_samples_long$Cluster <- factor(theta_samples_long$Cluster,
                                     levels = paste0("V", 1:K),
                                     labels = paste("G", 1:K))

theta_samples_long_subset <- theta_samples_long[theta_samples_long$Covariate %in% cov_names[1:4], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[1:4], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[7:12], ]
plot_data_subset <- plot_data[plot_data$Covariate %in% cov_names[13:ncov], ]


ggplot() +
  geom_point(data = theta_samples_long_subset,
             aes(x = X, y = Response, color = Cluster),
             size = 0.2, alpha = 0.4) +
  geom_line(data = plot_data_subset,
            aes(x = X, y = Response, group = Cluster),
            linewidth = 0.7, color = "black") +   # line on top
  facet_grid(Cluster ~ Covariate, scales = "free_x") +
  theme_classic() +
  labs(y = "Δ Response", x = NULL) +
  theme(
    strip.text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_color_grey(start = 0.6, end = 0.7)


