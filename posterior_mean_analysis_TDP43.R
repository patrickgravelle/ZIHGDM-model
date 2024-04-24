library(ggplot2)
library(gridExtra)
library(dplyr)
library(mgcv)
library(bayesplot)
library(DirichletReg)


load("ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus_sml.Rdata")

# load("ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p96_taus_sml.Rdata")


##############################################################################################################

# Drawing from the posterior distributions

##############################################################################################################





## Steps to sample new outcomes from the posterior distribution

nmouse <- 1
n_day <- 1
nint <- 1000
nobs <- nmouse * 2 * n_day * 24
# Step 0 
# Create a covariate matrix to make this all make sense
intercept <- c(rep(1, nobs*nint))
genotype_tdp <- c(rep(0, nobs/2*nint), rep(1, nobs/2*nint))
gender <- c(rep(0, nobs/4*nint), rep(1, nobs/4*nint), rep(0, nobs/4*nint), rep(1, nobs/4*nint))
hour <- rep(c(0:23), nmouse*2*nint)
# batch <- rep(c(rep(1, nobs/6),rep(2, nobs/6)), 3)

covariate_mat <- data.frame(intercept, genotype_tdp, gender, hour)

# build the cyclic spline basis matrix
sspec <- s(hour, bs = "cp")
sm <- smoothCon(sspec, data = covariate_mat)[[1]]
bspline <- sm$X



# Use the correct posterior matrix


posterior_data <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus_sml



nsims <- 1000




p180_TDP_mean_HOUR_diffs <- array(NA, c(24, nsims, 9))
p180_TDP_mean_HOUR_TDP <- array(NA, c(24, nsims, 9))
p180_TDP_mean_HOUR_WT <- array(NA, c(24, nsims, 9))
p180_TDP_mean_FOURHOUR_diffs <- array(NA, c(6, nsims, 9))






for (m in 1:nsims){
  
  # Step 1
  # Sample one tau for each of the taus we have
  tau_id_mu <- posterior_data$tau_id_mu[m + 1500]
  tau_days_mu <- posterior_data$tau_days_mu[m + 1500]
  # tau_batch_mu <- posterior_data$tau_batch_mu[m + 1500]
  tau_id_pi <- posterior_data$tau_id_pi[m + 1500]
  tau_days_pi <- posterior_data$tau_days_pi[m + 1500]
  # tau_batch_pi <- posterior_data$tau_batch_pi[m + 1500]
  
  
  
  
  # Step 2a
  # Draw 30 mice effects
  # Draw 5 day effect
  # Draw 2 batch effects 
  # each for pi and mu and all using the taus from above
  id_mu <- rnorm(nint, 0, tau_id_mu)
  days_mu <- rnorm(nint, 0, tau_days_mu)
  # batch_mu <- rnorm(nint, 0, tau_batch_mu)
  
  id_pi <- rnorm(nint, 0, tau_id_pi)
  days_pi <- rnorm(nint, 0, tau_days_pi)
  # batch_pi <- rnorm(nint, 0, tau_batch_pi)
  
  
  # Step 2b
  # Create vectors of these effects to add to the mean models
  id_mu_vec <- rep(rep(id_mu, each = nobs/2), 2)
  days_mu_vec <- rep(rep(days_mu, each = nobs/2), 2)
  # batch_mu_vec <- rep(rep(batch_mu, each = nobs/2), 2)
  
  id_pi_vec <- rep(rep(id_pi, each = nobs/2), 2)
  days_pi_vec <- rep(rep(days_pi, each = nobs/2), 2)
  # batch_pi_vec <- rep(rep(batch_pi, each = nobs/2), 2)
  
  
  # Step 3
  # Sample a beta, a vector of "a"'s, and an A for each of these different variables
  # Need to do this for each behaviour
  
  # Beta 0
  beta0_mu <- matrix(NA, nrow = 8, ncol = 1)
  beta0_sig <- matrix(NA, nrow = 8, ncol = 1)
  beta0_pi <- matrix(NA, nrow = 8, ncol = 1)
  
  beta0_mu[1,1] <- posterior_data$params_mu.1.1.[m + 1500]
  beta0_mu[2,1] <- posterior_data$params_mu.2.1.[m + 1500]
  beta0_mu[3,1] <- posterior_data$params_mu.3.1.[m + 1500]
  beta0_mu[4,1] <- posterior_data$params_mu.4.1.[m + 1500]
  beta0_mu[5,1] <- posterior_data$params_mu.5.1.[m + 1500]
  beta0_mu[6,1] <- posterior_data$params_mu.6.1.[m + 1500]
  beta0_mu[7,1] <- posterior_data$params_mu.7.1.[m + 1500]
  beta0_mu[8,1] <- posterior_data$params_mu.8.1.[m + 1500]
  
  
  beta0_sig[1,1] <- posterior_data$params_sig.1.1.[m + 1500]
  beta0_sig[2,1] <- posterior_data$params_sig.2.1.[m + 1500]
  beta0_sig[3,1] <- posterior_data$params_sig.3.1.[m + 1500]
  beta0_sig[4,1] <- posterior_data$params_sig.4.1.[m + 1500]
  beta0_sig[5,1] <- posterior_data$params_sig.5.1.[m + 1500]
  beta0_sig[6,1] <- posterior_data$params_sig.6.1.[m + 1500]
  beta0_sig[7,1] <- posterior_data$params_sig.7.1.[m + 1500]
  beta0_sig[8,1] <- posterior_data$params_sig.8.1.[m + 1500]
  
  
  beta0_pi[1,1] <- posterior_data$params_pi.1.1.[m + 1500]
  beta0_pi[2,1] <- posterior_data$params_pi.2.1.[m + 1500]
  beta0_pi[3,1] <- posterior_data$params_pi.3.1.[m + 1500]
  beta0_pi[4,1] <- posterior_data$params_pi.4.1.[m + 1500]
  beta0_pi[5,1] <- posterior_data$params_pi.5.1.[m + 1500]
  beta0_pi[6,1] <- posterior_data$params_pi.6.1.[m + 1500]
  beta0_pi[7,1] <- posterior_data$params_pi.7.1.[m + 1500]
  beta0_pi[8,1] <- posterior_data$params_pi.8.1.[m + 1500]
  
  
  
  # Beta 1
  beta1_mu <- matrix(NA, nrow = 8, ncol = 1)
  beta1_sig <- matrix(NA, nrow = 8, ncol = 1)
  beta1_pi <- matrix(NA, nrow = 8, ncol = 1)
  
  beta1_mu[1,1] <- posterior_data$params_mu.1.2.[m + 1500]
  beta1_mu[2,1] <- posterior_data$params_mu.2.2.[m + 1500]
  beta1_mu[3,1] <- posterior_data$params_mu.3.2.[m + 1500]
  beta1_mu[4,1] <- posterior_data$params_mu.4.2.[m + 1500]
  beta1_mu[5,1] <- posterior_data$params_mu.5.2.[m + 1500]
  beta1_mu[6,1] <- posterior_data$params_mu.6.2.[m + 1500]
  beta1_mu[7,1] <- posterior_data$params_mu.7.2.[m + 1500]
  beta1_mu[8,1] <- posterior_data$params_mu.8.2.[m + 1500]
  
  
  beta1_sig[1,1] <- posterior_data$params_sig.1.2.[m + 1500]
  beta1_sig[2,1] <- posterior_data$params_sig.2.2.[m + 1500]
  beta1_sig[3,1] <- posterior_data$params_sig.3.2.[m + 1500]
  beta1_sig[4,1] <- posterior_data$params_sig.4.2.[m + 1500]
  beta1_sig[5,1] <- posterior_data$params_sig.5.2.[m + 1500]
  beta1_sig[6,1] <- posterior_data$params_sig.6.2.[m + 1500]
  beta1_sig[7,1] <- posterior_data$params_sig.7.2.[m + 1500]
  beta1_sig[8,1] <- posterior_data$params_sig.8.2.[m + 1500]
  
  
  beta1_pi[1,1] <- posterior_data$params_pi.1.2.[m + 1500]
  beta1_pi[2,1] <- posterior_data$params_pi.2.2.[m + 1500]
  beta1_pi[3,1] <- posterior_data$params_pi.3.2.[m + 1500]
  beta1_pi[4,1] <- posterior_data$params_pi.4.2.[m + 1500]
  beta1_pi[5,1] <- posterior_data$params_pi.5.2.[m + 1500]
  beta1_pi[6,1] <- posterior_data$params_pi.6.2.[m + 1500]
  beta1_pi[7,1] <- posterior_data$params_pi.7.2.[m + 1500]
  beta1_pi[8,1] <- posterior_data$params_pi.8.2.[m + 1500]
  
  
  
  # Beta 2
  beta2_mu <- matrix(NA, nrow = 8, ncol = 1)
  beta2_sig <- matrix(NA, nrow = 8, ncol = 1)
  beta2_pi <- matrix(NA, nrow = 8, ncol = 1)
  
  beta2_mu[1,1] <- posterior_data$params_mu.1.3.[m + 1500]
  beta2_mu[2,1] <- posterior_data$params_mu.2.3.[m + 1500]
  beta2_mu[3,1] <- posterior_data$params_mu.3.3.[m + 1500]
  beta2_mu[4,1] <- posterior_data$params_mu.4.3.[m + 1500]
  beta2_mu[5,1] <- posterior_data$params_mu.5.3.[m + 1500]
  beta2_mu[6,1] <- posterior_data$params_mu.6.3.[m + 1500]
  beta2_mu[7,1] <- posterior_data$params_mu.7.3.[m + 1500]
  beta2_mu[8,1] <- posterior_data$params_mu.8.3.[m + 1500]
  
  
  beta2_sig[1,1] <- posterior_data$params_sig.1.3.[m + 1500]
  beta2_sig[2,1] <- posterior_data$params_sig.2.3.[m + 1500]
  beta2_sig[3,1] <- posterior_data$params_sig.3.3.[m + 1500]
  beta2_sig[4,1] <- posterior_data$params_sig.4.3.[m + 1500]
  beta2_sig[5,1] <- posterior_data$params_sig.5.3.[m + 1500]
  beta2_sig[6,1] <- posterior_data$params_sig.6.3.[m + 1500]
  beta2_sig[7,1] <- posterior_data$params_sig.7.3.[m + 1500]
  beta2_sig[8,1] <- posterior_data$params_sig.8.3.[m + 1500]
  
  
  beta2_pi[1,1] <- posterior_data$params_pi.1.3.[m + 1500]
  beta2_pi[2,1] <- posterior_data$params_pi.2.3.[m + 1500]
  beta2_pi[3,1] <- posterior_data$params_pi.3.3.[m + 1500]
  beta2_pi[4,1] <- posterior_data$params_pi.4.3.[m + 1500]
  beta2_pi[5,1] <- posterior_data$params_pi.5.3.[m + 1500]
  beta2_pi[6,1] <- posterior_data$params_pi.6.3.[m + 1500]
  beta2_pi[7,1] <- posterior_data$params_pi.7.3.[m + 1500]
  beta2_pi[8,1] <- posterior_data$params_pi.8.3.[m + 1500]
  
  
  
  # Spline Covariates Mu & Pi
  
  spline_a11_mu <- matrix(NA, nrow = 10, ncol = 8)
  spline_a13_pi <- matrix(NA, nrow = 10, ncol = 8)
  
  spline_a11_mu[1,1] <- posterior_data$a1_1.1.[m + 1500]
  spline_a11_mu[2,1] <- posterior_data$a1_1.2.[m + 1500]
  spline_a11_mu[3,1] <- posterior_data$a1_1.3.[m + 1500]
  spline_a11_mu[4,1] <- posterior_data$a1_1.4.[m + 1500]
  spline_a11_mu[5,1] <- posterior_data$a1_1.5.[m + 1500]
  spline_a11_mu[6,1] <- posterior_data$a1_1.6.[m + 1500]
  spline_a11_mu[7,1] <- posterior_data$a1_1.7.[m + 1500]
  spline_a11_mu[8,1] <- posterior_data$a1_1.8.[m + 1500]
  spline_a11_mu[9,1] <- posterior_data$a1_1.9.[m + 1500]
  spline_a11_mu[10,1] <- posterior_data$a1_1.10.[m + 1500]
  
  # eat 2
  spline_a11_mu[1,2] <- posterior_data$a2_1.1.[m + 1500]
  spline_a11_mu[2,2] <- posterior_data$a2_1.2.[m + 1500]
  spline_a11_mu[3,2] <- posterior_data$a2_1.3.[m + 1500]
  spline_a11_mu[4,2] <- posterior_data$a2_1.4.[m + 1500]
  spline_a11_mu[5,2] <- posterior_data$a2_1.5.[m + 1500]
  spline_a11_mu[6,2] <- posterior_data$a2_1.6.[m + 1500]
  spline_a11_mu[7,2] <- posterior_data$a2_1.7.[m + 1500]
  spline_a11_mu[8,2] <- posterior_data$a2_1.8.[m + 1500]
  spline_a11_mu[9,2] <- posterior_data$a2_1.9.[m + 1500]
  spline_a11_mu[10,2] <- posterior_data$a2_1.10.[m + 1500]
  
  # ebh 3
  spline_a11_mu[1,3] <- posterior_data$a3_1.1.[m + 1500]
  spline_a11_mu[2,3] <- posterior_data$a3_1.2.[m + 1500]
  spline_a11_mu[3,3] <- posterior_data$a3_1.3.[m + 1500]
  spline_a11_mu[4,3] <- posterior_data$a3_1.4.[m + 1500]
  spline_a11_mu[5,3] <- posterior_data$a3_1.5.[m + 1500]
  spline_a11_mu[6,3] <- posterior_data$a3_1.6.[m + 1500]
  spline_a11_mu[7,3] <- posterior_data$a3_1.7.[m + 1500]
  spline_a11_mu[8,3] <- posterior_data$a3_1.8.[m + 1500]
  spline_a11_mu[9,3] <- posterior_data$a3_1.9.[m + 1500]
  spline_a11_mu[10,3] <- posterior_data$a3_1.10.[m + 1500]
  
  # groom 4
  spline_a11_mu[1,4] <- posterior_data$a4_1.1.[m + 1500]
  spline_a11_mu[2,4] <- posterior_data$a4_1.2.[m + 1500]
  spline_a11_mu[3,4] <- posterior_data$a4_1.3.[m + 1500]
  spline_a11_mu[4,4] <- posterior_data$a4_1.4.[m + 1500]
  spline_a11_mu[5,4] <- posterior_data$a4_1.5.[m + 1500]
  spline_a11_mu[6,4] <- posterior_data$a4_1.6.[m + 1500]
  spline_a11_mu[7,4] <- posterior_data$a4_1.7.[m + 1500]
  spline_a11_mu[8,4] <- posterior_data$a4_1.8.[m + 1500]
  spline_a11_mu[9,4] <- posterior_data$a4_1.9.[m + 1500]
  spline_a11_mu[10,4] <- posterior_data$a4_1.10.[m + 1500]
  
  # hang 5
  spline_a11_mu[1,5] <- posterior_data$a5_1.1.[m + 1500]
  spline_a11_mu[2,5] <- posterior_data$a5_1.2.[m + 1500]
  spline_a11_mu[3,5] <- posterior_data$a5_1.3.[m + 1500]
  spline_a11_mu[4,5] <- posterior_data$a5_1.4.[m + 1500]
  spline_a11_mu[5,5] <- posterior_data$a5_1.5.[m + 1500]
  spline_a11_mu[6,5] <- posterior_data$a5_1.6.[m + 1500]
  spline_a11_mu[7,5] <- posterior_data$a5_1.7.[m + 1500]
  spline_a11_mu[8,5] <- posterior_data$a5_1.8.[m + 1500]
  spline_a11_mu[9,5] <- posterior_data$a5_1.9.[m + 1500]
  spline_a11_mu[10,5] <- posterior_data$a5_1.10.[m + 1500]
  
  # rear 6
  spline_a11_mu[1,6] <- posterior_data$a6_1.1.[m + 1500]
  spline_a11_mu[2,6] <- posterior_data$a6_1.2.[m + 1500]
  spline_a11_mu[3,6] <- posterior_data$a6_1.3.[m + 1500]
  spline_a11_mu[4,6] <- posterior_data$a6_1.4.[m + 1500]
  spline_a11_mu[5,6] <- posterior_data$a6_1.5.[m + 1500]
  spline_a11_mu[6,6] <- posterior_data$a6_1.6.[m + 1500]
  spline_a11_mu[7,6] <- posterior_data$a6_1.7.[m + 1500]
  spline_a11_mu[8,6] <- posterior_data$a6_1.8.[m + 1500]
  spline_a11_mu[9,6] <- posterior_data$a6_1.9.[m + 1500]
  spline_a11_mu[10,6] <- posterior_data$a6_1.10.[m + 1500]
  
  # rest 7
  spline_a11_mu[1,7] <- posterior_data$a7_1.1.[m + 1500]
  spline_a11_mu[2,7] <- posterior_data$a7_1.2.[m + 1500]
  spline_a11_mu[3,7] <- posterior_data$a7_1.3.[m + 1500]
  spline_a11_mu[4,7] <- posterior_data$a7_1.4.[m + 1500]
  spline_a11_mu[5,7] <- posterior_data$a7_1.5.[m + 1500]
  spline_a11_mu[6,7] <- posterior_data$a7_1.6.[m + 1500]
  spline_a11_mu[7,7] <- posterior_data$a7_1.7.[m + 1500]
  spline_a11_mu[8,7] <- posterior_data$a7_1.8.[m + 1500]
  spline_a11_mu[9,7] <- posterior_data$a7_1.9.[m + 1500]
  spline_a11_mu[10,7] <- posterior_data$a7_1.10.[m + 1500]
  
  # sniff 8
  spline_a11_mu[1,8] <- posterior_data$a8_1.1.[m + 1500]
  spline_a11_mu[2,8] <- posterior_data$a8_1.2.[m + 1500]
  spline_a11_mu[3,8] <- posterior_data$a8_1.3.[m + 1500]
  spline_a11_mu[4,8] <- posterior_data$a8_1.4.[m + 1500]
  spline_a11_mu[5,8] <- posterior_data$a8_1.5.[m + 1500]
  spline_a11_mu[6,8] <- posterior_data$a8_1.6.[m + 1500]
  spline_a11_mu[7,8] <- posterior_data$a8_1.7.[m + 1500]
  spline_a11_mu[8,8] <- posterior_data$a8_1.8.[m + 1500]
  spline_a11_mu[9,8] <- posterior_data$a8_1.9.[m + 1500]
  spline_a11_mu[10,8] <- posterior_data$a8_1.10.[m + 1500]
  
  
  # the pi splines
  # Drink 1 
  
  spline_a13_pi[1,1] <- posterior_data$a1_3.1.[m + 1500]
  spline_a13_pi[2,1] <- posterior_data$a1_3.2.[m + 1500]
  spline_a13_pi[3,1] <- posterior_data$a1_3.3.[m + 1500]
  spline_a13_pi[4,1] <- posterior_data$a1_3.4.[m + 1500]
  spline_a13_pi[5,1] <- posterior_data$a1_3.5.[m + 1500]
  spline_a13_pi[6,1] <- posterior_data$a1_3.6.[m + 1500]
  spline_a13_pi[7,1] <- posterior_data$a1_3.7.[m + 1500]
  spline_a13_pi[8,1] <- posterior_data$a1_3.8.[m + 1500]
  spline_a13_pi[9,1] <- posterior_data$a1_3.9.[m + 1500]
  spline_a13_pi[10,1] <- posterior_data$a1_3.10.[m + 1500]
  
  # eat 2
  spline_a13_pi[1,2] <- posterior_data$a2_3.1.[m + 1500]
  spline_a13_pi[2,2] <- posterior_data$a2_3.2.[m + 1500]
  spline_a13_pi[3,2] <- posterior_data$a2_3.3.[m + 1500]
  spline_a13_pi[4,2] <- posterior_data$a2_3.4.[m + 1500]
  spline_a13_pi[5,2] <- posterior_data$a2_3.5.[m + 1500]
  spline_a13_pi[6,2] <- posterior_data$a2_3.6.[m + 1500]
  spline_a13_pi[7,2] <- posterior_data$a2_3.7.[m + 1500]
  spline_a13_pi[8,2] <- posterior_data$a2_3.8.[m + 1500]
  spline_a13_pi[9,2] <- posterior_data$a2_3.9.[m + 1500]
  spline_a13_pi[10,2] <- posterior_data$a2_3.10.[m + 1500]
  
  # ebh 3
  spline_a13_pi[1,3] <- posterior_data$a3_3.1.[m + 1500]
  spline_a13_pi[2,3] <- posterior_data$a3_3.2.[m + 1500]
  spline_a13_pi[3,3] <- posterior_data$a3_3.3.[m + 1500]
  spline_a13_pi[4,3] <- posterior_data$a3_3.4.[m + 1500]
  spline_a13_pi[5,3] <- posterior_data$a3_3.5.[m + 1500]
  spline_a13_pi[6,3] <- posterior_data$a3_3.6.[m + 1500]
  spline_a13_pi[7,3] <- posterior_data$a3_3.7.[m + 1500]
  spline_a13_pi[8,3] <- posterior_data$a3_3.8.[m + 1500]
  spline_a13_pi[9,3] <- posterior_data$a3_3.9.[m + 1500]
  spline_a13_pi[10,3] <- posterior_data$a3_3.10.[m + 1500]
  
  # groom 4
  spline_a13_pi[1,4] <- posterior_data$a4_3.1.[m + 1500]
  spline_a13_pi[2,4] <- posterior_data$a4_3.2.[m + 1500]
  spline_a13_pi[3,4] <- posterior_data$a4_3.3.[m + 1500]
  spline_a13_pi[4,4] <- posterior_data$a4_3.4.[m + 1500]
  spline_a13_pi[5,4] <- posterior_data$a4_3.5.[m + 1500]
  spline_a13_pi[6,4] <- posterior_data$a4_3.6.[m + 1500]
  spline_a13_pi[7,4] <- posterior_data$a4_3.7.[m + 1500]
  spline_a13_pi[8,4] <- posterior_data$a4_3.8.[m + 1500]
  spline_a13_pi[9,4] <- posterior_data$a4_3.9.[m + 1500]
  spline_a13_pi[10,4] <- posterior_data$a4_3.10.[m + 1500]
  
  # hang 5
  spline_a13_pi[1,5] <- posterior_data$a5_3.1.[m + 1500]
  spline_a13_pi[2,5] <- posterior_data$a5_3.2.[m + 1500]
  spline_a13_pi[3,5] <- posterior_data$a5_3.3.[m + 1500]
  spline_a13_pi[4,5] <- posterior_data$a5_3.4.[m + 1500]
  spline_a13_pi[5,5] <- posterior_data$a5_3.5.[m + 1500]
  spline_a13_pi[6,5] <- posterior_data$a5_3.6.[m + 1500]
  spline_a13_pi[7,5] <- posterior_data$a5_3.7.[m + 1500]
  spline_a13_pi[8,5] <- posterior_data$a5_3.8.[m + 1500]
  spline_a13_pi[9,5] <- posterior_data$a5_3.9.[m + 1500]
  spline_a13_pi[10,5] <- posterior_data$a5_3.10.[m + 1500]
  
  # rear 6
  spline_a13_pi[1,6] <- posterior_data$a6_3.1.[m + 1500]
  spline_a13_pi[2,6] <- posterior_data$a6_3.2.[m + 1500]
  spline_a13_pi[3,6] <- posterior_data$a6_3.3.[m + 1500]
  spline_a13_pi[4,6] <- posterior_data$a6_3.4.[m + 1500]
  spline_a13_pi[5,6] <- posterior_data$a6_3.5.[m + 1500]
  spline_a13_pi[6,6] <- posterior_data$a6_3.6.[m + 1500]
  spline_a13_pi[7,6] <- posterior_data$a6_3.7.[m + 1500]
  spline_a13_pi[8,6] <- posterior_data$a6_3.8.[m + 1500]
  spline_a13_pi[9,6] <- posterior_data$a6_3.9.[m + 1500]
  spline_a13_pi[10,6] <- posterior_data$a6_3.10.[m + 1500]
  
  # rest 7
  spline_a13_pi[1,7] <- posterior_data$a7_3.1.[m + 1500]
  spline_a13_pi[2,7] <- posterior_data$a7_3.2.[m + 1500]
  spline_a13_pi[3,7] <- posterior_data$a7_3.3.[m + 1500]
  spline_a13_pi[4,7] <- posterior_data$a7_3.4.[m + 1500]
  spline_a13_pi[5,7] <- posterior_data$a7_3.5.[m + 1500]
  spline_a13_pi[6,7] <- posterior_data$a7_3.6.[m + 1500]
  spline_a13_pi[7,7] <- posterior_data$a7_3.7.[m + 1500]
  spline_a13_pi[8,7] <- posterior_data$a7_3.8.[m + 1500]
  spline_a13_pi[9,7] <- posterior_data$a7_3.9.[m + 1500]
  spline_a13_pi[10,7] <- posterior_data$a7_3.10.[m + 1500]
  
  # sniff 8
  spline_a13_pi[1,8] <- posterior_data$a8_3.1.[m + 1500]
  spline_a13_pi[2,8] <- posterior_data$a8_3.2.[m + 1500]
  spline_a13_pi[3,8] <- posterior_data$a8_3.3.[m + 1500]
  spline_a13_pi[4,8] <- posterior_data$a8_3.4.[m + 1500]
  spline_a13_pi[5,8] <- posterior_data$a8_3.5.[m + 1500]
  spline_a13_pi[6,8] <- posterior_data$a8_3.6.[m + 1500]
  spline_a13_pi[7,8] <- posterior_data$a8_3.7.[m + 1500]
  spline_a13_pi[8,8] <- posterior_data$a8_3.8.[m + 1500]
  spline_a13_pi[9,8] <- posterior_data$a8_3.9.[m + 1500]
  spline_a13_pi[10,8] <- posterior_data$a8_3.10.[m + 1500]
  
  
  
  # The interaction effect between genotype and spline
  
  # Interactions
  interaction1_mu <- matrix(NA, nrow = 8, ncol = 1)

  interaction3_pi <- matrix(NA, nrow = 8, ncol = 1)

  
  interaction1_mu[1,1] <- posterior_data$A1.1.[m + 1500]
  interaction1_mu[2,1] <- posterior_data$A1.2.[m + 1500]
  interaction1_mu[3,1] <- posterior_data$A1.3.[m + 1500]
  interaction1_mu[4,1] <- posterior_data$A1.4.[m + 1500]
  interaction1_mu[5,1] <- posterior_data$A1.5.[m + 1500]
  interaction1_mu[6,1] <- posterior_data$A1.6.[m + 1500]
  interaction1_mu[7,1] <- posterior_data$A1.7.[m + 1500]
  interaction1_mu[8,1] <- posterior_data$A1.8.[m + 1500]
  
  
  interaction3_pi[1,1] <- posterior_data$A3.1.[m + 1500]
  interaction3_pi[2,1] <- posterior_data$A3.2.[m + 1500]
  interaction3_pi[3,1] <- posterior_data$A3.3.[m + 1500]
  interaction3_pi[4,1] <- posterior_data$A3.4.[m + 1500]
  interaction3_pi[5,1] <- posterior_data$A3.5.[m + 1500]
  interaction3_pi[6,1] <- posterior_data$A3.6.[m + 1500]
  interaction3_pi[7,1] <- posterior_data$A3.7.[m + 1500]
  interaction3_pi[8,1] <- posterior_data$A3.8.[m + 1500]
  
  
  
  # Step 4a Compute the linear predictor for each mu, pi, and sigma
  
  logit_mu_matrix <- matrix(NA, nrow = nobs*nint, ncol = 8)
  logit_pi_matrix <- matrix(NA, nrow = nobs*nint, ncol = 8)
  logit_sig_matrix <- matrix(NA, nrow = nobs*nint, ncol =8)
  
  
  
  for (j in 1:8){
    logit_mu_matrix[,j] <- beta0_mu[j,]*covariate_mat$intercept + beta1_mu[j,]*covariate_mat$gender + 
      beta2_mu[j,]*covariate_mat$genotype_tdp + bspline %*% spline_a11_mu[,j] + 
      interaction1_mu[j,] * covariate_mat$genotype_tdp * bspline %*% spline_a11_mu[,j] + 
      id_mu_vec + days_mu_vec
    
    logit_pi_matrix[,j] <- beta0_pi[j,]*covariate_mat$intercept + beta1_pi[j,]*covariate_mat$gender + 
      beta2_pi[j,]*covariate_mat$genotype_tdp + bspline %*% spline_a13_pi[,j] + 
      interaction3_pi[j,] * covariate_mat$genotype_tdp * bspline %*% spline_a13_pi[,j] + 
      id_pi_vec + days_pi_vec
    
    logit_sig_matrix[,j] <- beta0_sig[j,]*covariate_mat$intercept + beta1_sig[j,]*covariate_mat$gender + 
      beta2_sig[j,]*covariate_mat$genotype_tdp
  }
  
  
  # need to take the means of each of these logits above
  # except the mean has to be separate for each genotype
  # so we should have two values afterwards for each genotype
  mu_matrix <- matrix(NA, nrow = nobs, ncol = 8)
  pi_matrix <- matrix(NA, nrow = nobs, ncol = 8)
  sig_matrix <- matrix(NA, nrow = nobs, ncol =8)
  
  
  logit_mu_matrix <- data.frame(logit_mu_matrix)
  logit_pi_matrix <- data.frame(logit_pi_matrix)
  logit_sig_matrix <- data.frame(logit_sig_matrix)
  
  logit_mu_matrix$genotype <- covariate_mat$genotype_tdp
  logit_mu_matrix$hour <- covariate_mat$hour
  
  logit_pi_matrix$genotype <- covariate_mat$genotype_tdp
  logit_pi_matrix$hour <- covariate_mat$hour
  
  logit_sig_matrix$genotype <- covariate_mat$genotype_tdp
  logit_sig_matrix$hour <- covariate_mat$hour
  
  
  for (g in 0:1){
    geno_mu <- logit_mu_matrix %>% filter(genotype == g)
    geno_pi <- logit_pi_matrix %>% filter(genotype == g)
    geno_sig <- logit_sig_matrix %>% filter(genotype == g)
    
    for (t in 0:23){
      hour_mu <- geno_mu %>% filter(hour == t)
      hour_pi <- geno_pi %>% filter(hour == t)
      hour_sig <- geno_sig %>% filter(hour == t)
      
      for(j in 1:8){
        mu_matrix[((t+1) + (g*24)),j] <- mean(hour_mu[,j])
        pi_matrix[((t+1) + (g*24)),j] <- mean(hour_pi[,j])
        sig_matrix[((t+1) + (g*24)),j] <- mean(hour_sig[,j])
      }
      
    }
    
  }
  
  mu_matrix <- exp(mu_matrix) / (1 + exp(mu_matrix))
  pi_matrix <- exp(pi_matrix) / (1 + exp(pi_matrix))
  sig_matrix <- exp(sig_matrix) / (1 + exp(sig_matrix))
  
  
  # compute the mean probability at each hour for each behavior and each genotype
  prob_matrix <- matrix(NA, nrow = nobs, ncol = 9)
  
  
  for (n in 1:nobs){
    prob_matrix[n,1] <- mu_matrix[n,1] * (1 - pi_matrix[n,1])
    prob_matrix[n,2] <- mu_matrix[n,2] * (1 - pi_matrix[n,2]) * (1 - mu_matrix[n,1])
    prob_matrix[n,3] <- mu_matrix[n,3] * (1 - pi_matrix[n,3]) * (1 - mu_matrix[n,1]) * (1 - mu_matrix[n,2])
    prob_matrix[n,4] <- mu_matrix[n,4] * (1 - pi_matrix[n,4]) * (1 - mu_matrix[n,1]) * (1 - mu_matrix[n,2]) * (1 - mu_matrix[n,3])
    prob_matrix[n,5] <- mu_matrix[n,5] * (1 - pi_matrix[n,5]) * (1 - mu_matrix[n,1]) * (1 - mu_matrix[n,2]) * (1 - mu_matrix[n,3]) * (1 - mu_matrix[n,4])
    prob_matrix[n,6] <- mu_matrix[n,6] * (1 - pi_matrix[n,6]) * (1 - mu_matrix[n,1]) * (1 - mu_matrix[n,2]) * (1 - mu_matrix[n,3]) * (1 - mu_matrix[n,4]) * (1 - mu_matrix[n,5])
    prob_matrix[n,7] <- mu_matrix[n,7] * (1 - pi_matrix[n,7]) * (1 - mu_matrix[n,1]) * (1 - mu_matrix[n,2]) * (1 - mu_matrix[n,3]) * (1 - mu_matrix[n,4]) * (1 - mu_matrix[n,5]) * (1 - mu_matrix[n,6])
    prob_matrix[n,8] <- mu_matrix[n,8] * (1 - pi_matrix[n,8]) * (1 - mu_matrix[n,1]) * (1 - mu_matrix[n,2]) * (1 - mu_matrix[n,3]) * (1 - mu_matrix[n,4]) * (1 - mu_matrix[n,5]) * (1 - mu_matrix[n,6]) * (1 - mu_matrix[n,7])
    prob_matrix[n,9] <- 1 - sum(prob_matrix[n,1:8])
  }
  
  
  # Step 10
  # compute the outcomes
  y_matrix <- matrix(NA, nrow = nobs, ncol = 9)
  
  y_matrix <- 3600 * prob_matrix
  
  
  
  # summary(y_matrix)
  
  
  
  
  
  # Step 11
  
  ### Make this an hourly difference matrix
  ### So first we just make this matrix the same dimensions as our y_matrix
  hourly_mat <- matrix(NA, nrow = (nobs), ncol = 9)
  
  hourly_mat <- y_matrix
  
  # Total Day Separated by Genotypes
  hourly_mat_WT <- hourly_mat[1:(nmouse*n_day*24),]
  hourly_mat_TDP <- hourly_mat[(1:(nmouse*n_day*24) + (nmouse*n_day*24)),]
  # Total Day Differences from WT by HOM and HET
  hourly_DIFF_WT_TDP <- hourly_mat_TDP - hourly_mat_WT
  
  
  hourVec <- 0:23
  
  hourly_DIFF_WT_TDP <- data.frame(hourly_DIFF_WT_TDP, hourVec)
  
  names(hourly_DIFF_WT_TDP) <- c("Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Hour")
  
  # TDP estimates
  hourly_EST_TDP <- data.frame(hourly_mat_TDP, hourVec)
  names(hourly_EST_TDP) <- c("Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Hour")
  
  # WT estimates
  hourly_EST_WT <- data.frame(hourly_mat_WT, hourVec)
  names(hourly_EST_WT) <- c("Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Hour")
  
  for (t in 0:23){
    hour_diff <- hourly_DIFF_WT_TDP %>% filter(Hour == t)
    hour_TDP <- hourly_EST_TDP %>% filter(Hour == t)
    hour_WT <- hourly_EST_WT %>% filter(Hour == t)
    
    for(j in 1:9){
      p180_TDP_mean_HOUR_diffs[(t+1),m,j] <- hour_diff[,j]
      p180_TDP_mean_HOUR_TDP[(t+1),m,j] <- hour_TDP[,j]
      p180_TDP_mean_HOUR_WT[(t+1),m,j] <- hour_WT[,j]
    }
    
  }
  
  
}




p180_TDP_mean_HOUR_diffs_all <- matrix(NA, nrow = (24*nsims), ncol = 9)
for (t in 1:24){
  p180_TDP_mean_HOUR_diffs_all[((1:nsims) + (t-1)*nsims), ] <- p180_TDP_mean_HOUR_diffs[t, , ]
}

# TDP Hourly
p180_TDP_mean_HOUR_TDP_all <- matrix(NA, nrow = (24*nsims), ncol = 9)
for (t in 1:24){
  p180_TDP_mean_HOUR_TDP_all[((1:nsims) + (t-1)*nsims), ] <- p180_TDP_mean_HOUR_TDP[t, , ]
}

# WT hourly
p180_TDP_mean_HOUR_WT_all <- matrix(NA, nrow = (24*nsims), ncol = 9)
for (t in 1:24){
  p180_TDP_mean_HOUR_WT_all[((1:nsims) + (t-1)*nsims), ] <- p180_TDP_mean_HOUR_WT[t, , ]
}





save(p180_TDP_mean_HOUR_diffs_all, file = "p180_TDP_mean_HOUR_diffs_all.Rdata")

# TDP save
save(p180_TDP_mean_HOUR_TDP_all, file = "p180_TDP_mean_HOUR_TDP_all.Rdata")

# WT save
save(p180_TDP_mean_HOUR_WT_all, file = "p180_TDP_mean_HOUR_WT_all.Rdata")




