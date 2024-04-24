library(ggplot2)
library(gridExtra)
library(dplyr)
library(bayesplot)
library(varhandle)
library(DirichletReg)
library(ggsci)
library(reshape2)
library(xtable)

load("p180_TDP_mean_HOUR_diffs_all.Rdata")

load("p180_TDP_mean_HOUR_diffs_Dirichlet_less_all_fix.Rdata")

load("p180_TDP_mean_HOUR_diffs_MVN_Splines_all.Rdata")



load("p180_TDP_mean_HOUR_TDP_all.Rdata")

load("p180_TDP_mean_HOUR_TDP_Dirichlet_less_all_fix.Rdata")

load("p180_TDP_mean_HOUR_TDP_MVN_Splines_all.Rdata")



load("p180_TDP_mean_HOUR_WT_all.Rdata")

load("p180_TDP_mean_HOUR_WT_Dirichlet_less_all_fix.Rdata")

load("p180_TDP_mean_HOUR_WT_MVN_Splines_all.Rdata")



############################################################
############################################################

# Latex Tables

############################################################
############################################################

# Fixing columns to be how I want to format it

estTDP_ZIHGDM_0 <- TDP_hourly_p180_TDPest_ZIHGDM %>% filter(Hour == 0)

estTDP_ZI_0_test <- estTDP_ZIHGDM_0

estTDP_ZI_0_test$right_brack <- " ("

estTDP_ZI_0_test$left_brack <- ")"

estTDP_ZI_0_test$comma <- ","

estTDP_ZI_0_test$Test <- paste0(estTDP_ZI_0_test$Mean, estTDP_ZI_0_test$right_brack, estTDP_ZI_0_test$Lower,
                                estTDP_ZI_0_test$comma, estTDP_ZI_0_test$Upper, estTDP_ZI_0_test$left_brack)


estTDP_ZI_12_test <- estTDP_ZI_0_test


# Do the same thing but with all hours included for a massive table
estTDP_ZIHGDM <- TDP_hourly_p180_TDPest_ZIHGDM

estTDP_ZI_test <- estTDP_ZIHGDM

estTDP_ZI_test$right_brack <- " ("

estTDP_ZI_test$left_brack <- ")"

estTDP_ZI_test$comma <- ","

estTDP_ZI_test$Test <- paste0(estTDP_ZI_test$Mean, estTDP_ZI_test$right_brack, estTDP_ZI_test$Lower,
                                estTDP_ZI_test$comma, estTDP_ZI_test$Upper, estTDP_ZI_test$left_brack)

############################################################
############################################################


# Hour 12

estTDP_ZIHGDM_12 <- TDP_hourly_p180_TDPest_ZIHGDM %>% filter(Hour == 12)
estTDP_ZIHGDM_12$CI_mean <- paste0(estTDP_ZIHGDM_12$Mean, estTDP_ZI_12_test$right_brack, estTDP_ZIHGDM_12$Lower,
                                  estTDP_ZI_12_test$comma, estTDP_ZIHGDM_12$Upper, estTDP_ZI_12_test$left_brack)

estTDP_Dirichlet_12 <- TDP_hourly_p180_TDPest_Dirichlet %>% filter(Hour == 12)
estTDP_Dirichlet_12$CI_mean <- paste0(estTDP_Dirichlet_12$Mean, estTDP_ZI_12_test$right_brack, estTDP_Dirichlet_12$Lower,
                                     estTDP_ZI_12_test$comma, estTDP_Dirichlet_12$Upper, estTDP_ZI_12_test$left_brack)

estTDP_MVNorm_12 <- TDP_hourly_p180_TDPest_MVNorm %>% filter(Hour == 12)
estTDP_MVNorm_12$CI_mean <- paste0(estTDP_MVNorm_12$Mean, estTDP_ZI_12_test$right_brack, estTDP_MVNorm_12$Lower,
                                  estTDP_ZI_12_test$comma, estTDP_MVNorm_12$Upper, estTDP_ZI_12_test$left_brack)



estWT_ZIHGDM_12 <- TDP_hourly_p180_WTest_ZIHGDM %>% filter(Hour == 12)
estWT_ZIHGDM_12$CI_mean <- paste0(estWT_ZIHGDM_12$Mean, estTDP_ZI_12_test$right_brack, estWT_ZIHGDM_12$Lower,
                                 estTDP_ZI_12_test$comma, estWT_ZIHGDM_12$Upper, estTDP_ZI_12_test$left_brack)

estWT_Dirichlet_12 <- TDP_hourly_p180_WTest_Dirichlet %>% filter(Hour == 12)
estWT_Dirichlet_12$CI_mean <- paste0(estWT_Dirichlet_12$Mean, estTDP_ZI_12_test$right_brack, estWT_Dirichlet_12$Lower,
                                    estTDP_ZI_12_test$comma, estWT_Dirichlet_12$Upper, estTDP_ZI_12_test$left_brack)

estWT_MVNorm_12 <- TDP_hourly_p180_WTest_MVNorm %>% filter(Hour == 12)
estWT_MVNorm_12$CI_mean <- paste0(estWT_MVNorm_12$Mean, estTDP_ZI_12_test$right_brack, estWT_MVNorm_12$Lower,
                                 estTDP_ZI_12_test$comma, estWT_MVNorm_12$Upper, estTDP_ZI_12_test$left_brack)



ests_12 <- cbind(estTDP_ZIHGDM_12[,c(1,2,6)], estWT_ZIHGDM_12$CI_mean, 
                estTDP_Dirichlet_12$CI_mean, estWT_Dirichlet_12$CI_mean, 
                estTDP_MVNorm_12$CI_mean, estWT_MVNorm_12$CI_mean)


ests_12_xtable <- xtable(ests_12)

print(ests_12_xtable, include.rownames = F, include.colnames = F)





# Hour 0

estTDP_ZIHGDM_0 <- TDP_hourly_p180_TDPest_ZIHGDM %>% filter(Hour == 0)
estTDP_ZIHGDM_0$CI_mean <- paste0(estTDP_ZIHGDM_0$Mean, estTDP_ZI_0_test$right_brack, estTDP_ZIHGDM_0$Lower,
                                  estTDP_ZI_0_test$comma, estTDP_ZIHGDM_0$Upper, estTDP_ZI_0_test$left_brack)

estTDP_Dirichlet_0 <- TDP_hourly_p180_TDPest_Dirichlet %>% filter(Hour == 0)
estTDP_Dirichlet_0$CI_mean <- paste0(estTDP_Dirichlet_0$Mean, estTDP_ZI_0_test$right_brack, estTDP_Dirichlet_0$Lower,
                                  estTDP_ZI_0_test$comma, estTDP_Dirichlet_0$Upper, estTDP_ZI_0_test$left_brack)

estTDP_MVNorm_0 <- TDP_hourly_p180_TDPest_MVNorm %>% filter(Hour == 0)
estTDP_MVNorm_0$CI_mean <- paste0(estTDP_MVNorm_0$Mean, estTDP_ZI_0_test$right_brack, estTDP_MVNorm_0$Lower,
                                  estTDP_ZI_0_test$comma, estTDP_MVNorm_0$Upper, estTDP_ZI_0_test$left_brack)



estWT_ZIHGDM_0 <- TDP_hourly_p180_WTest_ZIHGDM %>% filter(Hour == 0)
estWT_ZIHGDM_0$CI_mean <- paste0(estWT_ZIHGDM_0$Mean, estTDP_ZI_0_test$right_brack, estWT_ZIHGDM_0$Lower,
                                  estTDP_ZI_0_test$comma, estWT_ZIHGDM_0$Upper, estTDP_ZI_0_test$left_brack)

estWT_Dirichlet_0 <- TDP_hourly_p180_WTest_Dirichlet %>% filter(Hour == 0)
estWT_Dirichlet_0$CI_mean <- paste0(estWT_Dirichlet_0$Mean, estTDP_ZI_0_test$right_brack, estWT_Dirichlet_0$Lower,
                                 estTDP_ZI_0_test$comma, estWT_Dirichlet_0$Upper, estTDP_ZI_0_test$left_brack)

estWT_MVNorm_0 <- TDP_hourly_p180_WTest_MVNorm %>% filter(Hour == 0)
estWT_MVNorm_0$CI_mean <- paste0(estWT_MVNorm_0$Mean, estTDP_ZI_0_test$right_brack, estWT_MVNorm_0$Lower,
                                 estTDP_ZI_0_test$comma, estWT_MVNorm_0$Upper, estTDP_ZI_0_test$left_brack)



ests_0 <- cbind(estTDP_ZIHGDM_0[,c(1,2,6)], estWT_ZIHGDM_0$CI_mean, 
                 estTDP_Dirichlet_0$CI_mean, estWT_Dirichlet_0$CI_mean, 
                 estTDP_MVNorm_0$CI_mean, estWT_MVNorm_0$CI_mean)


ests_0_xtable <- xtable(ests_0)

print(ests_0_xtable, include.rownames = F, include.colnames = F)



# ALL HOURS BOTH TDP AND WT

estTDP_ZIHGDM <- TDP_hourly_p180_TDPest_ZIHGDM 
estTDP_ZIHGDM$CI_mean <- paste0(estTDP_ZIHGDM$Mean, estTDP_ZI_test$right_brack, estTDP_ZIHGDM$Lower,
                                  estTDP_ZI_test$comma, estTDP_ZIHGDM$Upper, estTDP_ZI_test$left_brack)

estTDP_Dirichlet <- TDP_hourly_p180_TDPest_Dirichlet 
estTDP_Dirichlet$CI_mean <- paste0(estTDP_Dirichlet$Mean, estTDP_ZI_test$right_brack, estTDP_Dirichlet$Lower,
                                     estTDP_ZI_test$comma, estTDP_Dirichlet$Upper, estTDP_ZI_test$left_brack)

estTDP_MVNorm <- TDP_hourly_p180_TDPest_MVNorm 
estTDP_MVNorm$CI_mean <- paste0(estTDP_MVNorm$Mean, estTDP_ZI_test$right_brack, estTDP_MVNorm$Lower,
                                  estTDP_ZI_test$comma, estTDP_MVNorm$Upper, estTDP_ZI_test$left_brack)



estWT_ZIHGDM <- TDP_hourly_p180_WTest_ZIHGDM 
estWT_ZIHGDM$CI_mean <- paste0(estWT_ZIHGDM$Mean, estTDP_ZI_test$right_brack, estWT_ZIHGDM$Lower,
                                 estTDP_ZI_test$comma, estWT_ZIHGDM$Upper, estTDP_ZI_test$left_brack)

estWT_Dirichlet <- TDP_hourly_p180_WTest_Dirichlet 
estWT_Dirichlet$CI_mean <- paste0(estWT_Dirichlet$Mean, estTDP_ZI_test$right_brack, estWT_Dirichlet$Lower,
                                    estTDP_ZI_test$comma, estWT_Dirichlet$Upper, estTDP_ZI_test$left_brack)

estWT_MVNorm <- TDP_hourly_p180_WTest_MVNorm
estWT_MVNorm$CI_mean <- paste0(estWT_MVNorm$Mean, estTDP_ZI_test$right_brack, estWT_MVNorm$Lower,
                                 estTDP_ZI_test$comma, estWT_MVNorm$Upper, estTDP_ZI_test$left_brack)



ests <- cbind(estTDP_ZIHGDM[,c(1,2,6)], estWT_ZIHGDM$CI_mean, 
                estTDP_Dirichlet$CI_mean, estWT_Dirichlet$CI_mean, 
                estTDP_MVNorm$CI_mean, estWT_MVNorm$CI_mean)


ests_xtable <- xtable(ests)

print(ests_xtable, include.rownames = F, include.colnames = F)




# DIFFERENCES ALL HOURS

diff_ZIHGDM <- TDP_hourly_p180_estimates_ZIHGDM 
diff_ZIHGDM$CI_mean <- paste0(diff_ZIHGDM$Mean, estTDP_ZI_test$right_brack, diff_ZIHGDM$Lower,
                                estTDP_ZI_test$comma, diff_ZIHGDM$Upper, estTDP_ZI_test$left_brack)

diff_Dirichlet <- TDP_hourly_p180_estimates_Dirichlet 
diff_Dirichlet$CI_mean <- paste0(diff_Dirichlet$Mean, estTDP_ZI_test$right_brack, diff_Dirichlet$Lower,
                                   estTDP_ZI_test$comma, diff_Dirichlet$Upper, estTDP_ZI_test$left_brack)

diff_MVNorm <- TDP_hourly_p180_estimates_MVNorm 
diff_MVNorm$CI_mean <- paste0(diff_MVNorm$Mean, estTDP_ZI_test$right_brack, diff_MVNorm$Lower,
                                estTDP_ZI_test$comma, diff_MVNorm$Upper, estTDP_ZI_test$left_brack)



diffs <- cbind(diff_ZIHGDM[,c(1,2,6)], 
                  diff_Dirichlet$CI_mean,
                  diff_MVNorm$CI_mean)


diffs_xtable <- xtable(diffs)

print(diffs_xtable, include.rownames = F, include.colnames = F, scalebox = 0.7) # tabular.environment = "longtable"



# Differences Hour 12


diff_ZIHGDM_12 <- TDP_hourly_p180_estimates_ZIHGDM %>% filter(Hour == 12)
diff_ZIHGDM_12$CI_mean <- paste0(diff_ZIHGDM_12$Mean, estTDP_ZI_12_test$right_brack, diff_ZIHGDM_12$Lower,
                                   estTDP_ZI_12_test$comma, diff_ZIHGDM_12$Upper, estTDP_ZI_12_test$left_brack)

diff_Dirichlet_12 <- TDP_hourly_p180_estimates_Dirichlet %>% filter(Hour == 12)
diff_Dirichlet_12$CI_mean <- paste0(diff_Dirichlet_12$Mean, estTDP_ZI_12_test$right_brack, diff_Dirichlet_12$Lower,
                                 estTDP_ZI_12_test$comma, diff_Dirichlet_12$Upper, estTDP_ZI_12_test$left_brack)

diff_MVNorm_12 <- TDP_hourly_p180_estimates_MVNorm %>% filter(Hour == 12)
diff_MVNorm_12$CI_mean <- paste0(diff_MVNorm_12$Mean, estTDP_ZI_12_test$right_brack, diff_MVNorm_12$Lower,
                                 estTDP_ZI_12_test$comma, diff_MVNorm_12$Upper, estTDP_ZI_12_test$left_brack)


diffs_12 <- cbind(diff_ZIHGDM_12[,c(1,2,6)], 
                  diff_Dirichlet_12$CI_mean,
                  diff_MVNorm_12$CI_mean)


diffs_12_xtable <- xtable(diffs_12)

print(diffs_12_xtable, include.rownames = F, include.colnames = F)


# Differences Hour 0


diff_ZIHGDM_0 <- TDP_hourly_p180_estimates_ZIHGDM %>% filter(Hour == 0)
diff_ZIHGDM_0$CI_mean <- paste0(diff_ZIHGDM_0$Mean, estTDP_ZI_0_test$right_brack, diff_ZIHGDM_0$Lower,
                                 estTDP_ZI_0_test$comma, diff_ZIHGDM_0$Upper, estTDP_ZI_0_test$left_brack)

diff_Dirichlet_0 <- TDP_hourly_p180_estimates_Dirichlet %>% filter(Hour == 0)
diff_Dirichlet_0$CI_mean <- paste0(diff_Dirichlet_0$Mean, estTDP_ZI_0_test$right_brack, diff_Dirichlet_0$Lower,
                                    estTDP_ZI_0_test$comma, diff_Dirichlet_0$Upper, estTDP_ZI_0_test$left_brack)

diff_MVNorm_0 <- TDP_hourly_p180_estimates_MVNorm %>% filter(Hour == 0)
diff_MVNorm_0$CI_mean <- paste0(diff_MVNorm_0$Mean, estTDP_ZI_0_test$right_brack, diff_MVNorm_0$Lower,
                                 estTDP_ZI_0_test$comma, diff_MVNorm_0$Upper, estTDP_ZI_0_test$left_brack)


diffs_0 <- cbind(diff_ZIHGDM_0[,c(1,2,6)], 
                  diff_Dirichlet_0$CI_mean,
                  diff_MVNorm_0$CI_mean)


diffs_0_xtable <- xtable(diffs_0)

print(diffs_0_xtable, include.rownames = F, include.colnames = F)


# Differences Hour 21


diff_ZIHGDM_21 <- TDP_hourly_p180_estimates_ZIHGDM %>% filter(Hour == 21)
diff_ZIHGDM_21$CI_mean <- paste0(diff_ZIHGDM_21$Mean, estTDP_ZI_12_test$right_brack, diff_ZIHGDM_21$Lower,
                                 estTDP_ZI_12_test$comma, diff_ZIHGDM_21$Upper, estTDP_ZI_12_test$left_brack)

diff_Dirichlet_21 <- TDP_hourly_p180_estimates_Dirichlet %>% filter(Hour == 21)
diff_Dirichlet_21$CI_mean <- paste0(diff_Dirichlet_21$Mean, estTDP_ZI_12_test$right_brack, diff_Dirichlet_21$Lower,
                                    estTDP_ZI_12_test$comma, diff_Dirichlet_21$Upper, estTDP_ZI_12_test$left_brack)

diff_MVNorm_21 <- TDP_hourly_p180_estimates_MVNorm %>% filter(Hour == 21)
diff_MVNorm_21$CI_mean <- paste0(diff_MVNorm_21$Mean, estTDP_ZI_12_test$right_brack, diff_MVNorm_21$Lower,
                                 estTDP_ZI_12_test$comma, diff_MVNorm_21$Upper, estTDP_ZI_12_test$left_brack)


diffs_21 <- cbind(diff_ZIHGDM_21[,c(1,2,6)], 
                  diff_Dirichlet_21$CI_mean,
                  diff_MVNorm_21$CI_mean)


diffs_21_xtable <- xtable(diffs_21)

print(diffs_21_xtable, include.rownames = F, include.colnames = F)

