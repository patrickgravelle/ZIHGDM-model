library(ggplot2)
library(gridExtra)
library(dplyr)
library(bayesplot)
library(varhandle)
library(DirichletReg)
library(readr)



mice_new <- read_csv("AllFallonTDP43Dec2021.csv", col_types = cols(`Day of Recording` = col_number(), 
                                                                   `Hour (0-23)` = col_number(), drink = col_number(), 
                                                                   eat = col_number(), groom = col_number(), 
                                                                   hang = col_number(), sniff = col_number(), 
                                                                   rear = col_number(), rest = col_number(), 
                                                                   walk = col_number(), eathand = col_number()))



mice_new <- read_csv("AllFallonTDP43Dec2021-4mo.csv", col_types = cols(`Day of Recording` = col_number(), 
                                                                   `Hour (0-23)` = col_number(), drink = col_number(), 
                                                                   eat = col_number(), groom = col_number(), 
                                                                   hang = col_number(), sniff = col_number(), 
                                                                   rear = col_number(), rest = col_number(), 
                                                                   walk = col_number(), eathand = col_number()))





mice_TDP43 <- mice_new %>% filter(`Day of Recording` > 0)

# should remove the observation with negative day of recording for WT_9 and WT_7
# the negative day corresponds to the number of days before the first recording took place in that batch
# so it's like the camera turned on accidentally or for a test so that data doesn't need to be there



# fix names of the columns
names(mice_TDP43) <- c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "Groom", "Hang", "Sniff", "Rear", "Rest", "Walk", "EBH")

# order the columns as in the original analysis
mice_TDP43 <- mice_TDP43[ , c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")]


# all fixed so now let's make a factor so that WT is baseline
mice_TDP43$GenotypeTDP <- ifelse(mice_TDP43$Genotype == "TDP43_Q331K", 1, 0)


mice_TDP43$AgeVec <- ifelse(mice_TDP43$Age == 4, "p120", 
                            ifelse(mice_TDP43$Age == 6, "p180", 
                                   ifelse(mice_TDP43$Age == 7.5, "p225", 
                                          ifelse(mice_TDP43$Age == 10, "p300", "p345"))))


mice_TDP43$GenderVec <- ifelse(mice_TDP43$Gender == "Male", 1, 0)

# only keeping mice aged p180 days 
mice_TDP43_p120 <- mice_TDP43 %>% filter(AgeVec == "p120")
mice_TDP43_p120$id <- match(mice_TDP43_p120$MouseID, unique(mice_TDP43_p120$MouseID))
mice_positive <- mice_TDP43_p120

mice_TDP43_p180 <- mice_TDP43 %>% filter(AgeVec == "p180")
mice_TDP43_p180$id <- match(mice_TDP43_p180$MouseID, unique(mice_TDP43_p180$MouseID))
mice_positive <- mice_TDP43_p180

mice_TDP43_p225 <- mice_TDP43 %>% filter(AgeVec == "p225")
mice_TDP43_p225$id <- match(mice_TDP43_p225$MouseID, unique(mice_TDP43_p225$MouseID))
mice_positive <- mice_TDP43_p225

mice_TDP43_p300 <- mice_TDP43 %>% filter(AgeVec == "p300")
mice_TDP43_p300$id <- match(mice_TDP43_p300$MouseID, unique(mice_TDP43_p300$MouseID))
mice_positive <- mice_TDP43_p300

mice_TDP43_p345 <- mice_TDP43 %>% filter(AgeVec == "p345")
mice_TDP43_p345$id <- match(mice_TDP43_p345$MouseID, unique(mice_TDP43_p345$MouseID))
mice_positive <- mice_TDP43_p345









smart_round <- function(x) {
  x <- as.matrix(x)
  y <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:nrow(x)){
    y[i,] <- floor(x[i,])
    indices <- tail(order(x[i,]-y[i,]), round(sum(x[i,])) - sum(y[i,]))
    y[i, indices] <- y[i, indices] + 1
  }
  y
  
}




# create a TDP vector to identify each mouse's observations
mice_tdp <- mice_positive$GenotypeTDP

table(mice_tdp)

# create a vector with how many observations we have for each specific TDP
tdp_reps <- c()

for(i in 1:2){
  tdp_reps <- c(tdp_reps, table(mice_tdp)[[i]])
}

sum(tdp_reps)

# create binary variables for each TDP identifying when they were observed
bnry_mice_tdp <- to.dummy(mice_tdp, "TDP")

# replicate and rbind the matrix above by 9 (one per behaviour)
bnry9times_tdp <- do.call(rbind, replicate(9, bnry_mice_tdp, simplify = F))






####################################################
####################################################
# Only run for Dirichlet Models

load("Dirichlet_adjusted_Xdat_Splines_Interaction_REs_less_p180.Rdata")

yreps_log2 <- Dirichlet_adjusted_Xdat_Splines_Interaction_REs_less_p180[,c(43588,65169)]
yreps_log2 <- Dirichlet_adjusted_Xdat_Splines_Interaction_REs_less_p180[,c(43588:65169)]


load("Dirichlet_adjusted_Xdat_Splines_Interaction_REs_less_LogLog_p180.Rdata")

yreps_log2 <- Dirichlet_adjusted_Xdat_Splines_Interaction_REs_less_LogLog_p180[,c(43588,65169)]
yreps_log2 <- Dirichlet_adjusted_Xdat_Splines_Interaction_REs_less_LogLog_p180[,c(43588:65169)]


# Need this for Dirichlet Models
mice_positive$Drink <- mice_positive$Drink + 10
mice_positive$Eat <- mice_positive$Eat + 10
mice_positive$EBH <- mice_positive$EBH + 10
mice_positive$Groom <- mice_positive$Groom + 10
mice_positive$Hang <- mice_positive$Hang + 10
mice_positive$Rear <- mice_positive$Rear + 10
mice_positive$Rest <- mice_positive$Rest + 10
mice_positive$Sniff <- mice_positive$Sniff + 10
mice_positive$Walk <- mice_positive$Walk + 10

mice_positive$Y <- DR_data(mice_positive[,8:16])
####################################################
####################################################



# ZIHGDM final model - stats paper
load("ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus_sml.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus_sml[,c(213727,235308)]
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus_sml[,c(213727:235308)]

# 
# 
# # ZIHGDM final model 
# load("ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p120_taus_sml.Rdata")
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p120_taus_sml[,c(213840,235430)]
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p120_taus_sml[,c(213840:235430)]
# 
# # ZIHGDM final model 
# load("ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_taus_sml.Rdata")
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_taus_sml[,c(213751,235332)]
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_taus_sml[,c(213751:235332)]
# 
# 
# # ZIHGDM final model 
# load("ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p225_taus_sml.Rdata")
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p225_taus_sml[,c(213751,235332)]
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p225_taus_sml[,c(213751:235332)]
# 
# # ZIHGDM final model 
# load("ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p300_taus_sml.Rdata")
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p300_taus_sml[,c(213840,235430)]
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p300_taus_sml[,c(213840:235430)]
# 
# # ZIHGDM final model 
# load("ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p345_taus_sml.Rdata")
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p345_taus_sml[,c(213306,234842)]
# yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p345_taus_sml[,c(213306:234842)]



# Need this for Multinomial Models
outcomes_int <- smart_round(mice_positive[,8:16])

mice_positive$Y <- outcomes_int







# matrix of each mouse's observations

hour <- rep((mice_positive$Hour),9)

yreps_mat <- rbind(yreps_log2, hour)


yreps_mat <- data.matrix(yreps_mat)







# Matrices with mice having only TDP or no TDP
TDP_no <- t(yreps_mat[ , as.logical(bnry9times_tdp[,1]), drop = F])
TDP_yes <- t(yreps_mat[ , as.logical(bnry9times_tdp[,2]), drop = F])






#########################
# p120
#########################

## TDP NO
posterior_pred1_dat1 <- data.frame(TDP_no[1:1200,])
posterior_pred2_dat1 <- data.frame(TDP_no[(1 + 1200*1):(1200 + 1200*1),])
posterior_pred3_dat1 <- data.frame(TDP_no[(1 + 1200*2):(1200 + 1200*2),])
posterior_pred4_dat1 <- data.frame(TDP_no[(1 + 1200*3):(1200 + 1200*3),])
posterior_pred5_dat1 <- data.frame(TDP_no[(1 + 1200*4):(1200 + 1200*4),])
posterior_pred6_dat1 <- data.frame(TDP_no[(1 + 1200*5):(1200 + 1200*5),])
posterior_pred7_dat1 <- data.frame(TDP_no[(1 + 1200*6):(1200 + 1200*6),])
posterior_pred8_dat1 <- data.frame(TDP_no[(1 + 1200*7):(1200 + 1200*7),])
posterior_pred9_dat1 <- data.frame(TDP_no[(1 + 1200*8):(1200 + 1200*8),])

mice_positive_sub_no <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(TDP_yes[1:1199,])
posterior_pred2_dat2 <- data.frame(TDP_yes[(1 + 1199*1):(1199 + 1199*1),])
posterior_pred3_dat2 <- data.frame(TDP_yes[(1 + 1199*2):(1199 + 1199*2),])
posterior_pred4_dat2 <- data.frame(TDP_yes[(1 + 1199*3):(1199 + 1199*3),])
posterior_pred5_dat2 <- data.frame(TDP_yes[(1 + 1199*4):(1199 + 1199*4),])
posterior_pred6_dat2 <- data.frame(TDP_yes[(1 + 1199*5):(1199 + 1199*5),])
posterior_pred7_dat2 <- data.frame(TDP_yes[(1 + 1199*6):(1199 + 1199*6),])
posterior_pred8_dat2 <- data.frame(TDP_yes[(1 + 1199*7):(1199 + 1199*7),])
posterior_pred9_dat2 <- data.frame(TDP_yes[(1 + 1199*8):(1199 + 1199*8),])

mice_positive_sub_yes <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]



#########################
# p180
#########################

## TDP NO
posterior_pred1_dat1 <- data.frame(TDP_no[1:1199,])
posterior_pred2_dat1 <- data.frame(TDP_no[(1 + 1199*1):(1199 + 1199*1),])
posterior_pred3_dat1 <- data.frame(TDP_no[(1 + 1199*2):(1199 + 1199*2),])
posterior_pred4_dat1 <- data.frame(TDP_no[(1 + 1199*3):(1199 + 1199*3),])
posterior_pred5_dat1 <- data.frame(TDP_no[(1 + 1199*4):(1199 + 1199*4),])
posterior_pred6_dat1 <- data.frame(TDP_no[(1 + 1199*5):(1199 + 1199*5),])
posterior_pred7_dat1 <- data.frame(TDP_no[(1 + 1199*6):(1199 + 1199*6),])
posterior_pred8_dat1 <- data.frame(TDP_no[(1 + 1199*7):(1199 + 1199*7),])
posterior_pred9_dat1 <- data.frame(TDP_no[(1 + 1199*8):(1199 + 1199*8),])

mice_positive_sub_no <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(TDP_yes[1:1199,])
posterior_pred2_dat2 <- data.frame(TDP_yes[(1 + 1199*1):(1199 + 1199*1),])
posterior_pred3_dat2 <- data.frame(TDP_yes[(1 + 1199*2):(1199 + 1199*2),])
posterior_pred4_dat2 <- data.frame(TDP_yes[(1 + 1199*3):(1199 + 1199*3),])
posterior_pred5_dat2 <- data.frame(TDP_yes[(1 + 1199*4):(1199 + 1199*4),])
posterior_pred6_dat2 <- data.frame(TDP_yes[(1 + 1199*5):(1199 + 1199*5),])
posterior_pred7_dat2 <- data.frame(TDP_yes[(1 + 1199*6):(1199 + 1199*6),])
posterior_pred8_dat2 <- data.frame(TDP_yes[(1 + 1199*7):(1199 + 1199*7),])
posterior_pred9_dat2 <- data.frame(TDP_yes[(1 + 1199*8):(1199 + 1199*8),])

mice_positive_sub_yes <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]



#########################
# p225
#########################

## TDP NO
posterior_pred1_dat1 <- data.frame(TDP_no[1:1199,])
posterior_pred2_dat1 <- data.frame(TDP_no[(1 + 1199*1):(1199 + 1199*1),])
posterior_pred3_dat1 <- data.frame(TDP_no[(1 + 1199*2):(1199 + 1199*2),])
posterior_pred4_dat1 <- data.frame(TDP_no[(1 + 1199*3):(1199 + 1199*3),])
posterior_pred5_dat1 <- data.frame(TDP_no[(1 + 1199*4):(1199 + 1199*4),])
posterior_pred6_dat1 <- data.frame(TDP_no[(1 + 1199*5):(1199 + 1199*5),])
posterior_pred7_dat1 <- data.frame(TDP_no[(1 + 1199*6):(1199 + 1199*6),])
posterior_pred8_dat1 <- data.frame(TDP_no[(1 + 1199*7):(1199 + 1199*7),])
posterior_pred9_dat1 <- data.frame(TDP_no[(1 + 1199*8):(1199 + 1199*8),])

mice_positive_sub_no <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(TDP_yes[1:1199,])
posterior_pred2_dat2 <- data.frame(TDP_yes[(1 + 1199*1):(1199 + 1199*1),])
posterior_pred3_dat2 <- data.frame(TDP_yes[(1 + 1199*2):(1199 + 1199*2),])
posterior_pred4_dat2 <- data.frame(TDP_yes[(1 + 1199*3):(1199 + 1199*3),])
posterior_pred5_dat2 <- data.frame(TDP_yes[(1 + 1199*4):(1199 + 1199*4),])
posterior_pred6_dat2 <- data.frame(TDP_yes[(1 + 1199*5):(1199 + 1199*5),])
posterior_pred7_dat2 <- data.frame(TDP_yes[(1 + 1199*6):(1199 + 1199*6),])
posterior_pred8_dat2 <- data.frame(TDP_yes[(1 + 1199*7):(1199 + 1199*7),])
posterior_pred9_dat2 <- data.frame(TDP_yes[(1 + 1199*8):(1199 + 1199*8),])

mice_positive_sub_yes <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]


#########################
# p300
#########################

## TDP NO
posterior_pred1_dat1 <- data.frame(TDP_no[1:1199,])
posterior_pred2_dat1 <- data.frame(TDP_no[(1 + 1199*1):(1199 + 1199*1),])
posterior_pred3_dat1 <- data.frame(TDP_no[(1 + 1199*2):(1199 + 1199*2),])
posterior_pred4_dat1 <- data.frame(TDP_no[(1 + 1199*3):(1199 + 1199*3),])
posterior_pred5_dat1 <- data.frame(TDP_no[(1 + 1199*4):(1199 + 1199*4),])
posterior_pred6_dat1 <- data.frame(TDP_no[(1 + 1199*5):(1199 + 1199*5),])
posterior_pred7_dat1 <- data.frame(TDP_no[(1 + 1199*6):(1199 + 1199*6),])
posterior_pred8_dat1 <- data.frame(TDP_no[(1 + 1199*7):(1199 + 1199*7),])
posterior_pred9_dat1 <- data.frame(TDP_no[(1 + 1199*8):(1199 + 1199*8),])

mice_positive_sub_no <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(TDP_yes[1:1200,])
posterior_pred2_dat2 <- data.frame(TDP_yes[(1 + 1200*1):(1200 + 1200*1),])
posterior_pred3_dat2 <- data.frame(TDP_yes[(1 + 1200*2):(1200 + 1200*2),])
posterior_pred4_dat2 <- data.frame(TDP_yes[(1 + 1200*3):(1200 + 1200*3),])
posterior_pred5_dat2 <- data.frame(TDP_yes[(1 + 1200*4):(1200 + 1200*4),])
posterior_pred6_dat2 <- data.frame(TDP_yes[(1 + 1200*5):(1200 + 1200*5),])
posterior_pred7_dat2 <- data.frame(TDP_yes[(1 + 1200*6):(1200 + 1200*6),])
posterior_pred8_dat2 <- data.frame(TDP_yes[(1 + 1200*7):(1200 + 1200*7),])
posterior_pred9_dat2 <- data.frame(TDP_yes[(1 + 1200*8):(1200 + 1200*8),])

mice_positive_sub_yes <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]



#########################
# p345
#########################

## TDP NO
posterior_pred1_dat1 <- data.frame(TDP_no[1:1197,])
posterior_pred2_dat1 <- data.frame(TDP_no[(1 + 1197*1):(1197 + 1197*1),])
posterior_pred3_dat1 <- data.frame(TDP_no[(1 + 1197*2):(1197 + 1197*2),])
posterior_pred4_dat1 <- data.frame(TDP_no[(1 + 1197*3):(1197 + 1197*3),])
posterior_pred5_dat1 <- data.frame(TDP_no[(1 + 1197*4):(1197 + 1197*4),])
posterior_pred6_dat1 <- data.frame(TDP_no[(1 + 1197*5):(1197 + 1197*5),])
posterior_pred7_dat1 <- data.frame(TDP_no[(1 + 1197*6):(1197 + 1197*6),])
posterior_pred8_dat1 <- data.frame(TDP_no[(1 + 1197*7):(1197 + 1197*7),])
posterior_pred9_dat1 <- data.frame(TDP_no[(1 + 1197*8):(1197 + 1197*8),])

mice_positive_sub_no <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(TDP_yes[1:1196,])
posterior_pred2_dat2 <- data.frame(TDP_yes[(1 + 1196*1):(1196 + 1196*1),])
posterior_pred3_dat2 <- data.frame(TDP_yes[(1 + 1196*2):(1196 + 1196*2),])
posterior_pred4_dat2 <- data.frame(TDP_yes[(1 + 1196*3):(1196 + 1196*3),])
posterior_pred5_dat2 <- data.frame(TDP_yes[(1 + 1196*4):(1196 + 1196*4),])
posterior_pred6_dat2 <- data.frame(TDP_yes[(1 + 1196*5):(1196 + 1196*5),])
posterior_pred7_dat2 <- data.frame(TDP_yes[(1 + 1196*6):(1196 + 1196*6),])
posterior_pred8_dat2 <- data.frame(TDP_yes[(1 + 1196*7):(1196 + 1196*7),])
posterior_pred9_dat2 <- data.frame(TDP_yes[(1 + 1196*8):(1196 + 1196*8),])

mice_positive_sub_yes <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]



################################################################################################################ 
################################################################################################################ 
################################################################################################################ 
################################################################################################################ 

# MEDIAN

################################################################################################################ 
################################################################################################################ 
################################################################################################################ 
################################################################################################################ 
my_colors <- c("#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")
############################ Drink ############################
obs_pred1 <- aggregate(x = mice_positive_sub_no$Y[,1:9], by = list(mice_positive_sub_no$Hour), FUN = median)
obs_pred2 <- aggregate(x = mice_positive_sub_yes$Y[,1:9], by = list(mice_positive_sub_yes$Hour), FUN = median)

names(obs_pred1) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")
names(obs_pred2) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")



# # Dirichlet Models Only
# posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:3500], by = list(posterior_pred1_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:3500], by = list(posterior_pred1_dat2$X3501), FUN = median)

posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:2500], by = list(posterior_pred1_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:2500], by = list(posterior_pred1_dat2$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p1 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() +  
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Drink Median Behavior Time") + xlab("Drink Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs")) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)

# checkpoint




############################ Eat ############################
posterior_pred11 <- aggregate(x = posterior_pred2_dat1[,1:2500], by = list(posterior_pred2_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred2_dat2[,1:2500], by = list(posterior_pred2_dat2$X25011), FUN = median)

# # Dirichlet Models Only
# posterior_pred11 <- aggregate(x = posterior_pred2_dat1[,1:3500], by = list(posterior_pred2_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred2_dat2[,1:3500], by = list(posterior_pred2_dat2$X3501), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Eat Median Behavior Time") + xlab("Eat Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs")) +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs")) +
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs"), alpha=0.75, size = 0.5) +
  scale_colour_manual(values = my_colors)


############################ EBH ############################

posterior_pred11 <- aggregate(x = posterior_pred3_dat1[,1:2500], by = list(posterior_pred3_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred3_dat2[,1:2500], by = list(posterior_pred3_dat2$X25011), FUN = median)

# # Dirichlet Models Only
# posterior_pred11 <- aggregate(x = posterior_pred3_dat1[,1:3500], by = list(posterior_pred3_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred3_dat2[,1:3500], by = list(posterior_pred3_dat2$X3501), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("EBH Median Behavior Time") + xlab("EBH Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs")) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Groom ############################

posterior_pred11 <- aggregate(x = posterior_pred4_dat1[,1:2500], by = list(posterior_pred4_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred4_dat2[,1:2500], by = list(posterior_pred4_dat2$X25011), FUN = median)

# # Dirichlet Models Only
# posterior_pred11 <- aggregate(x = posterior_pred4_dat1[,1:3500], by = list(posterior_pred4_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred4_dat2[,1:3500], by = list(posterior_pred4_dat2$X3501), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Groom Median Behavior Time") + xlab("Groom Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Hang ############################

posterior_pred11 <- aggregate(x = posterior_pred5_dat1[,1:2500], by = list(posterior_pred5_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred5_dat2[,1:2500], by = list(posterior_pred5_dat2$X25011), FUN = median)

# # Dirichlet Models Only
# posterior_pred11 <- aggregate(x = posterior_pred5_dat1[,1:3500], by = list(posterior_pred5_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred5_dat2[,1:3500], by = list(posterior_pred5_dat2$X3501), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Hang Median Behavior Time") + xlab("Hang Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rear ############################
posterior_pred11 <- aggregate(x = posterior_pred6_dat1[,1:2500], by = list(posterior_pred6_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred6_dat2[,1:2500], by = list(posterior_pred6_dat2$X25011), FUN = median)

# # Dirichlet Models Only
# posterior_pred11 <- aggregate(x = posterior_pred6_dat1[,1:3500], by = list(posterior_pred6_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred6_dat2[,1:3500], by = list(posterior_pred6_dat2$X3501), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Rear Median Behavior Time") + xlab("Rear Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rest ############################
posterior_pred11 <- aggregate(x = posterior_pred7_dat1[,1:2500], by = list(posterior_pred7_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred7_dat2[,1:2500], by = list(posterior_pred7_dat2$X25011), FUN = median)

# # Dirichlet Models Only
# posterior_pred11 <- aggregate(x = posterior_pred7_dat1[,1:3500], by = list(posterior_pred7_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred7_dat2[,1:3500], by = list(posterior_pred7_dat2$X3501), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Rest Median Behavior Time") + xlab("Rest Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)




############################ Sniff ############################
posterior_pred11 <- aggregate(x = posterior_pred8_dat1[,1:2500], by = list(posterior_pred8_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred8_dat2[,1:2500], by = list(posterior_pred8_dat2$X25011), FUN = median)

# # Dirichlet Models Only
# posterior_pred11 <- aggregate(x = posterior_pred8_dat1[,1:3500], by = list(posterior_pred8_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred8_dat2[,1:3500], by = list(posterior_pred8_dat2$X3501), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Sniff Median Behavior Time") + xlab("Sniff Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Walk ############################
posterior_pred11 <- aggregate(x = posterior_pred9_dat1[,1:2500], by = list(posterior_pred9_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred9_dat2[,1:2500], by = list(posterior_pred9_dat2$X25011), FUN = median)

# # Dirichlet Models Only
# posterior_pred11 <- aggregate(x = posterior_pred9_dat1[,1:3500], by = list(posterior_pred9_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred9_dat2[,1:3500], by = list(posterior_pred9_dat2$X3501), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Walk Median Behavior Time") + xlab("Walk Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("TDP - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Walk Median Behavior Time") +  xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



library(ggpubr)
# get a single legend at the bottom
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p9)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), 
                         p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), 
                         p5 + theme(legend.position="none"), p6 + theme(legend.position="none"),
                         p7 + theme(legend.position="none"), p8 + theme(legend.position="none"), 
                         p9 + theme(legend.position="none"),
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1), top=text_grob("TDP-43 vs. WT Posterior Predictive Median Behavior Time"))


block <- c()



###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

# MAX

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################



############################ Drink ############################
obs_pred1 <- aggregate(x = mice_positive_sub_no$Y[,1:9], by = list(mice_positive_sub_no$Hour), FUN = max)
obs_pred2 <- aggregate(x = mice_positive_sub_yes$Y[,1:9], by = list(mice_positive_sub_yes$Hour), FUN = max)

names(obs_pred1) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")
names(obs_pred2) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")


posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:2500], by = list(posterior_pred1_dat1$X25011), FUN = max)
posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:2500], by = list(posterior_pred1_dat2$X25011), FUN = max)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p1 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() +  
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Drink Max") + xlab("Drink Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Eat ############################
posterior_pred11 <- aggregate(x = posterior_pred2_dat1[,1:2500], by = list(posterior_pred2_dat1$X25011), FUN = max)
posterior_pred12 <- aggregate(x = posterior_pred2_dat2[,1:2500], by = list(posterior_pred2_dat2$X25011), FUN = max)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Eat Max") + xlab("Eat Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)


############################ EBH ############################

posterior_pred11 <- aggregate(x = posterior_pred3_dat1[,1:2500], by = list(posterior_pred3_dat1$X25011), FUN = max)
posterior_pred12 <- aggregate(x = posterior_pred3_dat2[,1:2500], by = list(posterior_pred3_dat2$X25011), FUN = max)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("EBH Max") + xlab("EBH Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Groom ############################

posterior_pred11 <- aggregate(x = posterior_pred4_dat1[,1:2500], by = list(posterior_pred4_dat1$X25011), FUN = max)
posterior_pred12 <- aggregate(x = posterior_pred4_dat2[,1:2500], by = list(posterior_pred4_dat2$X25011), FUN = max)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Groom Max") + xlab("Groom Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Hang ############################

posterior_pred11 <- aggregate(x = posterior_pred5_dat1[,1:2500], by = list(posterior_pred5_dat1$X25011), FUN = max)
posterior_pred12 <- aggregate(x = posterior_pred5_dat2[,1:2500], by = list(posterior_pred5_dat2$X25011), FUN = max)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Hang Max") + xlab("Hang Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rear ############################
posterior_pred11 <- aggregate(x = posterior_pred6_dat1[,1:2500], by = list(posterior_pred6_dat1$X25011), FUN = max)
posterior_pred12 <- aggregate(x = posterior_pred6_dat2[,1:2500], by = list(posterior_pred6_dat2$X25011), FUN = max)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Rear Max") + xlab("Rear Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rest ############################
posterior_pred11 <- aggregate(x = posterior_pred7_dat1[,1:2500], by = list(posterior_pred7_dat1$X25011), FUN = max)
posterior_pred12 <- aggregate(x = posterior_pred7_dat2[,1:2500], by = list(posterior_pred7_dat2$X25011), FUN = max)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Rest Max") + xlab("Rest Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)




############################ Sniff ############################
posterior_pred11 <- aggregate(x = posterior_pred8_dat1[,1:2500], by = list(posterior_pred8_dat1$X25011), FUN = max)
posterior_pred12 <- aggregate(x = posterior_pred8_dat2[,1:2500], by = list(posterior_pred8_dat2$X25011), FUN = max)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Sniff Max") + xlab("Sniff Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Walk ############################
posterior_pred11 <- aggregate(x = posterior_pred9_dat1[,1:2500], by = list(posterior_pred9_dat1$X25011), FUN = max)
posterior_pred12 <- aggregate(x = posterior_pred9_dat2[,1:2500], by = list(posterior_pred9_dat2$X25011), FUN = max)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk Max") + xlab("Walk Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("TDP - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk Max") + xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)


mylegend<-g_legend(p9)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), 
                         p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), 
                         p5 + theme(legend.position="none"), p6 + theme(legend.position="none"),
                         p7 + theme(legend.position="none"), p8 + theme(legend.position="none"), 
                         p9 + theme(legend.position="none"),
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1))



block <- c()



###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

# MINIMUM

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################



############################ Drink ############################
# obs_pred1 <- aggregate(x = mice_positive_sub_no[,8:16], by = list(mice_positive_sub_no$Hour), FUN = min)
# obs_pred2 <- aggregate(x = mice_positive_sub_yes[,8:16], by = list(mice_positive_sub_yes$Hour), FUN = min)

obs_pred1 <- aggregate(x = mice_positive_sub_no$Y[,1:9], by = list(mice_positive_sub_no$Hour), FUN = min)
obs_pred2 <- aggregate(x = mice_positive_sub_yes$Y[,1:9], by = list(mice_positive_sub_yes$Hour), FUN = min)

names(obs_pred1) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")
names(obs_pred2) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")


posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:2500], by = list(posterior_pred1_dat1$X25011), FUN = min)
posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:2500], by = list(posterior_pred1_dat2$X25011), FUN = min)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p1 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() +  
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Drink Min") + xlab("Drink Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Eat ############################
posterior_pred11 <- aggregate(x = posterior_pred2_dat1[,1:2500], by = list(posterior_pred2_dat1$X25011), FUN = min)
posterior_pred12 <- aggregate(x = posterior_pred2_dat2[,1:2500], by = list(posterior_pred2_dat2$X25011), FUN = min)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Eat Min") + xlab("Eat Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)


############################ EBH ############################

posterior_pred11 <- aggregate(x = posterior_pred3_dat1[,1:2500], by = list(posterior_pred3_dat1$X25011), FUN = min)
posterior_pred12 <- aggregate(x = posterior_pred3_dat2[,1:2500], by = list(posterior_pred3_dat2$X25011), FUN = min)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("EBH Min") + xlab("EBH Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Groom ############################

posterior_pred11 <- aggregate(x = posterior_pred4_dat1[,1:2500], by = list(posterior_pred4_dat1$X25011), FUN = min)
posterior_pred12 <- aggregate(x = posterior_pred4_dat2[,1:2500], by = list(posterior_pred4_dat2$X25011), FUN = min)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Groom Min") + xlab("Groom Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Hang ############################

posterior_pred11 <- aggregate(x = posterior_pred5_dat1[,1:2500], by = list(posterior_pred5_dat1$X25011), FUN = min)
posterior_pred12 <- aggregate(x = posterior_pred5_dat2[,1:2500], by = list(posterior_pred5_dat2$X25011), FUN = min)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Hang Min") + xlab("Hang Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rear ############################
posterior_pred11 <- aggregate(x = posterior_pred6_dat1[,1:2500], by = list(posterior_pred6_dat1$X25011), FUN = min)
posterior_pred12 <- aggregate(x = posterior_pred6_dat2[,1:2500], by = list(posterior_pred6_dat2$X25011), FUN = min)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Rear Min") + xlab("Rear Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs")) +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs"), alpha=0.75, size = 0.5) +
  scale_colour_manual(values = my_colors)



############################ Rest ############################
posterior_pred11 <- aggregate(x = posterior_pred7_dat1[,1:2500], by = list(posterior_pred7_dat1$X25011), FUN = min)
posterior_pred12 <- aggregate(x = posterior_pred7_dat2[,1:2500], by = list(posterior_pred7_dat2$X25011), FUN = min)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Rest Min") + xlab("Rest Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)




############################ Sniff ############################
posterior_pred11 <- aggregate(x = posterior_pred8_dat1[,1:2500], by = list(posterior_pred8_dat1$X25011), FUN = min)
posterior_pred12 <- aggregate(x = posterior_pred8_dat2[,1:2500], by = list(posterior_pred8_dat2$X25011), FUN = min)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Sniff Min") + xlab("Sniff Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Walk ############################
posterior_pred11 <- aggregate(x = posterior_pred9_dat1[,1:2500], by = list(posterior_pred9_dat1$X25011), FUN = min)
posterior_pred12 <- aggregate(x = posterior_pred9_dat2[,1:2500], by = list(posterior_pred9_dat2$X25011), FUN = min)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk Min") + xlab("Walk Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("TDP - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk Min") +  xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)


mylegend<-g_legend(p9)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), 
                         p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), 
                         p5 + theme(legend.position="none"), p6 + theme(legend.position="none"),
                         p7 + theme(legend.position="none"), p8 + theme(legend.position="none"), 
                         p9 + theme(legend.position="none"),
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1))



block <- c()



###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

# 25th PERCENTILE

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################



############################ Drink ############################
obs_pred1 <- aggregate(x = mice_positive_sub_no$Y[,1:9], by = list(mice_positive_sub_no$Hour), 
                       FUN = function(i) quantile(i, probs = 0.25))
obs_pred2 <- aggregate(x = mice_positive_sub_yes$Y[,1:9], by = list(mice_positive_sub_yes$Hour), 
                       FUN = function(i) quantile(i, probs = 0.25))

names(obs_pred1) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")
names(obs_pred2) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")


posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:2500], by = list(posterior_pred1_dat1$X25011), FUN = function(i) quantile(i, probs = 0.25))
posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:2500], by = list(posterior_pred1_dat2$X25011), FUN = function(i) quantile(i, probs = 0.25))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p1 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() +  
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Drink 25th Percentile") + xlab("Drink Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Eat ############################
posterior_pred11 <- aggregate(x = posterior_pred2_dat1[,1:2500], by = list(posterior_pred2_dat1$X25011), FUN = function(i) quantile(i, probs = 0.25))
posterior_pred12 <- aggregate(x = posterior_pred2_dat2[,1:2500], by = list(posterior_pred2_dat2$X25011), FUN = function(i) quantile(i, probs = 0.25))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Eat 25th Percentile") + xlab("Eat Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)


############################ EBH ############################

posterior_pred11 <- aggregate(x = posterior_pred3_dat1[,1:2500], by = list(posterior_pred3_dat1$X25011), FUN = function(i) quantile(i, probs = 0.25))
posterior_pred12 <- aggregate(x = posterior_pred3_dat2[,1:2500], by = list(posterior_pred3_dat2$X25011), FUN = function(i) quantile(i, probs = 0.25))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("EBH 25th Percentile") + xlab("EBH Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Groom ############################

posterior_pred11 <- aggregate(x = posterior_pred4_dat1[,1:2500], by = list(posterior_pred4_dat1$X25011), FUN = function(i) quantile(i, probs = 0.25))
posterior_pred12 <- aggregate(x = posterior_pred4_dat2[,1:2500], by = list(posterior_pred4_dat2$X25011), FUN = function(i) quantile(i, probs = 0.25))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Groom 25th Percentile") + xlab("Groom Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Hang ############################

posterior_pred11 <- aggregate(x = posterior_pred5_dat1[,1:2500], by = list(posterior_pred5_dat1$X25011), FUN = function(i) quantile(i, probs = 0.25))
posterior_pred12 <- aggregate(x = posterior_pred5_dat2[,1:2500], by = list(posterior_pred5_dat2$X25011), FUN = function(i) quantile(i, probs = 0.25))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Hang 25th Percentile") + xlab("Hang Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs")) +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs"), alpha=0.75, size = 0.5) +
  scale_colour_manual(values = my_colors)



############################ Rear ############################
posterior_pred11 <- aggregate(x = posterior_pred6_dat1[,1:2500], by = list(posterior_pred6_dat1$X25011), FUN = function(i) quantile(i, probs = 0.25))
posterior_pred12 <- aggregate(x = posterior_pred6_dat2[,1:2500], by = list(posterior_pred6_dat2$X25011), FUN = function(i) quantile(i, probs = 0.25))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Rear 25th Percentile") + xlab("Rear Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rest ############################
posterior_pred11 <- aggregate(x = posterior_pred7_dat1[,1:2500], by = list(posterior_pred7_dat1$X25011), FUN = function(i) quantile(i, probs = 0.25))
posterior_pred12 <- aggregate(x = posterior_pred7_dat2[,1:2500], by = list(posterior_pred7_dat2$X25011), FUN = function(i) quantile(i, probs = 0.25))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Rest 25th Percentile") + xlab("Rest Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)




############################ Sniff ############################
posterior_pred11 <- aggregate(x = posterior_pred8_dat1[,1:2500], by = list(posterior_pred8_dat1$X25011), FUN = function(i) quantile(i, probs = 0.25))
posterior_pred12 <- aggregate(x = posterior_pred8_dat2[,1:2500], by = list(posterior_pred8_dat2$X25011), FUN = function(i) quantile(i, probs = 0.25))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Sniff 25th Percentile") + xlab("Sniff Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Walk ############################
posterior_pred11 <- aggregate(x = posterior_pred9_dat1[,1:2500], by = list(posterior_pred9_dat1$X25011), FUN = function(i) quantile(i, probs = 0.25))
posterior_pred12 <- aggregate(x = posterior_pred9_dat2[,1:2500], by = list(posterior_pred9_dat2$X25011), FUN = function(i) quantile(i, probs = 0.25))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk 25th Percentile") + xlab("Walk Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("TDP - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk 25th Percentile") + xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)


mylegend<-g_legend(p9)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), 
                         p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), 
                         p5 + theme(legend.position="none"), p6 + theme(legend.position="none"),
                         p7 + theme(legend.position="none"), p8 + theme(legend.position="none"), 
                         p9 + theme(legend.position="none"),
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1))



block <- c()





###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

# 75th PERCENTILE

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################



############################ Drink ############################
obs_pred1 <- aggregate(x = mice_positive_sub_no$Y[,1:9], by = list(mice_positive_sub_no$Hour), 
                       FUN = function(i) quantile(i, probs = 0.75))
obs_pred2 <- aggregate(x = mice_positive_sub_yes$Y[,1:9], by = list(mice_positive_sub_yes$Hour), 
                       FUN = function(i) quantile(i, probs = 0.75))

names(obs_pred1) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")
names(obs_pred2) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")


posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:2500], by = list(posterior_pred1_dat1$X25011), FUN = function(i) quantile(i, probs = 0.75))
posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:2500], by = list(posterior_pred1_dat2$X25011), FUN = function(i) quantile(i, probs = 0.75))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p1 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() +  
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Drink 75th Percentile") + xlab("Drink Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Eat ############################
posterior_pred11 <- aggregate(x = posterior_pred2_dat1[,1:2500], by = list(posterior_pred2_dat1$X25011), FUN = function(i) quantile(i, probs = 0.75))
posterior_pred12 <- aggregate(x = posterior_pred2_dat2[,1:2500], by = list(posterior_pred2_dat2$X25011), FUN = function(i) quantile(i, probs = 0.75))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Eat 75th Percentile") + xlab("Eat Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)


############################ EBH ############################

posterior_pred11 <- aggregate(x = posterior_pred3_dat1[,1:2500], by = list(posterior_pred3_dat1$X25011), FUN = function(i) quantile(i, probs = 0.75))
posterior_pred12 <- aggregate(x = posterior_pred3_dat2[,1:2500], by = list(posterior_pred3_dat2$X25011), FUN = function(i) quantile(i, probs = 0.75))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("EBH 75th Percentile") + xlab("EBH Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Groom ############################

posterior_pred11 <- aggregate(x = posterior_pred4_dat1[,1:2500], by = list(posterior_pred4_dat1$X25011), FUN = function(i) quantile(i, probs = 0.75))
posterior_pred12 <- aggregate(x = posterior_pred4_dat2[,1:2500], by = list(posterior_pred4_dat2$X25011), FUN = function(i) quantile(i, probs = 0.75))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Groom 75th Percentile") + xlab("Groom Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Hang ############################

posterior_pred11 <- aggregate(x = posterior_pred5_dat1[,1:2500], by = list(posterior_pred5_dat1$X25011), FUN = function(i) quantile(i, probs = 0.75))
posterior_pred12 <- aggregate(x = posterior_pred5_dat2[,1:2500], by = list(posterior_pred5_dat2$X25011), FUN = function(i) quantile(i, probs = 0.75))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Hang 75th Percentile") + xlab("Hang Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rear ############################
posterior_pred11 <- aggregate(x = posterior_pred6_dat1[,1:2500], by = list(posterior_pred6_dat1$X25011), FUN = function(i) quantile(i, probs = 0.75))
posterior_pred12 <- aggregate(x = posterior_pred6_dat2[,1:2500], by = list(posterior_pred6_dat2$X25011), FUN = function(i) quantile(i, probs = 0.75))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Rear 75th Percentile") + xlab("Rear Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rest ############################
posterior_pred11 <- aggregate(x = posterior_pred7_dat1[,1:2500], by = list(posterior_pred7_dat1$X25011), FUN = function(i) quantile(i, probs = 0.75))
posterior_pred12 <- aggregate(x = posterior_pred7_dat2[,1:2500], by = list(posterior_pred7_dat2$X25011), FUN = function(i) quantile(i, probs = 0.75))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Rest 75th Percentile") + xlab("Rest Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)




############################ Sniff ############################
posterior_pred11 <- aggregate(x = posterior_pred8_dat1[,1:2500], by = list(posterior_pred8_dat1$X25011), FUN = function(i) quantile(i, probs = 0.75))
posterior_pred12 <- aggregate(x = posterior_pred8_dat2[,1:2500], by = list(posterior_pred8_dat2$X25011), FUN = function(i) quantile(i, probs = 0.75))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Sniff 75th Percentile") + xlab("Sniff Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Walk ############################
posterior_pred11 <- aggregate(x = posterior_pred9_dat1[,1:2500], by = list(posterior_pred9_dat1$X25011), FUN = function(i) quantile(i, probs = 0.75))
posterior_pred12 <- aggregate(x = posterior_pred9_dat2[,1:2500], by = list(posterior_pred9_dat2$X25011), FUN = function(i) quantile(i, probs = 0.75))

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12)

obs_sum1 <- rbind(obs_pred1, obs_pred2)

tdp_vec <- c(rep("WT", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk 75th Percentile") + xlab("Walk Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("TDP - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk 75th Percentile") +  xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "TDP - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)


mylegend<-g_legend(p9)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), 
                         p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), 
                         p5 + theme(legend.position="none"), p6 + theme(legend.position="none"),
                         p7 + theme(legend.position="none"), p8 + theme(legend.position="none"), 
                         p9 + theme(legend.position="none"),
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1))



block <- c()






