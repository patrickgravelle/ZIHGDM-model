library(ggplot2)
library(gridExtra)
library(dplyr)
library(bayesplot)
library(varhandle)
library(DirichletReg)
library(readr)

#### Subset the observed data

mice_new <- read_csv("ACBM_5XFAD_2month.csv", col_types = cols(`Day of Recording` = col_number(), 
                                                               `Hour (0-23)` = col_number(), drink = col_number(), 
                                                               eat = col_number(), groom = col_number(), 
                                                               hang = col_number(), sniff = col_number(), 
                                                               rear = col_number(), rest = col_number(), 
                                                               walk = col_number(), eathand = col_number()))



mice_new <- read_csv("ACBM_5XFAD_8month.csv", col_types = cols(`Day of Recording` = col_number(), 
                                                               `Hour (0-23)` = col_number(), drink = col_number(), 
                                                               eat = col_number(), groom = col_number(), 
                                                               hang = col_number(), sniff = col_number(), 
                                                               rear = col_number(), rest = col_number(), 
                                                               walk = col_number(), eathand = col_number()))



mice_5XFAD <- mice_new

# fix names of the columns
names(mice_5XFAD) <- c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "Groom", "Hang", "Sniff", "Rear", "Rest", "Walk", "EBH")

# order the columns as in the original analysis
mice_5XFAD <- mice_5XFAD[ , c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")]


# all fixed so now let's make a factor so that WT is baseline
mice_5XFAD$GenotypeFAD <- ifelse(mice_5XFAD$Genotype == "Hemi", 1, 0)


mice_5XFAD_p60 <- mice_5XFAD
mice_5XFAD_p60$id <- match(mice_5XFAD_p60$MouseID, unique(mice_5XFAD_p60$MouseID))
mice_positive <- mice_5XFAD_p60


mice_5XFAD_p180 <- mice_5XFAD
mice_5XFAD_p180$id <- match(mice_5XFAD_p180$MouseID, unique(mice_5XFAD_p180$MouseID))
mice_positive <- mice_5XFAD_p180








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



# mice_positive$Y
# 
# check <- rowSums(mice_positive$Y)
# 
# table(check)


# create a TDP vector to identify each mouse's observations
mice_tdp <- mice_positive$GenotypeFAD

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






# 5XFAD p60 model
load("ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p60_taus_sml.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p60_taus_sml[,c(249750,274976)]
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p60_taus_sml[,c(249750:274976)]
# OR can load the subset of data directly after saving it in CSS-HPC01
load("ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p60_taus_sml_draws.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p60_taus_sml_draws




# 5XFAD p180 model
load("ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus_sml.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus_sml[,c(198745,218814)]
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus_sml[,c(198745:218814)]
# OR can load the subset of data directly after saving it in CSS-HPC01
load("ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus_sml_draws.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus_sml_draws








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
# p60
#########################

## TDP NO
posterior_pred1_dat1 <- data.frame(TDP_no[1:1414,])
posterior_pred2_dat1 <- data.frame(TDP_no[(1 + 1414*1):(1414 + 1414*1),])
posterior_pred3_dat1 <- data.frame(TDP_no[(1 + 1414*2):(1414 + 1414*2),])
posterior_pred4_dat1 <- data.frame(TDP_no[(1 + 1414*3):(1414 + 1414*3),])
posterior_pred5_dat1 <- data.frame(TDP_no[(1 + 1414*4):(1414 + 1414*4),])
posterior_pred6_dat1 <- data.frame(TDP_no[(1 + 1414*5):(1414 + 1414*5),])
posterior_pred7_dat1 <- data.frame(TDP_no[(1 + 1414*6):(1414 + 1414*6),])
posterior_pred8_dat1 <- data.frame(TDP_no[(1 + 1414*7):(1414 + 1414*7),])
posterior_pred9_dat1 <- data.frame(TDP_no[(1 + 1414*8):(1414 + 1414*8),])

mice_positive_sub_no <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(TDP_yes[1:1389,])
posterior_pred2_dat2 <- data.frame(TDP_yes[(1 + 1389*1):(1389 + 1389*1),])
posterior_pred3_dat2 <- data.frame(TDP_yes[(1 + 1389*2):(1389 + 1389*2),])
posterior_pred4_dat2 <- data.frame(TDP_yes[(1 + 1389*3):(1389 + 1389*3),])
posterior_pred5_dat2 <- data.frame(TDP_yes[(1 + 1389*4):(1389 + 1389*4),])
posterior_pred6_dat2 <- data.frame(TDP_yes[(1 + 1389*5):(1389 + 1389*5),])
posterior_pred7_dat2 <- data.frame(TDP_yes[(1 + 1389*6):(1389 + 1389*6),])
posterior_pred8_dat2 <- data.frame(TDP_yes[(1 + 1389*7):(1389 + 1389*7),])
posterior_pred9_dat2 <- data.frame(TDP_yes[(1 + 1389*8):(1389 + 1389*8),])

mice_positive_sub_yes <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]



#########################
# p180
#########################

## TDP NO
posterior_pred1_dat1 <- data.frame(TDP_no[1:1120,])
posterior_pred2_dat1 <- data.frame(TDP_no[(1 + 1120*1):(1120 + 1120*1),])
posterior_pred3_dat1 <- data.frame(TDP_no[(1 + 1120*2):(1120 + 1120*2),])
posterior_pred4_dat1 <- data.frame(TDP_no[(1 + 1120*3):(1120 + 1120*3),])
posterior_pred5_dat1 <- data.frame(TDP_no[(1 + 1120*4):(1120 + 1120*4),])
posterior_pred6_dat1 <- data.frame(TDP_no[(1 + 1120*5):(1120 + 1120*5),])
posterior_pred7_dat1 <- data.frame(TDP_no[(1 + 1120*6):(1120 + 1120*6),])
posterior_pred8_dat1 <- data.frame(TDP_no[(1 + 1120*7):(1120 + 1120*7),])
posterior_pred9_dat1 <- data.frame(TDP_no[(1 + 1120*8):(1120 + 1120*8),])

mice_positive_sub_no <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(TDP_yes[1:1110,])
posterior_pred2_dat2 <- data.frame(TDP_yes[(1 + 1110*1):(1110 + 1110*1),])
posterior_pred3_dat2 <- data.frame(TDP_yes[(1 + 1110*2):(1110 + 1110*2),])
posterior_pred4_dat2 <- data.frame(TDP_yes[(1 + 1110*3):(1110 + 1110*3),])
posterior_pred5_dat2 <- data.frame(TDP_yes[(1 + 1110*4):(1110 + 1110*4),])
posterior_pred6_dat2 <- data.frame(TDP_yes[(1 + 1110*5):(1110 + 1110*5),])
posterior_pred7_dat2 <- data.frame(TDP_yes[(1 + 1110*6):(1110 + 1110*6),])
posterior_pred8_dat2 <- data.frame(TDP_yes[(1 + 1110*7):(1110 + 1110*7),])
posterior_pred9_dat2 <- data.frame(TDP_yes[(1 + 1110*8):(1110 + 1110*8),])

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


# Before all of these used to be this below:
# posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:2500], by = list(posterior_pred1_dat1$X25011), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:2500], by = list(posterior_pred1_dat2$X25011), FUN = median)


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

tdp_vec <- c(rep("WT", 24), rep("HEMI", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p1 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() +  
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Drink Median") + xlab("Drink Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs")) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) +
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) + 
  scale_colour_manual(values = my_colors)

# checkpoint




############################ Eat ############################
posterior_pred11 <- aggregate(x = posterior_pred2_dat1[,1:2500], by = list(posterior_pred2_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred2_dat2[,1:2500], by = list(posterior_pred2_dat2$X25011), FUN = median)

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

tdp_vec <- c(rep("WT", 24), rep("HEMI", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Eat Median") + xlab("Eat Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs")) +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs")) +
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) +
  scale_colour_manual(values = my_colors)


############################ EBH ############################

posterior_pred11 <- aggregate(x = posterior_pred3_dat1[,1:2500], by = list(posterior_pred3_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred3_dat2[,1:2500], by = list(posterior_pred3_dat2$X25011), FUN = median)

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

tdp_vec <- c(rep("WT", 24), rep("HEMI", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("EBH Median") + xlab("EBH Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs")) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Groom ############################

posterior_pred11 <- aggregate(x = posterior_pred4_dat1[,1:2500], by = list(posterior_pred4_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred4_dat2[,1:2500], by = list(posterior_pred4_dat2$X25011), FUN = median)

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

tdp_vec <- c(rep("WT", 24), rep("HEMI", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Groom Median") + xlab("Groom Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Hang ############################

posterior_pred11 <- aggregate(x = posterior_pred5_dat1[,1:2500], by = list(posterior_pred5_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred5_dat2[,1:2500], by = list(posterior_pred5_dat2$X25011), FUN = median)

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

tdp_vec <- c(rep("WT", 24), rep("HEMI", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Hang Median") + xlab("Hang Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rear ############################
posterior_pred11 <- aggregate(x = posterior_pred6_dat1[,1:2500], by = list(posterior_pred6_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred6_dat2[,1:2500], by = list(posterior_pred6_dat2$X25011), FUN = median)

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

tdp_vec <- c(rep("WT", 24), rep("HEMI", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Rear Median") + xlab("Rear Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Rest ############################
posterior_pred11 <- aggregate(x = posterior_pred7_dat1[,1:2500], by = list(posterior_pred7_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred7_dat2[,1:2500], by = list(posterior_pred7_dat2$X25011), FUN = median)

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

tdp_vec <- c(rep("WT", 24), rep("HEMI", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Rest Median") + xlab("Rest Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) + 
  scale_colour_manual(values = my_colors)




############################ Sniff ############################
posterior_pred11 <- aggregate(x = posterior_pred8_dat1[,1:2500], by = list(posterior_pred8_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred8_dat2[,1:2500], by = list(posterior_pred8_dat2$X25011), FUN = median)

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

tdp_vec <- c(rep("WT", 24), rep("HEMI", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Sniff Median") + xlab("Sniff Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) + 
  scale_colour_manual(values = my_colors)



############################ Walk ############################
posterior_pred11 <- aggregate(x = posterior_pred9_dat1[,1:2500], by = list(posterior_pred9_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred9_dat2[,1:2500], by = list(posterior_pred9_dat2$X25011), FUN = median)

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

tdp_vec <- c(rep("WT", 24), rep("HEMI", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Walk Median") + xlab("Walk Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("HEMI - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, linewidth=1) + ylab("Walk Median") +  xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, linewidth = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, linewidth = 0.5) + 
  scale_colour_manual(values = my_colors)



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
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1), top=("5XFAD - HEMI vs. WT Aged 2 Months Posterior Predictive Median Behavior"))


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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("HEMI - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk Max") + xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs")) +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) +
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("HEMI - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk Min") +  xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs")) +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) +
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("HEMI - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk 25th Percentile") + xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
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
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)



# labels = c("HEMI - Obs", "WT", "TDP", "WT - Obs")


p9 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + ylab("Walk 75th Percentile") +  xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HEMI - Obs"), alpha=0.75, size = 0.5) + 
  scale_colour_manual(values = my_colors)


mylegend<-g_legend(p9)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), 
                         p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), 
                         p5 + theme(legend.position="none"), p6 + theme(legend.position="none"),
                         p7 + theme(legend.position="none"), p8 + theme(legend.position="none"), 
                         p9 + theme(legend.position="none"),
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1))



block <- c()








