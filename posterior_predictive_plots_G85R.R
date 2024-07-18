library(ggplot2)
library(gridExtra)
library(dplyr)
library(bayesplot)
library(varhandle)
library(DirichletReg)
library(readr)
library(grid)

#### Subset the observed data


mice_new <- read_csv("AllFallonDataOct2021.csv", col_types = cols(`Day of Recording` = col_number(), 
                                                                  `Hour (0-23)` = col_number(), drink = col_number(), 
                                                                  eat = col_number(), groom = col_number(), 
                                                                  hang = col_number(), sniff = col_number(), 
                                                                  rear = col_number(), rest = col_number(), 
                                                                  walk = col_number(), eathand = col_number()))



# get dataset with the messed up genotype name and fix it
check_micenew <- mice_new %>% filter(Genotype != "G85R WT", Genotype != "G85R hom", Genotype != "G85R het")

check_micenew$Genotype <- rep("G85R WT", length(check_micenew$Genotype))

# unique(check_micenew$Genotype)

# get dataset without the messed up genotype name
check_micenew2 <- mice_new %>% filter(Genotype %in% c("G85R WT","G85R hom","G85R het"))

# unique(check_micenew2$Genotype)

# merge the fixed dataset with the subset of already good values
mice_G85R <- rbind(check_micenew2, check_micenew)


# fix names of the columns
names(mice_G85R) <- c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "Groom", "Hang", "Sniff", "Rear", "Rest", "Walk", "EBH")

# order the columns as in the original analysis
mice_G85R <- mice_G85R[ , c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")]


# all fixed so now let's make a factor so that WT is baseline
mice_G85R$GenotypeHOM <- ifelse(mice_G85R$Genotype == "G85R hom", 1, 0)
mice_G85R$GenotypeHET <- ifelse(mice_G85R$Genotype == "G85R het", 1, 0)


# only keeping mice aged 30 days to compare against original analysis
mice_G85R_p30 <- mice_G85R %>% filter(Age == "p30")

mice_G85R_p96 <- mice_G85R %>% filter(Age == "p96")

mice_G85R_p180 <- mice_G85R %>% filter(Age == "p180")

mice_G85R_p270 <- mice_G85R %>% filter(Age == "p270")

# define the age to run model for
mice_positive <- mice_G85R_p30
mice_positive <- mice_G85R_p96
mice_positive <- mice_G85R_p180
mice_positive <- mice_G85R_p270

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

outcomes_int <- smart_round(mice_positive[,8:16])

mice_positive$Y <- outcomes_int




# create a TDP vector to identify each mouse's observations
mice_tdp <- mice_positive$Genotype

table(mice_tdp)

# create a vector with how many observations we have for each specific TDP
tdp_reps <- c()

for(i in 1:3){
  tdp_reps <- c(tdp_reps, table(mice_tdp)[[i]])
}

sum(tdp_reps)

# create binary variables for each TDP identifying when they were observed
bnry_mice_tdp <- to.dummy(mice_tdp, "p30")

# replicate and rbind the matrix above by 9 (one per behaviour)
bnry9times_tdp <- do.call(rbind, replicate(9, bnry_mice_tdp, simplify = F))


load("ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p30_taus_Day_sml.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p30_taus_Day_sml[,c(275826,303680)]
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p30_taus_Day_sml[,c(275826:303680)]


load("ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p30_taus_Day_sml_shorter.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p30_taus_Day_sml_shorter[,c(275826,303680)]
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p30_taus_Day_sml_shorter[,c(275826:303680)]


load("ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p96_taus_Day_sml.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p96_taus_Day_sml[,c(286504,315438)]
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p96_taus_Day_sml[,c(286504:315438)]


load("ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p180_taus_Day_sml.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p180_taus_Day_sml[,c(298076,328180)]
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p180_taus_Day_sml[,c(298076:328180)]

load("ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p270_taus_Day_sml.Rdata")
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p270_taus_Day_sml[,c(310091,341410)]
yreps_log2 <- ZIGDM_Sep13_RSIGQ_Simple_posterior_G85R_p270_taus_Day_sml[,c(310091:341410)]


# matrix of each mouse's observations

hour <- rep((mice_positive$Hour),9)

yreps_mat <- rbind(yreps_log2, hour)


yreps_mat <- data.matrix(yreps_mat)







# Matrices with mice having only TDP or no TDP
# TDP No
G85R_het <- t(yreps_mat[ , as.logical(bnry9times_tdp[,1]), drop = F])
# TDP Yes
G85R_hom <- t(yreps_mat[ , as.logical(bnry9times_tdp[,2]), drop = F])
# NEW TDP
G85R_WT <- t(yreps_mat[ , as.logical(bnry9times_tdp[,3]), drop = F])


#################### 
##### p30 #####
#################### 

## TDP NO
posterior_pred1_dat1 <- data.frame(G85R_het[1:1086,])
posterior_pred2_dat1 <- data.frame(G85R_het[(1 + 1086*1):(1086 + 1086*1),])
posterior_pred3_dat1 <- data.frame(G85R_het[(1 + 1086*2):(1086 + 1086*2),])
posterior_pred4_dat1 <- data.frame(G85R_het[(1 + 1086*3):(1086 + 1086*3),])
posterior_pred5_dat1 <- data.frame(G85R_het[(1 + 1086*4):(1086 + 1086*4),])
posterior_pred6_dat1 <- data.frame(G85R_het[(1 + 1086*5):(1086 + 1086*5),])
posterior_pred7_dat1 <- data.frame(G85R_het[(1 + 1086*6):(1086 + 1086*6),])
posterior_pred8_dat1 <- data.frame(G85R_het[(1 + 1086*7):(1086 + 1086*7),])
posterior_pred9_dat1 <- data.frame(G85R_het[(1 + 1086*8):(1086 + 1086*8),])

mice_positive_sub_HET <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(G85R_hom[1:1053,])
posterior_pred2_dat2 <- data.frame(G85R_hom[(1 + 1053*1):(1053 + 1053*1),])
posterior_pred3_dat2 <- data.frame(G85R_hom[(1 + 1053*2):(1053 + 1053*2),])
posterior_pred4_dat2 <- data.frame(G85R_hom[(1 + 1053*3):(1053 + 1053*3),])
posterior_pred5_dat2 <- data.frame(G85R_hom[(1 + 1053*4):(1053 + 1053*4),])
posterior_pred6_dat2 <- data.frame(G85R_hom[(1 + 1053*5):(1053 + 1053*5),])
posterior_pred7_dat2 <- data.frame(G85R_hom[(1 + 1053*6):(1053 + 1053*6),])
posterior_pred8_dat2 <- data.frame(G85R_hom[(1 + 1053*7):(1053 + 1053*7),])
posterior_pred9_dat2 <- data.frame(G85R_hom[(1 + 1053*8):(1053 + 1053*8),])

mice_positive_sub_HOM <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]


## TDP YES
posterior_pred1_dat3 <- data.frame(G85R_WT[1:956,])
posterior_pred2_dat3 <- data.frame(G85R_WT[(1 + 956*1):(956 + 956*1),])
posterior_pred3_dat3 <- data.frame(G85R_WT[(1 + 956*2):(956 + 956*2),])
posterior_pred4_dat3 <- data.frame(G85R_WT[(1 + 956*3):(956 + 956*3),])
posterior_pred5_dat3 <- data.frame(G85R_WT[(1 + 956*4):(956 + 956*4),])
posterior_pred6_dat3 <- data.frame(G85R_WT[(1 + 956*5):(956 + 956*5),])
posterior_pred7_dat3 <- data.frame(G85R_WT[(1 + 956*6):(956 + 956*6),])
posterior_pred8_dat3 <- data.frame(G85R_WT[(1 + 956*7):(956 + 956*7),])
posterior_pred9_dat3 <- data.frame(G85R_WT[(1 + 956*8):(956 + 956*8),])

mice_positive_sub_WT <- mice_positive[ as.logical(bnry_mice_tdp[,3]) , , drop = F]





#################### 
##### p96 #####
#################### 

## TDP NO
posterior_pred1_dat1 <- data.frame(G85R_het[1:1175,])
posterior_pred2_dat1 <- data.frame(G85R_het[(1 + 1175*1):(1175 + 1175*1),])
posterior_pred3_dat1 <- data.frame(G85R_het[(1 + 1175*2):(1175 + 1175*2),])
posterior_pred4_dat1 <- data.frame(G85R_het[(1 + 1175*3):(1175 + 1175*3),])
posterior_pred5_dat1 <- data.frame(G85R_het[(1 + 1175*4):(1175 + 1175*4),])
posterior_pred6_dat1 <- data.frame(G85R_het[(1 + 1175*5):(1175 + 1175*5),])
posterior_pred7_dat1 <- data.frame(G85R_het[(1 + 1175*6):(1175 + 1175*6),])
posterior_pred8_dat1 <- data.frame(G85R_het[(1 + 1175*7):(1175 + 1175*7),])
posterior_pred9_dat1 <- data.frame(G85R_het[(1 + 1175*8):(1175 + 1175*8),])

mice_positive_sub_HET <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(G85R_hom[1:1080,])
posterior_pred2_dat2 <- data.frame(G85R_hom[(1 + 1080*1):(1080 + 1080*1),])
posterior_pred3_dat2 <- data.frame(G85R_hom[(1 + 1080*2):(1080 + 1080*2),])
posterior_pred4_dat2 <- data.frame(G85R_hom[(1 + 1080*3):(1080 + 1080*3),])
posterior_pred5_dat2 <- data.frame(G85R_hom[(1 + 1080*4):(1080 + 1080*4),])
posterior_pred6_dat2 <- data.frame(G85R_hom[(1 + 1080*5):(1080 + 1080*5),])
posterior_pred7_dat2 <- data.frame(G85R_hom[(1 + 1080*6):(1080 + 1080*6),])
posterior_pred8_dat2 <- data.frame(G85R_hom[(1 + 1080*7):(1080 + 1080*7),])
posterior_pred9_dat2 <- data.frame(G85R_hom[(1 + 1080*8):(1080 + 1080*8),])

mice_positive_sub_HOM <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]


## TDP YES
posterior_pred1_dat3 <- data.frame(G85R_WT[1:960,])
posterior_pred2_dat3 <- data.frame(G85R_WT[(1 + 960*1):(960 + 960*1),])
posterior_pred3_dat3 <- data.frame(G85R_WT[(1 + 960*2):(960 + 960*2),])
posterior_pred4_dat3 <- data.frame(G85R_WT[(1 + 960*3):(960 + 960*3),])
posterior_pred5_dat3 <- data.frame(G85R_WT[(1 + 960*4):(960 + 960*4),])
posterior_pred6_dat3 <- data.frame(G85R_WT[(1 + 960*5):(960 + 960*5),])
posterior_pred7_dat3 <- data.frame(G85R_WT[(1 + 960*6):(960 + 960*6),])
posterior_pred8_dat3 <- data.frame(G85R_WT[(1 + 960*7):(960 + 960*7),])
posterior_pred9_dat3 <- data.frame(G85R_WT[(1 + 960*8):(960 + 960*8),])

mice_positive_sub_WT <- mice_positive[ as.logical(bnry_mice_tdp[,3]) , , drop = F]



#################### 
##### p180 #####
#################### 

## TDP NO
posterior_pred1_dat1 <- data.frame(G85R_het[1:1073,])
posterior_pred2_dat1 <- data.frame(G85R_het[(1 + 1073*1):(1073 + 1073*1),])
posterior_pred3_dat1 <- data.frame(G85R_het[(1 + 1073*2):(1073 + 1073*2),])
posterior_pred4_dat1 <- data.frame(G85R_het[(1 + 1073*3):(1073 + 1073*3),])
posterior_pred5_dat1 <- data.frame(G85R_het[(1 + 1073*4):(1073 + 1073*4),])
posterior_pred6_dat1 <- data.frame(G85R_het[(1 + 1073*5):(1073 + 1073*5),])
posterior_pred7_dat1 <- data.frame(G85R_het[(1 + 1073*6):(1073 + 1073*6),])
posterior_pred8_dat1 <- data.frame(G85R_het[(1 + 1073*7):(1073 + 1073*7),])
posterior_pred9_dat1 <- data.frame(G85R_het[(1 + 1073*8):(1073 + 1073*8),])

mice_positive_sub_HET <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(G85R_hom[1:1199,])
posterior_pred2_dat2 <- data.frame(G85R_hom[(1 + 1199*1):(1199 + 1199*1),])
posterior_pred3_dat2 <- data.frame(G85R_hom[(1 + 1199*2):(1199 + 1199*2),])
posterior_pred4_dat2 <- data.frame(G85R_hom[(1 + 1199*3):(1199 + 1199*3),])
posterior_pred5_dat2 <- data.frame(G85R_hom[(1 + 1199*4):(1199 + 1199*4),])
posterior_pred6_dat2 <- data.frame(G85R_hom[(1 + 1199*5):(1199 + 1199*5),])
posterior_pred7_dat2 <- data.frame(G85R_hom[(1 + 1199*6):(1199 + 1199*6),])
posterior_pred8_dat2 <- data.frame(G85R_hom[(1 + 1199*7):(1199 + 1199*7),])
posterior_pred9_dat2 <- data.frame(G85R_hom[(1 + 1199*8):(1199 + 1199*8),])

mice_positive_sub_HOM <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]


## TDP YES
posterior_pred1_dat3 <- data.frame(G85R_WT[1:1073,])
posterior_pred2_dat3 <- data.frame(G85R_WT[(1 + 1073*1):(1073 + 1073*1),])
posterior_pred3_dat3 <- data.frame(G85R_WT[(1 + 1073*2):(1073 + 1073*2),])
posterior_pred4_dat3 <- data.frame(G85R_WT[(1 + 1073*3):(1073 + 1073*3),])
posterior_pred5_dat3 <- data.frame(G85R_WT[(1 + 1073*4):(1073 + 1073*4),])
posterior_pred6_dat3 <- data.frame(G85R_WT[(1 + 1073*5):(1073 + 1073*5),])
posterior_pred7_dat3 <- data.frame(G85R_WT[(1 + 1073*6):(1073 + 1073*6),])
posterior_pred8_dat3 <- data.frame(G85R_WT[(1 + 1073*7):(1073 + 1073*7),])
posterior_pred9_dat3 <- data.frame(G85R_WT[(1 + 1073*8):(1073 + 1073*8),])

mice_positive_sub_WT <- mice_positive[ as.logical(bnry_mice_tdp[,3]) , , drop = F]




#################### 
##### p270 #####
#################### 

## TDP NO
posterior_pred1_dat1 <- data.frame(G85R_het[1:1200,])
posterior_pred2_dat1 <- data.frame(G85R_het[(1 + 1200*1):(1200 + 1200*1),])
posterior_pred3_dat1 <- data.frame(G85R_het[(1 + 1200*2):(1200 + 1200*2),])
posterior_pred4_dat1 <- data.frame(G85R_het[(1 + 1200*3):(1200 + 1200*3),])
posterior_pred5_dat1 <- data.frame(G85R_het[(1 + 1200*4):(1200 + 1200*4),])
posterior_pred6_dat1 <- data.frame(G85R_het[(1 + 1200*5):(1200 + 1200*5),])
posterior_pred7_dat1 <- data.frame(G85R_het[(1 + 1200*6):(1200 + 1200*6),])
posterior_pred8_dat1 <- data.frame(G85R_het[(1 + 1200*7):(1200 + 1200*7),])
posterior_pred9_dat1 <- data.frame(G85R_het[(1 + 1200*8):(1200 + 1200*8),])

mice_positive_sub_HET <- mice_positive[ as.logical(bnry_mice_tdp[,1]) , , drop = F]



## TDP YES
posterior_pred1_dat2 <- data.frame(G85R_hom[1:1320,])
posterior_pred2_dat2 <- data.frame(G85R_hom[(1 + 1320*1):(1320 + 1320*1),])
posterior_pred3_dat2 <- data.frame(G85R_hom[(1 + 1320*2):(1320 + 1320*2),])
posterior_pred4_dat2 <- data.frame(G85R_hom[(1 + 1320*3):(1320 + 1320*3),])
posterior_pred5_dat2 <- data.frame(G85R_hom[(1 + 1320*4):(1320 + 1320*4),])
posterior_pred6_dat2 <- data.frame(G85R_hom[(1 + 1320*5):(1320 + 1320*5),])
posterior_pred7_dat2 <- data.frame(G85R_hom[(1 + 1320*6):(1320 + 1320*6),])
posterior_pred8_dat2 <- data.frame(G85R_hom[(1 + 1320*7):(1320 + 1320*7),])
posterior_pred9_dat2 <- data.frame(G85R_hom[(1 + 1320*8):(1320 + 1320*8),])

mice_positive_sub_HOM <- mice_positive[ as.logical(bnry_mice_tdp[,2]) , , drop = F]


## TDP YES
posterior_pred1_dat3 <- data.frame(G85R_WT[1:960,])
posterior_pred2_dat3 <- data.frame(G85R_WT[(1 + 960*1):(960 + 960*1),])
posterior_pred3_dat3 <- data.frame(G85R_WT[(1 + 960*2):(960 + 960*2),])
posterior_pred4_dat3 <- data.frame(G85R_WT[(1 + 960*3):(960 + 960*3),])
posterior_pred5_dat3 <- data.frame(G85R_WT[(1 + 960*4):(960 + 960*4),])
posterior_pred6_dat3 <- data.frame(G85R_WT[(1 + 960*5):(960 + 960*5),])
posterior_pred7_dat3 <- data.frame(G85R_WT[(1 + 960*6):(960 + 960*6),])
posterior_pred8_dat3 <- data.frame(G85R_WT[(1 + 960*7):(960 + 960*7),])
posterior_pred9_dat3 <- data.frame(G85R_WT[(1 + 960*8):(960 + 960*8),])

mice_positive_sub_WT <- mice_positive[ as.logical(bnry_mice_tdp[,3]) , , drop = F]




################################################################################################################ 
################################################################################################################ 
################################################################################################################ 
################################################################################################################ 

# MEDIAN WITH ERROR BAR PLOT NOT BANDS

################################################################################################################ 
################################################################################################################ 
################################################################################################################ 
################################################################################################################ 

############################ Drink ############################
obs_pred1 <- aggregate(x = mice_positive_sub_HET$Y, by = list(mice_positive_sub_HET$Hour), FUN = median)
obs_pred2 <- aggregate(x = mice_positive_sub_HOM$Y, by = list(mice_positive_sub_HOM$Hour), FUN = median)
obs_pred3 <- aggregate(x = mice_positive_sub_WT$Y, by = list(mice_positive_sub_WT$Hour), FUN = median)

names(obs_pred1) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")
names(obs_pred2) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")
names(obs_pred3) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")


# Before all of these used to be this below:
# posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:3500], by = list(posterior_pred1_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:3500], by = list(posterior_pred1_dat2$X3501), FUN = median)


posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:2500], by = list(posterior_pred1_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:2500], by = list(posterior_pred1_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred1_dat3[,1:2500], by = list(posterior_pred1_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p1 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("Drink Median") + xlab("Drink Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "HET - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=Drink, colour = "WT - Obs")) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")


############################ Eat ############################
posterior_pred11 <- aggregate(x = posterior_pred2_dat1[,1:2500], by = list(posterior_pred2_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred2_dat2[,1:2500], by = list(posterior_pred2_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred2_dat3[,1:2500], by = list(posterior_pred2_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("Eat Median") + xlab("Eat Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "HET - Obs"), size = 1) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HOM - Obs"), size = 1) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=Eat, colour = "WT - Obs"), size = 1) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=Eat, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")



############################ EBH ############################

posterior_pred11 <- aggregate(x = posterior_pred3_dat1[,1:2500], by = list(posterior_pred3_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred3_dat2[,1:2500], by = list(posterior_pred3_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred3_dat3[,1:2500], by = list(posterior_pred3_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p3 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("EBH Median") + xlab("EBH Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "HET - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=EBH, colour = "WT - Obs")) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=EBH, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")


############################ Groom ############################

posterior_pred11 <- aggregate(x = posterior_pred4_dat1[,1:2500], by = list(posterior_pred4_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred4_dat2[,1:2500], by = list(posterior_pred4_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred4_dat3[,1:2500], by = list(posterior_pred4_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p4 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("Groom Median") + xlab("Groom Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "HET - Obs"), size = 1) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HOM - Obs"), size = 1) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=Groom, colour = "WT - Obs"), size = 1) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=Groom, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")

############################ Hang ############################

posterior_pred11 <- aggregate(x = posterior_pred5_dat1[,1:2500], by = list(posterior_pred5_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred5_dat2[,1:2500], by = list(posterior_pred5_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred5_dat3[,1:2500], by = list(posterior_pred5_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p5 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("Hang Median") + xlab("Hang Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "HET - Obs"), size = 1) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HOM - Obs"), size = 1) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=Hang, colour = "WT - Obs"), size = 1) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=Hang, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")

############################ Rear ############################
posterior_pred11 <- aggregate(x = posterior_pred6_dat1[,1:2500], by = list(posterior_pred6_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred6_dat2[,1:2500], by = list(posterior_pred6_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred6_dat3[,1:2500], by = list(posterior_pred6_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p6 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("Rear Median") + xlab("Rear Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "HET - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=Rear, colour = "WT - Obs")) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=Rear, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")


############################ Rest ############################
posterior_pred11 <- aggregate(x = posterior_pred7_dat1[,1:2500], by = list(posterior_pred7_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred7_dat2[,1:2500], by = list(posterior_pred7_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred7_dat3[,1:2500], by = list(posterior_pred7_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p7 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("Rest Median") + xlab("Rest Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "HET - Obs"), size = 1) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HOM - Obs"), size = 1) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=Rest, colour = "WT - Obs"), size = 1) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=Rest, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")



############################ Sniff ############################
posterior_pred11 <- aggregate(x = posterior_pred8_dat1[,1:2500], by = list(posterior_pred8_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred8_dat2[,1:2500], by = list(posterior_pred8_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred8_dat3[,1:2500], by = list(posterior_pred8_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p8 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("Sniff Median") + xlab("Sniff Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "HET - Obs"), size = 1) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HOM - Obs"), size = 1) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), size = 1) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=Sniff, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")

############################ Walk ############################
posterior_pred11 <- aggregate(x = posterior_pred9_dat1[,1:2500], by = list(posterior_pred9_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred9_dat2[,1:2500], by = list(posterior_pred9_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred9_dat3[,1:2500], by = list(posterior_pred9_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")




p9 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("Walk Median") + xlab("Walk Hour") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=Walk, colour = "WT - Obs")) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")


# labels = c("TDP - Obs", "CTRL", "TDP", "CTRL - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
  geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
  ylab("Walk Median") + xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_point(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
  geom_point(data = obs_pred3, aes(x=Group.1, y=Walk, colour = "WT - Obs")) +
  geom_line(data = obs_pred3, aes(x=Group.1, y=Walk, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
  scale_color_brewer(palette = "Paired")


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
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1), top=textGrob("G85R HOM vs. G85R HET vs. WT Aged 1 Month Posterior Predictive Median Behavior"))

block <- c()



################################################################################################################ 
################################################################################################################ 
################################################################################################################ 
################################################################################################################ 

# MEDIAN

################################################################################################################ 
################################################################################################################ 
################################################################################################################ 
################################################################################################################ 

############################ Drink ############################
obs_pred1 <- aggregate(x = mice_positive_sub_HET$Y, by = list(mice_positive_sub_HET$Hour), FUN = median)
obs_pred2 <- aggregate(x = mice_positive_sub_HOM$Y, by = list(mice_positive_sub_HOM$Hour), FUN = median)
obs_pred3 <- aggregate(x = mice_positive_sub_WT$Y, by = list(mice_positive_sub_WT$Hour), FUN = median)

names(obs_pred1) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")
names(obs_pred2) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")
names(obs_pred3) <- c("Group.1","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")


# Before all of these used to be this below:
# posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:3500], by = list(posterior_pred1_dat1$X3501), FUN = median)
# posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:3500], by = list(posterior_pred1_dat2$X3501), FUN = median)


posterior_pred11 <- aggregate(x = posterior_pred1_dat1[,1:2500], by = list(posterior_pred1_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred1_dat2[,1:2500], by = list(posterior_pred1_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred1_dat3[,1:2500], by = list(posterior_pred1_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p1 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Drink Median") + xlab("Drink Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=Drink, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")


# p1 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + 
#   geom_errorbar(mapping = aes(x=Hour, ymin = Lower, ymax = Upper), width=0.2, size=1) + 
#   ylab("Drink Median") + xlab("Drink Hour") +
#   geom_point(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "HET - Obs")) + 
#   geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "HET - Obs"), alpha=0.75, size = 0.5) +
#   geom_point(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HOM - Obs")) + 
#   geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HOM - Obs"), alpha=0.75, size = 0.5) +
#   geom_point(data = obs_pred3, aes(x=Group.1, y=Drink, colour = "WT - Obs")) +
#   geom_line(data = obs_pred3, aes(x=Group.1, y=Drink, colour = "WT - Obs"), alpha=0.75, size = 0.5) +
#   scale_color_brewer(palette = "Paired")


############################ Eat ############################
posterior_pred11 <- aggregate(x = posterior_pred2_dat1[,1:2500], by = list(posterior_pred2_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred2_dat2[,1:2500], by = list(posterior_pred2_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred2_dat3[,1:2500], by = list(posterior_pred2_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Eat Median") + xlab("Eat Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=Eat, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")


############################ EBH ############################

posterior_pred11 <- aggregate(x = posterior_pred3_dat1[,1:2500], by = list(posterior_pred3_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred3_dat2[,1:2500], by = list(posterior_pred3_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred3_dat3[,1:2500], by = list(posterior_pred3_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("EBH Median") + xlab("EBH Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=EBH, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")



############################ Groom ############################

posterior_pred11 <- aggregate(x = posterior_pred4_dat1[,1:2500], by = list(posterior_pred4_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred4_dat2[,1:2500], by = list(posterior_pred4_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred4_dat3[,1:2500], by = list(posterior_pred4_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p4 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Groom Median") + xlab("Groom Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=Groom, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")



############################ Hang ############################

posterior_pred11 <- aggregate(x = posterior_pred5_dat1[,1:2500], by = list(posterior_pred5_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred5_dat2[,1:2500], by = list(posterior_pred5_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred5_dat3[,1:2500], by = list(posterior_pred5_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p5 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Hang Median") + xlab("Hang Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=Hang, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")



############################ Rear ############################
posterior_pred11 <- aggregate(x = posterior_pred6_dat1[,1:2500], by = list(posterior_pred6_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred6_dat2[,1:2500], by = list(posterior_pred6_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred6_dat3[,1:2500], by = list(posterior_pred6_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p6 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rear Median") + xlab("Rear Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=Rear, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")



############################ Rest ############################
posterior_pred11 <- aggregate(x = posterior_pred7_dat1[,1:2500], by = list(posterior_pred7_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred7_dat2[,1:2500], by = list(posterior_pred7_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred7_dat3[,1:2500], by = list(posterior_pred7_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p7 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rest Median") + xlab("Rest Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=Rest, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")




############################ Sniff ############################
posterior_pred11 <- aggregate(x = posterior_pred8_dat1[,1:2500], by = list(posterior_pred8_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred8_dat2[,1:2500], by = list(posterior_pred8_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred8_dat3[,1:2500], by = list(posterior_pred8_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")


p8 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Sniff Median") + xlab("Sniff Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=Sniff, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")



############################ Walk ############################
posterior_pred11 <- aggregate(x = posterior_pred9_dat1[,1:2500], by = list(posterior_pred9_dat1$X25011), FUN = median)
posterior_pred12 <- aggregate(x = posterior_pred9_dat2[,1:2500], by = list(posterior_pred9_dat2$X25011), FUN = median)
posterior_pred13 <- aggregate(x = posterior_pred9_dat3[,1:2500], by = list(posterior_pred9_dat3$X25011), FUN = median)

pred_summary11 <- data.frame(t(posterior_pred11))
pred_summary11 <- pred_summary11[-1,]
names(pred_summary11) <- c(0:23)

pred_summary12 <- data.frame(t(posterior_pred12))
pred_summary12 <- pred_summary12[-1,]
names(pred_summary12) <- c(0:23)

pred_summary13 <- data.frame(t(posterior_pred13))
pred_summary13 <- pred_summary13[-1,]
names(pred_summary13) <- c(0:23)

pp_summary11 <- cbind(apply(pred_summary11,2,median), apply(pred_summary11,2,quantile,c(0.025)), 
                      apply(pred_summary11,2,quantile,c(0.975)))

pp_summary12 <- cbind(apply(pred_summary12,2,median), apply(pred_summary12,2,quantile,c(0.025)), 
                      apply(pred_summary12,2,quantile,c(0.975)))

pp_summary13 <- cbind(apply(pred_summary13,2,median), apply(pred_summary13,2,quantile,c(0.025)), 
                      apply(pred_summary13,2,quantile,c(0.975)))

pp_summary1 <- rbind(pp_summary11, pp_summary12, pp_summary13)

obs_sum1 <- rbind(obs_pred1, obs_pred2, obs_pred3)

tdp_vec <- c(rep("HET", 24), rep("HOM", 24), rep("WT", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Median", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")




p9 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk Median") + xlab("Walk Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")



# labels = c("TDP - Obs", "CTRL", "TDP", "CTRL - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Median, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk Median") +  xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  geom_line(data = obs_pred3, aes(x=Group.1, y=Walk, colour = "WT - Obs")) + 
  scale_color_brewer(palette = "Paired")




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
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1), top=textGrob("p270 G85R Behavior Analysis"))

block <- c()





p1
p2
p3
p4
p5
p6
p7
p8
p9

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
obs_pred1 <- aggregate(x = mice_positive_sub_HET$Y, by = list(mice_positive_sub_HET$Hour), FUN = max)
obs_pred2 <- aggregate(x = mice_positive_sub_HOM$Y, by = list(mice_positive_sub_HOM$Hour), FUN = max)

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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p1 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Drink Max") + xlab("Drink Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Eat Max") + xlab("Eat Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))


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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("EBH Max") + xlab("EBH Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Groom Max") + xlab("Groom Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Hang Max") + xlab("Hang Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rear Max") + xlab("Rear Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rest Max") + xlab("Rest Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))




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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Sniff Max") + xlab("Sniff Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Max", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk Max") + xlab("Walk Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



# labels = c("TDP - Obs", "CTRL", "TDP", "CTRL - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Max, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk Max") + xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))


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
# obs_pred1 <- aggregate(x = mice_positive_sub_HET[,8:16], by = list(mice_positive_sub_HET$Hour), FUN = min)
# obs_pred2 <- aggregate(x = mice_positive_sub_HOM[,8:16], by = list(mice_positive_sub_HOM$Hour), FUN = min)

obs_pred1 <- aggregate(x = mice_positive_sub_HET$Y, by = list(mice_positive_sub_HET$Hour), FUN = min)
obs_pred2 <- aggregate(x = mice_positive_sub_HOM$Y, by = list(mice_positive_sub_HOM$Hour), FUN = min)

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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p1 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Drink Min") + xlab("Drink Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Eat Min") + xlab("Eat Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))


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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("EBH Min") + xlab("EBH Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Groom Min") + xlab("Groom Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Hang Min") + xlab("Hang Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rear Min") + xlab("Rear Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rest Min") + xlab("Rest Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))




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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Sniff Min") + xlab("Sniff Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("Min", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk Min") + xlab("Walk Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



# labels = c("TDP - Obs", "CTRL", "TDP", "CTRL - Obs")


p9 <- ggplot(pp_summary, aes(Hour, Min, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk Min") +  xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))


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
obs_pred1 <- aggregate(x = mice_positive_sub_HET$Y, by = list(mice_positive_sub_HET$Hour), 
                       FUN = function(i) quantile(i, probs = 0.25))
obs_pred2 <- aggregate(x = mice_positive_sub_HOM$Y, by = list(mice_positive_sub_HOM$Hour), 
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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p1 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Drink 25th Percentile") + xlab("Drink Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Eat 25th Percentile") + xlab("Eat Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))


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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("EBH 25th Percentile") + xlab("EBH Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Groom 25th Percentile") + xlab("Groom Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Hang 25th Percentile") + xlab("Hang Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rear 25th Percentile") + xlab("Rear Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rest 25th Percentile") + xlab("Rest Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))




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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Sniff 25th Percentile") + xlab("Sniff Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("25th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk 25th Percentile") + xlab("Walk Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



# labels = c("TDP - Obs", "CTRL", "TDP", "CTRL - Obs")


p9 <- ggplot(pp_summary, aes(Hour, `25th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk 25th Percentile") + xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))


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
obs_pred1 <- aggregate(x = mice_positive_sub_HET$Y, by = list(mice_positive_sub_HET$Hour), 
                       FUN = function(i) quantile(i, probs = 0.75))
obs_pred2 <- aggregate(x = mice_positive_sub_HOM$Y, by = list(mice_positive_sub_HOM$Hour), 
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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p1 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Drink 75th Percentile") + xlab("Drink Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Drink, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Drink, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p2 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Eat 75th Percentile") + xlab("Eat Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Eat, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Eat, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))


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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p3 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("EBH 75th Percentile") + xlab("EBH Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=EBH, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=EBH, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p4 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Groom 75th Percentile") + xlab("Groom Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Groom, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Groom, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p5 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Hang 75th Percentile") + xlab("Hang Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Hang, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Hang, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p6 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rear 75th Percentile") + xlab("Rear Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rear, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rear, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p7 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Rest 75th Percentile") + xlab("Rest Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Rest, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Rest, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))




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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")

p8 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Sniff 75th Percentile") + xlab("Sniff Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Sniff, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Sniff, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



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

tdp_vec <- c(rep("CTRL", 24), rep("TDP", 24))

pp_summary <- data.frame(pp_summary1, obs_sum1, tdp_vec)

names(pp_summary) <- c("75th Percentile", "Lower", "Upper", "Hour", 
                       "Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk", "Genotype")



p9 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk 75th Percentile") + xlab("Walk Hour") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))



# labels = c("TDP - Obs", "CTRL", "TDP", "CTRL - Obs")


p9 <- ggplot(pp_summary, aes(Hour, `75th Percentile`, colour = Genotype)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), linetype=2, alpha=0.1) + ylab("Walk 75th Percentile") +  xlab("Walk Hour") +
  theme(legend.title=element_blank(), legend.position="bottom") +
  geom_line(data = obs_pred1, aes(x=Group.1, y=Walk, colour = "HET - Obs")) + 
  geom_line(data = obs_pred2, aes(x=Group.1, y=Walk, colour = "HOM - Obs")) + 
  scale_color_manual(values = c("#F8766D", "violetred1", "#00BFC4","blue"))


mylegend<-g_legend(p9)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), 
                         p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), 
                         p5 + theme(legend.position="none"), p6 + theme(legend.position="none"),
                         p7 + theme(legend.position="none"), p8 + theme(legend.position="none"), 
                         p9 + theme(legend.position="none"),
                         nrow = 3), mylegend, nrow = 2, heights=c(14,1))



block <- c()







