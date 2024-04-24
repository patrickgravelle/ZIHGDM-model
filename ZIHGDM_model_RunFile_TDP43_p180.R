library(rstan)
library(reshape2)
library(dplyr)
library(splines)
library(mgcv)
library(DirichletReg)
library(readr)
library(mondate)

######################################################################
# TDP-43 Mouse Trial
######################################################################
#loading data file



mice_new <- read_csv("AllFallonTDP43Dec2021.csv", col_types = cols(`Day of Recording` = col_number(), 
                                                                  `Hour (0-23)` = col_number(), drink = col_number(), 
                                                                  eat = col_number(), groom = col_number(), 
                                                                  hang = col_number(), sniff = col_number(), 
                                                                  rear = col_number(), rest = col_number(), 
                                                                  walk = col_number(), eathand = col_number()))


# # how many mice do we have
# length(unique(mice_new$`Mouse ID`))
# # 
# # how many genotypes do we have
# unique(mice_new$Genotype)
# 
# # how many mice do we have per genotype
# table(mice_new$Genotype)


# table(mice_new$`Day of Recording`)
# 
# check_tdp <- mice_new %>% filter(`Day of Recording` %in% c(-9, -4))
# 
# check_tdp_mouse <- mice_new %>% filter(`Mouse ID` == "WT_7")



mice_TDP43 <- mice_new %>% filter(`Day of Recording` > 0)

# should remove the observation with negative day of recording for WT_9 and WT_7
# the negative day corresponds to the number of days before the first recording took place in that batch
# so it's like the camera turned on accidently or for a test so that data doesn't need to be there



# fix names of the columns
names(mice_TDP43) <- c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "Groom", "Hang", "Sniff", "Rear", "Rest", "Walk", "EBH")

# order the columns as in the original analysis
mice_TDP43 <- mice_TDP43[ , c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")]


# all fixed so now let's make a factor so that WT is baseline
mice_TDP43$GenotypeTDP <- ifelse(mice_TDP43$Genotype == "TDP43_Q331K", 1, 0)


mice_TDP43$AgeVec <- ifelse(mice_TDP43$Age == 6, "p180", 
                            ifelse(mice_TDP43$Age == 7.5, "p225", 
                                   ifelse(mice_TDP43$Age == 10, "p300", "p345")))


mice_TDP43$GenderVec <- ifelse(mice_TDP43$Gender == "Male", 1, 0)

# only keeping mice aged p180 days 
mice_TDP43_p180 <- mice_TDP43 %>% filter(AgeVec == "p180")

mice_TDP43_p225 <- mice_TDP43 %>% filter(AgeVec == "p225")

mice_TDP43_p300 <- mice_TDP43 %>% filter(AgeVec == "p300")

mice_TDP43_p345 <- mice_TDP43 %>% filter(AgeVec == "p345")


# # Every age was just one batch
# # were all the mice recorded in one batch
# table(mice_TDP43_p180$Date)
# 
# table(mice_TDP43_p225$Date)
# 
# table(mice_TDP43_p300$Date)
# 
# table(mice_TDP43_p345$Date)
# 
# # no so we have to make a batch vector 
# check_ymd <- ymd(mice_TDP43_p180$Date)
# 
# batch_vec <- ifelse(check_ymd[,2] == 9, 1, 
#                     ifelse(check_ymd[,2] == 11, 2, 3))
# 
# # # checking to make sure it worked
# # table(batch_vec)
# # table(check_ymd[,2])
# 
# mice_TDP43_p180 <- data.frame(mice_TDP43_p180, batch_vec)


# Have to make the mouse IDs 1:unique(mouseID)

mice_TDP43_p180$id <- match(mice_TDP43_p180$MouseID, unique(mice_TDP43_p180$MouseID))


# table(mice_TDP43_p180$id)
# 
# table(mice_TDP43_p180$MouseID)
# 
# unique(mice_TDP43_p180$MouseID)



mice_positive <- mice_TDP43_p180



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

Y <- outcomes_int

# check <- rowSums(outcomes_int)
# table(check)

X <- as.matrix(model.matrix(lm(Drink ~ GenderVec + GenotypeTDP, data = mice_positive)))

X <- matrix(nrow = nrow(X), ncol = ncol(X), data = as.numeric(X))

# Y <- mice_positive$Y





mice_id <- mice_TDP43_p180$id

nd <- length(unique(mice_id))

days <- mice_TDP43_p180$Day

ndays <- length(unique(days))

# batch <- mice_TDP43_p180$batch_vec
# 
# nbatch <- length(unique(batch))



hour_dat <- mice_TDP43_p180$Hour

num_data <- length(hour_dat)




sspec <- s(Hour, bs = "cp")

sm <- smoothCon(sspec, data = mice_TDP43_p180)[[1]]

num_basis <- nrow(t(sm$X))


# cbind(sm$X, hour_dat)

data4Stan <- list(N = nrow(Y), ncolY = ncol(Y), ncolX = ncol(X),
                  X = X, Y = Y,
                  num_data = num_data, # length of hour vector
                  num_basis = num_basis, # nrow of bsplines
                  Bsplines = t(sm$X),
                  nd = nd,
                  mice_id = mice_id, # mouse id vector
                  days = days, # days id vector
                  ndays = ndays)
                  # batch = batch,
                  # nbatch = nbatch) # ,delta = delta_mat)



#running STAN with wishart prior for covariance matrix
fitCalc <- stan(
        file = "ZIHGDM_model_Stan_file_TDP43.stan",  # Stan program
        # file = "Dirichlet_adjusted_Xdat_Splines_Interaction_REs_less_GDM.stan",  # Stan program
        data = data4Stan,    # named list of data
        chains = 1,             # number of Markov chains
        warmup = 2000,          # number of warmup iterations per chain
        iter = 12000,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0,             # no progress shown
        control = list(max_treedepth = 15, adapt_delta = 0.95))



saveRDS(fitCalc, "ZIGDM_Sep13_RSIGQ_Simple_TDP43_p180_noDay_taus.rds")

#obtaining posterior draws  
ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus <- data.frame(as.matrix(fitCalc))


index <- seq(1, 10000, by = 4)

ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus_sml <- ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus[index,]

save(ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus_sml, file = "ZIGDM_Sep13_RSIGQ_Simple_posterior_TDP43_p180_noDay_taus_sml.Rdata")


