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




mice_new <- read_csv("ACBM_5XFAD_8month.csv", col_types = cols(`Day of Recording` = col_number(), 
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
# 
# 
# table(mice_new$`Day of Recording`)
# 
# table(mice_new$Gender)


mice_5XFAD <- mice_new

# fix names of the columns
names(mice_5XFAD) <- c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "Groom", "Hang", "Sniff", "Rear", "Rest", "Walk", "EBH")

# order the columns as in the original analysis
mice_5XFAD <- mice_5XFAD[ , c("MouseID", "Genotype", "Gender", "Age", "Day", "Hour", "Date","Drink", "Eat", "EBH", "Groom", "Hang", "Rear", "Rest", "Sniff", "Walk")]


# all fixed so now let's make a factor so that WT is baseline
mice_5XFAD$GenotypeFAD <- ifelse(mice_5XFAD$Genotype == "Hemi", 1, 0)

mice_5XFAD_p180 <- mice_5XFAD

# mice_5XFAD$AgeVec <- ifelse(mice_5XFAD$Age == 6, "p180", 
#                             ifelse(mice_5XFAD$Age == 7.5, "p225", 
#                                    ifelse(mice_5XFAD$Age == 10, "p300", "p345")))
# 
# 
# mice_5XFAD$GenderVec <- ifelse(mice_5XFAD$Gender == "Male", 1, 0)
# 
# # only keeping mice aged p180 days 
# mice_5XFAD_p180 <- mice_5XFAD %>% filter(AgeVec == "p180")
# 
# mice_5XFAD_p225 <- mice_5XFAD %>% filter(AgeVec == "p225")
# 
# mice_5XFAD_p300 <- mice_5XFAD %>% filter(AgeVec == "p300")
# 
# mice_5XFAD_p345 <- mice_5XFAD %>% filter(AgeVec == "p345")


# # Every age was just one batch
# # were all the mice recorded in one batch
# table(mice_5XFAD_p180$Date)
# 
# table(mice_5XFAD_p225$Date)
# 
# table(mice_5XFAD_p300$Date)
# 
# table(mice_5XFAD_p345$Date)
# 
# # no so we have to make a batch vector 
# check_ymd <- ymd(mice_5XFAD_p180$Date)
# 
# batch_vec <- ifelse(check_ymd[,2] == 9, 1, 
#                     ifelse(check_ymd[,2] == 11, 2, 3))
# 
# # # checking to make sure it worked
# # table(batch_vec)
# # table(check_ymd[,2])
# 
# mice_5XFAD_p180 <- data.frame(mice_5XFAD_p180, batch_vec)


# Have to make the mouse IDs 1:unique(mouseID)

mice_5XFAD_p180$id <- match(mice_5XFAD_p180$MouseID, unique(mice_5XFAD_p180$MouseID))


# table(mice_5XFAD_p180$id)
# 
# table(mice_5XFAD_p180$MouseID)
# 
# unique(mice_5XFAD_p180$MouseID)





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

outcomes_int <- smart_round(mice_positive[,8:16])

Y <- outcomes_int

# check <- rowSums(outcomes_int)
# table(check)

# All mice are male so just have intercept and Genotype for fixed effects
X <- as.matrix(model.matrix(lm(Drink ~ GenotypeFAD, data = mice_positive)))

X <- matrix(nrow = nrow(X), ncol = ncol(X), data = as.numeric(X))





mice_id <- mice_5XFAD_p180$id

nd <- length(unique(mice_id))

days <- mice_5XFAD_p180$Day

ndays <- length(unique(days))

# batch <- mice_5XFAD_p180$batch_vec
# 
# nbatch <- length(unique(batch))



hour_dat <- mice_5XFAD_p180$Hour

num_data <- length(hour_dat)



sspec <- s(Hour, bs = "cp")

sm <- smoothCon(sspec, data = mice_5XFAD_p180)[[1]]

num_basis <- nrow(t(sm$X))


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
        file = "ZIHGDM_model_Stan_file_5XFAD.stan",  # Stan program
        data = data4Stan,    # named list of data
        chains = 1,             # number of Markov chains
        warmup = 2000,          # number of warmup iterations per chain
        iter = 12000,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0,             # no progress shown
        control = list(max_treedepth = 15, adapt_delta = 0.95))


# Give new names (these are old) but make sure to use the new names in the corresponding files

saveRDS(fitCalc, "ZIGDM_Sep13_RSIGQ_Simple_5XFAD_p180_taus.rds")

#obtaining posterior draws  
ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus <- data.frame(as.matrix(fitCalc))

index <- seq(1, 10000, by = 4)

ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus_sml <- ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus[index,]

save(ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus_sml, file = "ZIGDM_Sep13_RSIGQ_Simple_posterior_5XFAD_p180_taus_sml.Rdata")


