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


mice_new <- read_csv("AllFallonDataOct2021.csv", col_types = cols(`Day of Recording` = col_number(), 
                                                                  `Hour (0-23)` = col_number(), drink = col_number(), 
                                                                  eat = col_number(), groom = col_number(), 
                                                                  hang = col_number(), sniff = col_number(), 
                                                                  rear = col_number(), rest = col_number(), 
                                                                  walk = col_number(), eathand = col_number()))


# get dataset with the messed up genotype name and fix it
check_micenew <- mice_new %>% filter(Genotype != "G85R WT", Genotype != "G85R hom", Genotype != "G85R het")

check_micenew$Genotype <- rep("G85R WT", length(check_micenew$Genotype))


# get dataset without the messed up genotype name
check_micenew2 <- mice_new %>% filter(Genotype %in% c("G85R WT","G85R hom","G85R het"))


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
mice_G85R_p180 <- mice_G85R %>% filter(Age == "p180")

# no so we have to make a batch vector 
check_ymd <- ymd(mice_G85R_p180$Date)

batch_vec <- ifelse(check_ymd[,2] == 2, 1, 
                    ifelse(check_ymd[,2] == 4, 2, 
                           ifelse(check_ymd[,2] == 5, 3, 4)))



mice_G85R_p180 <- data.frame(mice_G85R_p180, batch_vec)


# Have to make the mouse IDs 1:unique(mouseID)

mice_G85R_p180$id <- match(mice_G85R_p180$MouseID, unique(mice_G85R_p180$MouseID))



mice_positive <- mice_G85R_p180



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


X <- as.matrix(model.matrix(lm(Drink ~ GenotypeHOM + GenotypeHET + Day, data = mice_positive)))

X <- matrix(nrow = nrow(X), ncol = ncol(X), data = as.numeric(X))





mice_id <- mice_G85R_p180$id

nd <- length(unique(mice_id))

days <- mice_G85R_p180$Day

ndays <- length(unique(days))

batch <- mice_G85R_p180$batch_vec

nbatch <- length(unique(batch))



hour_dat <- mice_G85R_p180$Hour

num_data <- length(hour_dat)




sspec <- s(Hour, bs = "cp")

sm <- smoothCon(sspec, data = mice_G85R_p180)[[1]]

num_basis <- nrow(t(sm$X))


data4Stan <- list(N = nrow(Y), ncolY = ncol(Y), ncolX = ncol(X),
                  X = X, Y = Y,
                  num_data = num_data, # length of hour vector
                  num_basis = num_basis, # nrow of bsplines
                  Bsplines = t(sm$X),
                  nd = nd,
                  mice_id = mice_id, # mouse id vector
                  days = days, # days id vector
                  ndays = ndays,
                  batch = batch,
                  nbatch = nbatch) # ,delta = delta_mat)




#running STAN with wishart prior for covariance matrix
fitCalc <- stan(
        file = "ZIHGDM_model_Stan_file_G85R.stan",  # Stan program
        # file = "Dirichlet_adjusted_Xdat_Splines_Interaction_REs_less_GDM.stan",  # Stan program
        data = data4Stan,    # named list of data
        chains = 1,             # number of Markov chains
        warmup = 3000,          # number of warmup iterations per chain
        iter = 13000,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0,             # no progress shown
        control = list(max_treedepth = 15, adapt_delta = 0.99))


# fitCalc

saveRDS(fitCalc, "ZIHGDM_model_Stan_file_G85R.rds")

#obtaining posterior draws  
ZIHGDM_posterior_G85R <- data.frame(as.matrix(fitCalc))

#save the posterior draws
index <- seq(1, 10000, by = 4)

ZIHGDM_posterior_G85R_sml <- ZIHGDM_posterior_G85R[index,]

save(ZIHGDM_posterior_G85R_sml, file = "ZIHGDM_posterior_G85R_sml.Rdata")



######################################################################
######################################################################
