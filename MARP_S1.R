######################################################
##### Many Analysts Religion Project: Stage 1  #######
######################################################

### Team 008
### Script written by Theiss Bendixen and Benjamin G. Purzycki
### Contact email:  tb@cas.au.dk

### Last updated: 22nd December, 2020

library(truncnorm)
library(brms)
library(brmstools)

##################
### Simulation ###
##################

N <- 9600 # Number of participants, 400 individuals from 24 countries.

ID <- as.factor(rep(1:N))
Country <- rep(1:24, each=400)
Gender <- as.numeric(paste(sample(c("1", "2", "3"), N, replace=T)))
Age <- round(rtruncnorm(N, a=18, b=90, mean = 30, sd = 10), 0)
SES <- rtruncnorm(N, a=1, b=10, mean = 5, sd = 2)
Edu <- rtruncnorm(N, a=1, b=7, mean = 4, sd = 2)
AC <- round(rtruncnorm(N, a=0, b=1, mean = .9, sd = .25), 0)
Rel <- rtruncnorm(N, a=0, b=1, mean = .5, sd = .7)
Norm <- rtruncnorm(N, a=0, b=1, mean = .5, sd = .7)

# Setting arbitrary simulated effects (very crude!)
alpha <- rnorm(N, mean = 30, sd = 2) # Intercept
beta1 <- -1 # Religiosity
beta2 <- 1 # Cultural Norms
beta3 <- 2 # Interaction btw. Rel and Norm
beta4 <- 1 # Attention Check
beta5 <- 1 # Gender
beta6 <- 0 # Age
beta7 <- 1 # Education
beta8 <- -1 # SES

# Constructing the well-being response variable (also very crude!)
WB <- alpha + beta1*Rel + beta2*Norm + beta3*Rel*Norm + beta4*AC + beta5*Gender + beta6*Age + beta7*Edu + beta8*SES

# Binding the simulated data in a data frame
data <- data.frame(ID, Country, Gender, Age, SES, Edu, AC, Rel, Norm, WB)
data$Gender <- as.integer(data$Gender) # gender as index
data$Country <- as.integer(data$Country) # country as index
str(data)

# Group-mean centering
data$Age.grp.mean <- ave(data$Age, data$Country)
data$Age.c <- data$Age - data$Age.grp.mean
data$SES.grp.mean <- ave(data$SES, data$Country)
data$SES.c <- data$SES - data$SES.grp.mean 
data$Edu.grp.mean <- ave(data$Edu, data$Country)
data$Edu.c <- data$Edu - data$Edu.grp.mean

# Reverse-scoring attention check, so that "passed" (= 0) is the reference.
# This also reverses the posterior estimate.
data$AC.r <- abs(data$AC - 1)
data$AC.r <- as.integer(data$AC.r)

#############################
### Modelling: Scenario 1 ###
#############################

### Model 1 ###

# Setting priors
aprior <- set_prior("normal(30,2)", class = "Intercept")
bprior <- set_prior("normal(0,1)", class = "b")
lkjprior <- set_prior("lkj_corr_cholesky(2)", class = "cor")
sdprior <-set_prior("exponential(1)", class = "sd")
sigmaprior <-set_prior("exponential(1)", class = "sigma")

# Quick prior predictive check
m1.pp <- brm(WB ~ Rel*Norm + (Rel*Norm | Country) 
             + Gender + Age.c + SES.c + Edu.c + AC.r,
              data = data, 
              family = gaussian(), 
              prior = aprior + bprior + lkjprior + sdprior + sigmaprior, 
              sample_prior = "only")

pp_check(m1.pp, nsamples = 30)

# Fit the model
m1 <- brm(WB ~ Rel*Norm + (Rel*Norm | Country) 
          + Gender + Age.c + SES.c + Edu.c + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 2 ###

# Quick prior predictive check
m2.pp <- brm(WB ~ Rel*Norm + (Rel*Norm | Country) + AC.r,
             data = data, 
             family = gaussian(), 
             prior = aprior + bprior + lkjprior + sdprior + sigmaprior, 
             sample_prior = "only")

pp_check(m2.pp, nsamples = 30)

# Fit the model
m2 <- brm(WB ~ Rel*Norm + (Rel*Norm | Country) + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 3 ###

# Quick prior predictive check
m3.pp <- brm(WB ~ Rel*Norm + (1 | Country) 
             + Gender + Age.c + SES.c + Edu.c + AC.r,
             data = data, 
             family = gaussian(), 
             prior = aprior + bprior + sdprior + sigmaprior, 
             sample_prior = "only")

pp_check(m3.pp, nsamples = 30)

# Fit the model
m3 <- brm(WB ~ Rel * Norm + (1 | Country) 
          + Gender + Age.c + SES.c + Edu.c + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 4 ###

# Quick prior predictive check
m4.pp <- brm(WB ~ Rel*Norm + (1 | Country) + AC.r,
             data = data, 
             family = gaussian(), 
             prior = aprior + bprior + sdprior + sigmaprior, 
             sample_prior = "only")

pp_check(m4.pp, nsamples = 30)

# Fit the model
m4 <- brm(WB ~ Rel*Norm + (1 | Country) + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 5 ###

# Quick prior predictive check
m5.pp <- brm(WB ~ Rel + Norm + (1 | Country) 
             + Gender + Age.c + SES.c + Edu.c + AC.r,
             data = data, 
             family = gaussian(), 
             prior = aprior + bprior + sdprior + sigmaprior, 
             sample_prior = "only")

pp_check(m5.pp, nsamples = 30)

# Fit the model
m5 <- brm(WB ~ Rel + Norm + (1 | Country) 
          + Gender + Age.c + SES.c + Edu.c + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 6 ###

# Quick prior predictive check
m6.pp <- brm(WB ~ Rel + Norm + (1 | Country) + AC.r,
             data = data, 
             family = gaussian(), 
             prior = aprior + bprior + sdprior + sigmaprior, 
             sample_prior = "only")

pp_check(m6.pp, nsamples = 30)

# Fit the model
m6 <- brm(WB ~ Rel + Norm + (1 | Country) + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# Quick summaries
summary(m1)
ranef(m1)
conditional_effects(m1)
summary(m2)
ranef(m2)
conditional_effects(m2)
summary(m3)
ranef(m3)
conditional_effects(m3)
summary(m4)
ranef(m4)
conditional_effects(m4)
summary(m5)
ranef(m5)
conditional_effects(m5)
summary(m6)
ranef(m6)
conditional_effects(m6)

# Quick diagnostic plots
plot(m1)
plot(m2)
plot(m3)
plot(m4)
plot(m5)
plot(m6)

# Quick visualization option: forest plot of coefficients across countries
forest(m1)
forest(m2)
forest(m3)
forest(m4)
forest(m5)
forest(m6)

# Quick posterior predictive checks with 30 draws
pp_check(m1, nsamples = 30)
pp_check(m2, nsamples = 30)
pp_check(m3, nsamples = 30)
pp_check(m4, nsamples = 30)
pp_check(m5, nsamples = 30)
pp_check(m6, nsamples = 30)

### Model Comparison
(loo_m1 <- loo(m1))
(loo_m2 <- loo(m2))
(loo_m3 <- loo(m3))
(loo_m4 <- loo(m4))
(loo_m5 <- loo(m5))
(loo_m6 <- loo(m6))

loo_compare(loo_m1, loo_m2, loo_m3, loo_m4, loo_m5, loo_m6)
round(model_weights(m1, m2, m3, m4, m5, m6, weights = "loo"), 3)

#############################
### Modelling: Scenario 2 ###
#############################

# Collapsing Rel and Norm onto the same 0-1 scale
data$Rel2 <- (data$Rel + data$Norm)/2

# Group-mean center Religiosity
data$Rel2.grp.mean <- ave(data$Rel2, data$Country)
data$Rel2.c <- data$Rel2 - data$Rel2.grp.mean

### Model 7 ###

# Quick prior predictive check
m7.pp <- brm(WB ~ Rel2.c + (Rel2.c | Country) 
             + Gender + Age.c + SES.c + Edu.c + AC.r,
             data = data, 
             family = gaussian(), 
             prior = aprior + bprior + lkjprior + sdprior + sigmaprior, 
             sample_prior = "only")

pp_check(m7.pp, nsamples = 30)

# Fit the model
m7 <- brm(WB ~ Rel2.c + (Rel2.c | Country) 
          + Gender + Age.c + SES.c + Edu.c + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 8 ###

# Quick prior predictive check
m8.pp <- brm(WB ~ Rel2.c + (Rel2.c | Country) + AC.r,
             data = data, 
             family = gaussian(), 
             prior = aprior + bprior + lkjprior + sdprior + sigmaprior, 
             sample_prior = "only")

pp_check(m8.pp, nsamples = 30)

# Fit the model
m8 <- brm(WB ~ Rel2.c + (Rel2.c | Country) + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 9 ###

# Quick prior predictive check
m9.pp <- brm(WB ~ Rel2.c + (1 | Country) 
             + Gender + Age.c + SES.c + Edu.c + AC.r,
             data = data, 
             family = gaussian(), 
             prior = aprior + bprior + lkjprior + sdprior + sigmaprior, 
             sample_prior = "only")

pp_check(m9.pp, nsamples = 30)

# Fit the model
m9 <- brm(WB ~ Rel2.c + (1 | Country) 
          + Gender + Age.c + SES.c + Edu.c + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 10 ###

# Quick prior predictive check
m10.pp <- brm(WB ~ Rel.c + (1 | Country) + AC.r,
             data = data, 
             family = gaussian(), 
             prior = aprior + bprior + sdprior + sigmaprior, 
             sample_prior = "only")

pp_check(m10.pp, nsamples = 30)

# Fit the model
m10 <- brm(WB ~ Rel2.c + (1 | Country) + AC.r,
          data = data,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# See plotting, dignostics, posterior predictive checks, LOO,
# etc., from Scenario 1.

########################
### Plotting Model 1 ###
########################

labs <- c("Gender", "Age.c", "SES.c", "Edu.c", "AC.r",
          "Rel.", "Norm.", "Rel.*Norm.")
x <- 1:8
COEF <- c(1.01, -0.00, -1.02, 1.02, 
          -0.99, -0.77  , 1.13 , 1.64)
LL <- c(0.96 , -0.01 , -1.04, 0.99 ,
        -1.14, -1.05, 0.83, 1.14)
UL <- c(1.05 , 0.00,  -0.99, 1.05,
        -0.84, -0.48,  1.42, 2.14)
LS <- COEF-LL
US <- UL-COEF

tab <- data.frame(cbind(labs, x, COEF, LL, UL, LS, US))

par(mar = c(2, 6, 0, 0))
plot(COEF, x, pch = 16, xlim = c(-2, 3), ylim = c(0.5, 8), xlab = "Estimate", 
     ylab = NA, yaxt = "n", frame.plot = F)
arrows(x0=COEF-LS, y0=x, x1=US+COEF, y1=x, code=3, angle=90, length=0.1)
abline(v = 0, lty = 2)
axis(2, at = x, labels = labs, las = 2)

####

### End ###