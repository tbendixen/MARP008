######################################################
##### Many Analysts Religion Project: Stage 2  #######
######################################################

### Team 008
### Script written by Theiss Bendixen and Benjamin G. Purzycki
### Contact email:  tb@cas.au.dk

### Last updated: 24th Feb 2021

### Setting up

library(brms)
library(Bayesrel)
library(brmstools) # this package is deprecated but it serves our specific purpose (forest plots)
library(mice)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(tidybayes)
library(tidyr)
library(dplyr)
library(modelr)
library(RColorBrewer)

setwd("") # set working directory

data <- read.csv("MARP_data.csv", sep = ",") # load MARP data
d <- data # store data in object d

############################
### Reliability Analyses ###
############################

### Well-Being
WBcol <- c("wb_psych_1", "wb_psych_2", "wb_psych_3", "wb_psych_4", "wb_psych_5", "wb_psych_6",
           "wb_soc_1", "wb_soc_2", "wb_gen_1", "wb_gen_2") # select variables to test
WBrel <- d[WBcol] # create new dataframe with selected test variables
WBstrel <- strel(WBrel, item.dropped = TRUE) # store the reliability test in an object, include drop-item analysis
summary(WBstrel) # summarize the reliability test

### Religiosity
RELcol <- c("rel_1", "rel_2", "rel_5", "rel_6", "rel_7", "rel_8", "rel_9")
RELrel <- d[RELcol]
RELstrel <- strel(RELrel, item.dropped = TRUE)
summary(RELstrel)

### Cultural Norm
CNcol <- c("cnorm_1", "cnorm_2")
CNrel <- d[CNcol]
CNstrel <- strel(CNrel, item.dropped = FALSE) # no drop-item analysis, since there are only two items
summary(CNstrel)

### Religiosity and Cultural Norms
RCNrel <- cbind(RELrel, CNrel)
RCNstrel <- strel(RCNrel, item.dropped = TRUE)
summary(RCNstrel)

#############################################
### Checking clustering of missing values ###
#############################################

md.pattern(d) # missingness pattern of entire dataset

View(table(data$country, data$age, exclude = FALSE)) # view missingness of age; some clustering -- 11 NAs in Morocco
table(data$country, data$ses, exclude = FALSE) # view missingness of SES; no obvious clustering

#########################
### Prepare variables ###
#########################

### Prepare variables

# Complete cases 
d0 <- data[complete.cases(data$age),] # exclude participants with missing age
d1 <- d0[complete.cases(d0$ses),] # exclude participants with missing SES
d1 <- d1[d1$age >= 18,] # exclude participants with age below 18

# Well-Being
attach(d1)
W <- wb_psych_1 + wb_psych_2 + wb_psych_3 + wb_psych_4 + wb_psych_5 + wb_psych_6 +
 wb_soc_1 + wb_soc_2 + wb_gen_1 + wb_gen_2 # select and sum the well-being score
detach(d1)

d1$W <- W # store the well-being score in a new column

# Religiosity + Cultural Norms, group-mean centered
attach(d1)
R <- rel_1 + rel_2 + rel_5 + rel_6 + rel_7 + rel_8 + rel_9 + cnorm_1 + cnorm_2 # select and sum the religiosity score
detach(d1)

d1$R <- R # store the religiosity score in a new column

d1$R.grp.mean <- ave(d1$R, d1$country) # calculate group-mean of religiosity
d1$R.cf <- d1$R - d1$R.grp.mean # group-mean centering 
d1$R.c <- d1$R.cf / max(d1$R.cf) # transform group-mean centered var onto scale with 1 = maximum possible score on the scale

# Demographic variables, group-mean centered
d1$age.grp.mean <- ave(d1$age, d1$country) # calculate group-mean of age
d1$age.c <- d1$age - d1$age.grp.mean # group-mean centering 
d1$ses.grp.mean <- ave(d1$ses, d1$country) # calculate group-mean of SES
d1$ses.c <- d1$ses - d1$ses.grp.mean # group-mean centering 
d1$education.grp.mean <- ave(d1$education, d1$country) # calculate group-mean of education
d1$education.c <- d1$education - d1$education.grp.mean # group-mean centering

# set attention check to zero = passed, 1 = otherwise
d1$attention_check.r <- as.integer(ifelse(d1$attention_check == 1, 0, 1))

###########################
### Modeling Scenario 2 ###
###########################

# Create list of variables
d_list <- list(
  WB = d1$W,
  Rel = d1$R.c,
  Age = d1$age.c,
  Gender = d1$gender,
  SES = d1$ses.c,
  Edu = d1$education.c,
  AC = d1$attention_check.r,
  Country = d1$country
)

# Setting priors
aprior <- set_prior("normal(30,2)", class = "Intercept")
bprior <- set_prior("normal(0,1)", class = "b")
lkjprior <- set_prior("lkj_corr_cholesky(2)", class = "cor")
sdprior <-set_prior("exponential(1)", class = "sd")
sigmaprior <-set_prior("exponential(1)", class = "sigma")

# Fit m7
m7 <- brm(WB ~ Rel + (Rel | Country) 
          + Gender + Age + SES + Edu + AC,
          data = d_list,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# m7 with brms default priors -- qualitatively similar to m7
m0 <- brm(WB ~ Rel + (Rel | Country) 
          + Gender + Age + SES + Edu + AC,
          data = d_list, 
          family = gaussian(),
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# Quick summaries
summary(m7)
ranef(m7)
conditional_effects(m7)

# Fit m8
m8 <- brm(WB ~ Rel + (Rel | Country) + AC,
          data = d_list,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# Fit m9
m9 <- brm(WB ~ Rel + (1 | Country) 
          + Gender + Age + SES + Edu + AC,
          data = d_list,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# Fit the model
m10 <- brm(WB ~ Rel + (1 | Country) + AC,
           data = d_list,
           family = gaussian(), 
           prior = aprior + bprior + sdprior + sigmaprior,
           sample_prior = TRUE,
           control = list(adapt_delta = 0.99),
           chains = 4, iter = 2000)

########################
### Model Comparison ###
########################

# Table 1

(loo_m7 <- loo(m7)) # approximate leave-one-out cross-validation
(loo_m8 <- loo(m8)) # approximate leave-one-out cross-validation
(loo_m9 <- loo(m9)) # approximate leave-one-out cross-validation
(loo_m10 <- loo(m10)) # approximate leave-one-out cross-validation

loo_compare(loo_m7, loo_m8, loo_m9, loo_m10) # compare models
model_weights(m7, m8, m9, m10, weights = "loo") # compare models

########################
### Plotting Model 7 ###
########################

# Global posterior predictive check (not in the write-up, but if interested, unblock and run next line)
#ppc_m7.global <- pp_check(m7, nsamples = 20) + ylim(0,.067) + xlim(0, 60) + theme(legend.position="none")

### Figure 1: posterior predictive check (ppc) across countries
y <- d1$W # observed (raw) well-being scores
yrep <- brms::posterior_predict(m7) # extract posterior predictions

ppc_m7.country <- ppc_dens_overlay_grouped(y, yrep[1:20,], # overlay 20 predictive draws on observed scores grouped by countries
                                           group = d_list$Country) + theme(legend.position="none") + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

(fig1 <- ppc_m7.country + theme(text = element_text(size=40))) + ylim (0, 0.1) + xlim(5,60) + geom_vline(xintercept=35.49, linetype="dotted", color = "darkblue") # with dotted lines at grand mean

### Figure 2: forest plots of across countries
m7_inter <- forest(m7, pars = c("Intercept"))
m7_rel <- forest(m7, pars = c("Rel")) + geom_vline(xintercept=c(0), linetype="dotted")
(m7_inter) + (m7_rel) + patchwork::plot_layout(ncol = 1, nrow = 2)

### Figure 3: scatterplot across countries with posterior predictive draws
scatdf <- as.data.frame(d_list) # convert d_list to dataframe

(scatplotfd = scatdf %>%
    dplyr::group_by(Country) %>% # group dataframe by countries
    modelr::data_grid(Rel = seq_range(-1:1, by = 0.1), Gender="man", SES=0, Edu=0, AC=0, Age=0) %>% # values to be predicted from
    add_linpred_draws(newdata=., model=m7, n = 100) %>% # add 100 linear posterior predictions to dataframe (for each country)
    ggplot(aes(x = Rel, y = WB)) + # set up panels
    geom_jitter(data = scatdf, width = 0.01, shape = 1, color = "black") + # add raw data points with jitter
    geom_line(aes(y = .value, group = paste(Country, .draw)), alpha = 1/20, color = "#08519C") + # add predictive lines
    facet_wrap(~Country, nrow = 4) + # break up the plot for each country
    theme_bw() + # remove grey panel background
    theme(strip.background = element_blank(), # remove grey background for title
          axis.title = element_text( size = 30, face = "bold"), # change size and bold font of global axis titles (well-being and religiosity)
          axis.text = element_text( size = 15), # axes text
          strip.text = element_text(size = 40)) + # country text size
    scale_y_continuous(limits=c(10,50), name = "Well-Being") + # scale and text on global Y-axis
    scale_x_continuous(name = "Religiosity (country-mean centered)")) # scale and text on global X-axis

########################################################################################
### Modelling Scenario 2 again, with absolute deviance from country-mean Religiosity ###
########################################################################################

# data list
d_list2 <- list(
  WB = d1$W,
  Rel2 = abs(d1$R.c), # absolute deviance from group-mean centered Religiosity
  Age = d1$age.c,
  Gender = d1$gender,
  SES = d1$ses.c,
  Edu = d1$education.c,
  AC = d1$attention_check.r,
  Country = d1$country
)

# Setting priors
aprior <- set_prior("normal(30,2)", class = "Intercept")
bprior <- set_prior("normal(0,1)", class = "b")
lkjprior <- set_prior("lkj_corr_cholesky(2)", class = "cor")
sdprior <-set_prior("exponential(1)", class = "sd")
sigmaprior <-set_prior("exponential(1)", class = "sigma")

# Fit m11
m11 <- brm(WB ~ Rel2 + (Rel2 | Country) 
          + Gender + Age + SES + Edu + AC,
          data = d_list2,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# Predicted draws, Model 11 (does not feature in report -- for notation, see above)
scatdf2 <- as.data.frame(d_list2)

(scatplotfd2 = scatdf2 %>%
    dplyr::group_by(Country) %>%
    modelr::data_grid(Rel2 = seq_range(0:1, by = 0.1), Gender="man", SES=0, Edu=0, AC=0, Age=0) %>%
    add_linpred_draws(newdata=., model=m11, n = 100) %>%
    ggplot(aes(x = Rel2, y = WB)) +
    geom_jitter(data = scatdf2, width = 0.01, shape = 1, color = "black") +
    geom_line(aes(y = .value, group = paste(Country, .draw)), alpha = 1/20, color = "#08519C") +
    facet_wrap(~Country, nrow = 4) + 
    theme_bw() +
    theme(strip.background = element_blank(),
          axis.title = element_text( size = 30, face = "bold"),
          axis.text = element_text( size = 15),
          strip.text = element_text(size = 40)) + 
    scale_y_continuous(limits=c(10,50), name = "Well-Being") + 
    scale_x_continuous(limits=c(0,1), name = "Religiosity (absolute deviance from country-mean)"))

# Quick summaries
summary(m11)
ranef(m11)
conditional_effects(m11)

### Figure 4: forest plots of Model 11
m11_inter <- forest(m11, pars = c("Intercept"))
m11_rel <- forest(m11, pars = c("Rel2")) + geom_vline(xintercept=c(0), linetype="dotted")
(m11_inter) + (m11_rel) + patchwork::plot_layout(ncol = 1, nrow = 2)

# Fit m12
m12 <- brm(WB ~ Rel2 + (Rel2 | Country) + AC,
          data = d_list2,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# Quick summaries
summary(m12)
ranef(m12)
conditional_effects(m12)
pp_check(m12)
forest(m12)

# Fit m13
m13 <- brm(WB ~ Rel2 + (1 | Country) 
          + Gender + Age + SES + Edu + AC,
          data = d_list2,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# Fit m14
m14 <- brm(WB ~ Rel2 + (1 | Country) + AC,
           data = d_list2,
           family = gaussian(), 
           prior = aprior + bprior + sdprior + sigmaprior,
           sample_prior = TRUE,
           control = list(adapt_delta = 0.99),
           chains = 4, iter = 4000)


###############################################
### EXTRAs:                                 ###
### - truncated normal outcome distribution ###
### - missing values imputed                ###
###############################################

### Modelling Model 7 with truncated normal -- samples very slowly, so blocked out ###

# Fit m7 with a truncated normal distribution
# m7.t <- brm(WB | trunc(lb = 10, ub = 50) ~ Rel + (Rel | Country) 
#          + Gender + Age + SES + Edu + AC,
#          data = d_list,
#          family = gaussian(), 
#          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
#          sample_prior = TRUE,
#          control = list(adapt_delta = 0.99),
#          chains = 4, iter = 2000)

# Quick summaries
#summary(m7.t)
#ranef(m7.t)
#conditional_effects(m7.t)
#pp_check(m7.t)
#forest(m7.t)

### Missing age values imputation ###
data.imp <- subset(data, age >= 18 | is.na(age)) # exclude ages below 18
d.i <- as.data.frame(is.na(data.imp)) # new data frame to do imputation

# set irrelevant vars with NAs to FALSE (so that they are not imputed)
d.i$wb_soc_3 <- FALSE
d.i$ethnicity <- FALSE
d.i$denomination <- FALSE
d.i$ses <- FALSE

ini <- mice(data.imp, maxit = 0) # create imputation matrix
pred <- ini$pred # create imputation matrix (which variables are used to impute)
pred["subject",] <- 0 # do not use subject ID to impute
pred["age",] <- 0 # do not use age to impute, since we are imputing missing age values

# create 20 imputed datasets using "predictive mean matching" (pmm)
d.imp <- mice(data.imp, pred = pred, where = d.i, meth = "pmm", m = 20, maxit = 20)

# diagnostics: check convergence of imputation
plot(d.imp, c("age"))

# convenience functions to compare the imputed values with observed data
densityplot(d.imp, ~ age) # red lines = imputed values, blue lines = observed
xyplot(d.imp, age ~ country) # red points = imputed values, blue points = observed
summary(complete(d.imp)) # summarize imputed data descriptives
summary(data) # summarize original data descriptives
hist(data.imp[data.imp$country == "Morocco",]$age) # zoom in on Morocco before imputation
hist(complete(d.imp)[complete(d.imp)$country == "Morocco",]$age) # zoom in on Morocco after imputation
hist(data.imp$age) # histogram of age in original data
hist(complete(d.imp)$age) # histogram of age in imputed data

# extract an imputed data set
df.imp <- complete(d.imp)
df.imp <- df.imp[complete.cases(df.imp$ses),] # exclude the few non-clustered missing values in SES

# The following code mirrors the preparation of variables for Models 7-10 (see above), 
# but now we use a randomly drawn imputed dataset.
 
# Demographic variables, group-mean centered
df.imp$age.grp.mean <- ave(df.imp$age, df.imp$country)
df.imp$age.c <- df.imp$age - df.imp$age.grp.mean
df.imp$ses.grp.mean <- ave(df.imp$ses, df.imp$country)
df.imp$ses.c <- df.imp$ses - df.imp$ses.grp.mean 
df.imp$education.grp.mean <- ave(df.imp$education, df.imp$country)
df.imp$education.c <- df.imp$education - df.imp$education.grp.mean

# Well-Being
attach(df.imp)
W.imp <- wb_psych_1 + wb_psych_2 + wb_psych_3 + wb_psych_4 + wb_psych_5 + wb_psych_6 +
  wb_soc_1 + wb_soc_2 + wb_gen_1 + wb_gen_2
detach(df.imp)

df.imp$W <- W.imp

# Religiosity + Cultural Norms, group-mean centered
attach(df.imp)
R.imp <- rel_1 + rel_2 + rel_5 + rel_6 + rel_7 + rel_8 + rel_9 + cnorm_1 + cnorm_2
detach(df.imp)

df.imp$R <- R.imp

df.imp$R.grp.mean <- ave(df.imp$R, df.imp$country) # group-mean centering
df.imp$R.cf <- df.imp$R - df.imp$R.grp.mean # group-mean centering
df.imp$R.c <- df.imp$R.cf / max(df.imp$R.cf) # transform onto -1 - 1 scale

# Attention check
df.imp$attention_check.r <- as.integer(ifelse(df.imp$attention_check == 1, 0, 1)) # set attention check to zero = passed, 1 = otherwise

# data list
d_list.imp <- list(
  WB = df.imp$W,
  Rel = df.imp$R.c,
  Age = df.imp$age.c,
  Gender = df.imp$gender,
  SES = df.imp$ses.c,
  Edu = df.imp$education.c,
  AC = df.imp$attention_check.r,
  Country = df.imp$country
)

# Fit m7
m7.imp <- brm(WB ~ Rel + (Rel | Country) 
            + Gender + Age + SES + Edu + AC,
            data = d_list.imp,
            family = gaussian(), 
            prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
            sample_prior = TRUE,
            control = list(adapt_delta = 0.99),
            chains = 4, iter = 2000)

# Quick summaries
summary(m7.imp) # results are very similar to Model 7
ranef(m7.imp)
conditional_effects(m7.imp)
pp_check(m7.imp)
forest(m7.imp)

################################
### Appendix A - Scenario 1  ###
################################

# Religiosity, without Cultural Norms items
attach(d1)
R_app <- rel_1 + rel_2 + rel_5 + rel_6 + rel_7 + rel_8 + rel_9
detach(d1)

d1$R_app <- R_app
d1$R_app <- d1$R_app / 7 # transform onto 0 - 1 scale, with 1 = maximum possible score on the scale (7)

# Cultural Norms
attach(d1)
C <- cnorm_1 + cnorm_2
detach(d1)

d1$C <- C
d1$C <- d1$C / 2 # transform onto 0-1 scale, with 1 = maximum possible score on the scale (2)

# data list
d_list_app <- list(
  WB = d1$W,
  Rel = d1$R_app,
  Age = d1$age.c,
  Gender = d1$gender,
  SES = d1$ses.c,
  Edu = d1$education.c,
  AC = d1$attention_check.r,
  Country = d1$country,
  Norm = d1$C
)

# Fit the model
m1 <- brm(WB ~ Rel*Norm + (Rel*Norm | Country) 
          + Gender + Age + SES + Edu + AC,
          data = d_list_app,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model Comparison (Table A.3)
(loo_m1 <- loo(m1))
(loo_m2 <- loo(m2))
(loo_m3 <- loo(m3))
(loo_m4 <- loo(m4))
(loo_m5 <- loo(m5))
(loo_m5.extra <- loo(m5.extra))
(loo_m5.extra.norm <- loo(m5.extra.norm))
(loo_m5.extra.rel <- loo(m5.extra.rel))
(loo_m6 <- loo(m6))

loo_compare(loo_m1, loo_m2, loo_m3, loo_m4, loo_m5, loo_m5.extra, loo_m5.extra.norm, loo_m5.extra.rel, loo_m6)
round(model_weights(m1, m2, m3, m4, m5, m5.extra,m5.extra.norm,m5.extra.rel, m6, weights = "loo"), 3)

### Figure A.5: Predicted draws, Model 1, Cultural Norms on X-axis
scatdf_app <- as.data.frame(d_list_app)

# add posterior predictive draws for the interaction, when Religiosity and the interaction term equals 1
m1xnorm = scatdf_app %>% 
  dplyr::group_by(Country) %>%
  modelr::data_grid(Rel = 1, Norm = seq_range(0:1, by = .5), "Rel:Norm" = 1, Gender="man", SES=0, Edu=0, AC=0, Age=0) %>%
  add_linpred_draws(newdata=., model = m1, value = ".int", n = 100)

# add posterior predictive draws for the interaction, when Religiosity and the interaction term equals 0,
# and then combine both set of interaction lines in the same plot. 
(scatplotm1norm = scatdf_app %>%
    dplyr::group_by(Country) %>%
    modelr::data_grid(Rel = 0, Norm = seq_range(0:1, by = .5), "Rel:Norm" = 0, Gender="man", SES=0, Edu=0, AC=0, Age=0) %>%
    add_linpred_draws(newdata=., model = m1, value = ".value", n = 100) %>%
    ggplot(aes(x = Norm, y = WB, color = factor(ifelse(Rel > 0.5, 1, 0)))) + # color the raw points according to value of Religiosity
    scale_colour_manual(values = c("black", "#08519C")) + # specify colors (blue points = Religiosity > .5, black otherwise
    geom_jitter(data = scatdf_app, width = 0.01, shape = 1) +
    geom_line(aes(y = .value, group = paste(Country, .draw)), alpha = 1/20, color = "black") + # black predictive lines for low religiosity
    geom_line(aes(y = m1xnorm$.int, group = paste(m1xnorm$Country, m1xnorm$.draw)), alpha = 1/20, color = "#08519C") + # blue predictive lines for high religiosity
    facet_wrap(~Country, nrow = 4) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "none", # remove legend
          axis.title = element_text( size = 30, face = "bold"),
          axis.text = element_text(  size = 15),
          strip.text = element_text(size = 30)) + 
    scale_y_continuous(limits=c(10,50), name = "Well-Being") + 
    scale_x_continuous(breaks=c(0,.5,1), name = "Cultural Norm"))

### Predicted draws, Model 1, Religiosity on X-axis (does not feature in the report)

# add posterior predictive draws for the interaction, when Cultural Norms and the interaction term equals 1
m1xrel = scatdf_app %>%
  dplyr::group_by(Country) %>%
  modelr::data_grid(Rel = seq_range(0:1, by = .5), Norm = 1, "Rel:Norm" = 1, Gender="man", SES=0, Edu=0, AC=0, Age=0) %>%
  add_linpred_draws(newdata=., model = m1, value = ".int", n = 100)

# add posterior predictive draws for the interaction, when Cultural Norms and the interaction term equals 0,
# and then combine both set of interaction lines in the same plot. 
(scatplotm1 = scatdf_app %>%
    dplyr::group_by(Country) %>%
    modelr::data_grid(Rel = seq_range(0:1, by = 0.5), Norm = 0, "Rel:Norm" = 0, Gender="man", SES=0, Edu=0, AC=0, Age=0) %>%
    add_linpred_draws(newdata=., model = m1, value = ".value", n = 100) %>%
    ggplot(aes(x = Rel, y = WB, color = factor(ifelse(Norm > 0.5, 1, 0)))) + # color the raw points according to value of Cultural Norms
    scale_colour_manual(values = c("black", "#08519C")) + # specify colors (blue points = Cultural Norms > .5, black otherwise
    geom_jitter(data = scatdf_app, width = 0.01, shape = 1) +
    geom_line(aes(y = .value, group = paste(Country, .draw)), alpha = 1/20, color = "black") + # blue predictive lines for low "Cultural Norm"-ness
    geom_line(aes(y = m1xrel$.int, group = paste(m1xrel$Country, m1xrel$.draw)), alpha = 1/20, color = "#08519C") + # blue predictive lines for high "Cultural Norm"-ness
    facet_wrap(~Country, nrow = 4) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "none",
          axis.title = element_text( size = 30, face = "bold"),
          axis.text = element_text( size = 15),
          strip.text = element_text(size = 30)) +
    scale_y_continuous(limits=c(10,50), name = "Well-Being") + 
    scale_x_continuous(breaks=c(0,.5,1), name = "Religiosity"))

# Combine the two interaction plots (does not feature in the report)
(scatplotm1) + (scatplotm1norm) + patchwork::plot_layout(ncol = 2, nrow = 1)

### Figure A.6: forest plots of Model 1
m1_inter <- forest(m1, pars = c("Intercept"))
m1_rel <- forest(m1, pars = c("Rel")) + geom_vline(xintercept=c(0), linetype="dotted") + xlim(-4,8)
m1_norm <- forest(m1, pars = c("Norm")) + geom_vline(xintercept=c(0), linetype="dotted") + xlim(-4,8)
m1_relnorm <- forest(m1, pars = c("Rel:Norm")) + geom_vline(xintercept=c(0), linetype="dotted") + xlim(-4,8)

(m1_inter) + (m1_rel) + (m1_norm) + (m1_relnorm) + patchwork::plot_layout(ncol = 2, nrow = 2)

### Model 2 ###

# Fit the model
m2 <- brm(WB ~ Rel*Norm + (Rel*Norm | Country) + AC,
          data = d_list_app,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 3 ###

# Fit the model
m3 <- brm(WB ~ Rel*Norm + (1 | Country) 
          + Gender + Age + SES + Edu + AC,
          data = d_list_app,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 4 ###

# Fit the model
m4 <- brm(WB ~ Rel*Norm + (1 | Country) + AC,
          data = d_list_app,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

### Model 5 ###

# Fit the model
m5 <- brm(WB ~ Rel + Norm + (1 | Country) 
          + Gender + Age + SES + Edu + AC,
          data = d_list_app,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# Fit the model
m5.extra <- brm(WB ~ Rel + Norm + (Rel + Norm | Country) 
          + Gender + Age + SES + Edu + AC,
          data = d_list_app,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# Fit model with only Norm, to illustrate multicollinearity
m5.extra.norm <- brm(WB ~ Norm + (Norm | Country) 
                     + Gender + Age + SES + Edu + AC,
                     data = d_list_app,
                     family = gaussian(), 
                     prior = aprior + bprior + sdprior + sigmaprior,
                     sample_prior = TRUE,
                     control = list(adapt_delta = 0.99),
                     chains = 4, iter = 2000)

# Fit model with only Rel, to illustrate multicollinearity
m5.extra.rel <- brm(WB ~ Rel + (Rel | Country) 
                     + Gender + Age + SES + Edu + AC,
                     data = d_list_app,
                     family = gaussian(), 
                     prior = aprior + bprior + sdprior + sigmaprior,
                     sample_prior = TRUE,
                     control = list(adapt_delta = 0.99),
                     chains = 4, iter = 2000)

### Model 6 ###

# Fit the model
m6 <- brm(WB ~ Rel + Norm + (1 | Country) + AC,
          data = d_list_app,
          family = gaussian(), 
          prior = aprior + bprior + sdprior + sigmaprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

# For summary and plotting methods, see above.

##############################################
### Appendix B - Social Desirability Model ###
##############################################

# Create list of variables
d_list.phys <- list(
  phys_WB = (d1$wb_phys_1 + d1$wb_phys_2 + d1$wb_phys_3 + d1$wb_phys_4 + d1$wb_phys_5 + d1$wb_phys_6 + d1$wb_phys_7)*(50/35), # getting the summation on the same scale as the main models (10 to 50)
  Rel = as.integer((d1$rel_1*6)+1), # getting the item back in original categorical form, in order to model it monotonically.
  Age = d1$age.c,
  Gender = d1$gender,
  SES = d1$ses.c,
  Edu = d1$education.c,
  AC = d1$attention_check.r,
  Country = d1$country
)

# Setting priors
aprior <- set_prior("normal(30,2)", class = "Intercept")
bprior <- set_prior("normal(0,1)", class = "b")
lkjprior <- set_prior("lkj_corr_cholesky(2)", class = "cor")
sdprior <-set_prior("exponential(1)", class = "sd")
sigmaprior <-set_prior("exponential(1)", class = "sigma")
diriprior <- set_prior("dirichlet(2)", class = "simo", coef= "moRel1") # prior for monotonic function

# Fit m7.phys
m7.phys <- brm(phys_WB ~ mo(Rel) + (mo(Rel) | Country) # modelling Religiosity as monotonic function, mo(Rel)
          + Gender + Age + SES + Edu + AC,
          data = d_list.phys,
          family = gaussian(), 
          prior = aprior + bprior + lkjprior + sdprior + sigmaprior + diriprior,
          sample_prior = TRUE,
          control = list(adapt_delta = 0.99),
          chains = 4, iter = 2000)

summary(m7.phys)

### Figure B.7: forest plot of social desirability model
m7.phys_inter <- forest(m7.phys, pars = c("Intercept"))
m7.phys_rel <- forest(m7.phys, pars = c("moRel")) + geom_vline(xintercept=c(0), linetype="dotted")

(m7.phys_inter) + (m7.phys_rel) + patchwork::plot_layout(ncol = 1, nrow = 2)

###########
### END ###
###########