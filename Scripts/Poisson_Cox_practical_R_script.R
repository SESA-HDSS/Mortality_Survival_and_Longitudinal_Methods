###############################################################################
# SURVIVAL ANALYSIS FOR INCIDENCE RATES AND RATE RATIOS (Poisson / Lexis)
# and COX PH MODELS
# NIM, 15/11/2025
#
# This script for Poisson and cox models using survival analysis
#  time-at-risk: each subject contributes person-time between `enter`
# (start of follow-up) and `exit` (end of follow-up).
#
# The failure indicator is:
#    died = 1 event (death)
#    died = 0 censored
#
# The Surv(p_years, died) counting-process syntax or
# The Surv(enter, exit, died) counting-process syntax can be used
###############################################################################

# Initial set-up ----------------------------------------------------------------

setwd("C:/Users/.../poisson_cox_survival_model")

rm(list=ls())


###############################################################################
# 1. Basic incidence rate calculation
###############################################################################

# install.packages(c("haven", "tidyverse", "survival", "epiR", "Epi", "sandwich", "lmtest"))
# library(haven)
# library(tidyverse)
# library(survival)
# library(Epi)
# library(epiR)
# library(sandwich)
# library(lmtest)

pacman::p_load(pacman, haven, tidyverse, ggplot2, epiDisplay, sandwich, sandwich,lmtest,epiR) 


# Import data
mortality <- read_dta("mortality17.dta") |> as_factor()


###############################################################################
# 2. Basic incidence rate calculation
###############################################################################

# Total deaths
sum(mortality$died)

# Total person-years (or person-days depending on unit of enter/exit)
total_p_years <- sum(mortality$exit - mortality$enter)/365.25
total_p_years

# Crude incidence rate
crude_rate <- sum(mortality$died) / as.numeric(total_p_years)
crude_rate

# Per 1000 person-years (common epidemiological scale)
crude_rate_1000 <- crude_rate * 1000
crude_rate_1000


###############################################################################
# 3. Incidence rates by exposure (vimp)
###############################################################################

# Lexis splitting not needed unless rates vary with time;
# simple summarisation:

inc_vimp <- mortality %>%
  group_by(vimp) %>%
  summarise(
    deaths = sum(died),
    p_years = as.numeric(sum(exit - enter)/365.25),
    rate = deaths / p_years,
    rate_1000 = rate * 1000
  )
inc_vimp


###############################################################################
# 4. Poisson regression with survival time as offset
###############################################################################
# (This replaces the earlier glm(died ~ exposure, family=poisson) without time.)

# offset = log(person-time)
mortality$p_years <-as.numeric (mortality$exit - mortality$enter)/365.25
mortality$log_p_years <- log(mortality$p_years)

# Poisson model for rate ratio
prm_died_vimp <- glm(died ~ vimp + offset(log_p_years),
                     family = poisson(),
                     data = mortality)

coeftest(prm_died_vimp)

# Rate ratios and CIs
est <- coef(prm_died_vimp)
pvalue <-round(coeftest(prm_died_vimp)[,4],4)
se  <- round(coeftest(prm_died_vimp)[,2],8)
cbind(round(cbind(IRR = exp(est), LCL = exp(est - 1.96*se), UCL = exp(est + 1.96*se)),3), pvalue=pvalue)

###############################################################################
# 5.  Poisson regression for mfgrp
###############################################################################

prm_died_mfgrp <- glm(died ~ mfgrp + offset(log_p_years),
                      family = poisson(),
                      data = mortality)

coeftest(prm_died_mfgrp)

est <- coef(prm_died_mfgrp)
pvalue <-round(coeftest(prm_died_mfgrp)[,4],4)
se  <- round(coeftest(prm_died_mfgrp)[,2],8)
cbind(round(cbind(IRR = exp(est), LCL = exp(est - 1.96*se), UCL = exp(est + 1.96*se)),3), pvalue=pvalue)


###############################################################################
# 6. Likelihood ratio tests for nested survival Poisson models
###############################################################################

# Null model (intercept only with offset)
prm_null <- glm(died ~ offset(log_p_years),
                family = poisson(),
                data = mortality)

# Compare null vs vimp:
lrtest(prm_null, prm_died_vimp)

# Compare null vs mfgrp:
lrtest(prm_null, prm_died_mfgrp)


###############################################################################
# 7. Poisson model with age
###############################################################################

prm_died_agegrp <- glm(died ~ agegrp + offset(log_p_years),
                       family = poisson(),
                       data = mortality)

coeftest(prm_died_agegrp)

est <- coef(prm_died_agegrp)
pvalue <-round(coeftest(prm_died_agegrp)[,4],4)
se  <- round(coeftest(prm_died_agegrp)[,2],8)
cbind(round(cbind(IRR = exp(est), LCL = exp(est - 1.96*se), UCL = exp(est + 1.96*se)),3), pvalue=pvalue)

# Compare null vs age:
lrtest(prm_null, prm_died_agegrp)

###############################################################################
# 8. Cox model alternative (for hazard ratios)
###############################################################################

### vimp
#cox_vimp <- coxph(SurvObj ~ vimp, data = mortality)
cox_vimp <- coxph(Surv(p_years,died) ~ vimp, data = mortality)
summary(cox_vimp)

coeftest(cox_vimp)
est <- coef(cox_vimp)
pvalue <-round(coeftest(cox_vimp)[,4],4)
se  <- round(coeftest(cox_vimp)[,2],8)
cbind(round(cbind(HR = exp(est), LCL = exp(est - 1.96*se), UCL = exp(est + 1.96*se)),3), pvalue=pvalue)


### mfgrp
#cox_mfgrp <- coxph(SurvObj ~ mfgrp, data = mortality)
cox_mfgrp <- coxph(Surv(p_years,died) ~ mfgrp, data = mortality)
summary(cox_mfgrp)

coeftest(cox_mfgrp)
est <- coef(cox_mfgrp)
pvalue <-round(coeftest(cox_mfgrp)[,4],4)
se  <- round(coeftest(cox_mfgrp)[,2],8)
cbind(round(cbind(HR = exp(est), LCL = exp(est - 1.96*se), UCL = exp(est + 1.96*se)),3), pvalue=pvalue)


# Null model (intercept only with offset)
cox_null <- coxph(Surv(p_years,died)  ~ 1, data = mortality)
summary(cox_null)

# Compare null vs vimp:
lrtest(cox_null, cox_vimp)

# Compare null vs mfgrp:
lrtest(cox_null, cox_mfgrp)


###############################################################################
# 9. Multivariable Poisson model: 
###############################################################################

### (A) vimp + agegrp (confounding analysis)

prm_died_vimp_agegrp <- glm(died ~ vimp + agegrp + offset(log_p_years),
                            family = poisson(),
                            data = mortality)

coeftest(prm_died_vimp_agegrp)

est <- coef(prm_died_vimp_agegrp)
pvalue <-round(coeftest(prm_died_vimp_agegrp)[,4],4)
se  <- round(coeftest(prm_died_vimp_agegrp)[,2],8)
cbind(round(cbind(IRR = exp(est), LCL = exp(est - 1.96*se), UCL = exp(est + 1.96*se)),3), pvalue=pvalue)

# Compare vimp vs vim+age:
lrtest(prm_died_vimp, prm_died_vimp_agegrp)

########################################################
# (B) vimp + agegrp + mfgrp (confounding analysis)

prm_died_vimp_agegrp_mfgrp <- glm(died ~ vimp + agegrp + mfgrp + offset(log_p_years),
                                  family = poisson(),
                                  data = mortality)

coeftest(prm_died_vimp_agegrp_mfgrp)

est <- coef(prm_died_vimp_agegrp_mfgrp)
pvalue <-round(coeftest(prm_died_vimp_agegrp_mfgrp)[,4],4)
se  <- round(coeftest(prm_died_vimp_agegrp_mfgrp)[,2],8)
cbind(round(cbind(IRR = exp(est), LCL = exp(est - 1.96*se), UCL = exp(est + 1.96*se)),3), pvalue=pvalue)

# Compare vimp vs vim+age:
lrtest(prm_died_vimp_agegrp, prm_died_vimp_agegrp_mfgrp)

###############################################################################
# 10. Multivariable Cox model: vimp + agegrp (confounding analysis)
###############################################################################

### (A) vimp + agegrp (confounding analysis)

#cox_vimp_age <- coxph(SurvObj ~ vimp + agegrp, data = mortality)
cox_vimp_age <- coxph(Surv(p_years,died)  ~ vimp + agegrp, data = mortality)
summary(cox_vimp_age)

coeftest(cox_vimp_age)
est <- coef(cox_vimp_age)
pvalue <-round(coeftest(cox_vimp_age)[,4],4)
se  <- round(coeftest(cox_vimp_age)[,2],8)
cbind(round(cbind(HR = exp(est), LCL = exp(est - 1.96*se), UCL = exp(est + 1.96*se)),3), pvalue=pvalue)

# Compare vimp vs agegrp:
lrtest(cox_vimp, cox_vimp_age)

########################################################
### (B) vimp + agegrp + mfgrp (confounding analysis)

#cox_vimp_age_mfgrp <- coxph(SurvObj ~ vimp + agegrp + mfgrp , data = mortality)
cox_vimp_age_mfgrp <- coxph(Surv(p_years,died)  ~ vimp + agegrp + mfgrp , data = mortality)
summary(cox_vimp_age_mfgrp)

coeftest(cox_vimp_age_mfgrp)
est <- coef(cox_vimp_age_mfgrp)
pvalue <-round(coeftest(cox_vimp_age_mfgrp)[,4],4)
se  <- round(coeftest(cox_vimp_age_mfgrp)[,2],8)
cbind(round(cbind(HR = exp(est), LCL = exp(est - 1.96*se), UCL = exp(est + 1.96*se)),3), pvalue=pvalue)

# Compare vimp vs agegrp:
lrtest(cox_vimp_age, cox_vimp_age_mfgrp)


# (Cox results provide hazard ratios; Poisson gives incidence rate ratios.)



###############################################################################
# END OF SURVIVAL ANALYSIS SCRIPT
###############################################################################
