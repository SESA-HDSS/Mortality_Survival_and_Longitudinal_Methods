

## Still birth


pacman::p_load(tidyverse, rio, ggplot2, data.table, Hmisc, dplyr, lme4, survival, lubridate, patchwork, margins, broom, car, DescTools, ez, emmeans) ## all at once





setwd("C:\\Users\\user pc\\Desktop\\Multilevel model HDSS\\stata")



still_birth<- import("C:\\Users\\user pc\\Desktop\\Multilevel model HDSS\\stata\\Pregnancy_outcomeFinal.dta")


## change to Factor 


still_birth <- still_birth %>%
  mutate(stillbirth = factor(
    stillbirth,
    levels = 0:1,
    labels = c("No", "Yes")
  ))

table(still_birth$stillbirth)



still_birth <- still_birth %>%
  mutate(sba = factor(
    BIRTH_ATTE,
    levels = 1:9,
    labels = c(
      "Doctor",
      "Nurse",
      "Midwife",
      "Health Officer",
      "Health Extension Worker",
      "Traditional Birth Attendant",
      "Family/Friend",
      "Other",
      "No one"
    )
  ))

## Recode 


still_birth <- still_birth %>%
  mutate(atten = case_when(
    BIRTH_ATTE %in% 1:5 ~ 1,
    BIRTH_ATTE %in% 6:9 ~ 0
  ),
  atten = factor(atten, levels = c(0,1),
                 labels = c(" skilled", "Not Skilled")))


str(still_birth$stillbirth)
table(still_birth$stillbirth)
table(still_birth$atten)



#### simulat age 

set.seed(123)   # for reproducibility

still_birth <- still_birth %>%
  mutate(
    age = runif(n(), min = 15, max = 49)   # continuous age 15–49
  )




### MLM

# Empty model 

null_mod <- glmer(
  stillbirth ~ 1 +
    (1 | LOCATIONID) +
    (1 | INDIVIDID),
  data = still_birth,
  family = binomial(link="logit"),
  control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(null_mod)



###################################### Variance and ICC


vc <- as.data.frame(VarCorr(null_mod))
var_loc   <- vc$vcov[vc$grp == "LOCATIONID"]   # Level 3
var_woman <- vc$vcov[vc$grp == "INDIVIDID"]    # Level 2
var_resid <- (pi^2)/3                          # Level 1 logistic residual

total_var <- var_loc + var_woman + var_resid

ICC_loc   <- var_loc   / total_var*100
ICC_loc
ICC_woman <- var_woman / total_var*100
ICC_woman
ICC_resid <- var_resid / total_var*100
ICC_resid




#### ML with factors



fullmodel <- glmer(
  stillbirth ~ atten + 
    (1 | LOCATIONID) + 
    (1 | INDIVIDID),
  data = still_birth,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5),
    calc.derivs = FALSE
  )
)

summary(fullmodel) 



############## Variance and ICC


vc2 <- as.data.frame(VarCorr(fullmodel))
var_loc2   <- vc2$vcov[vc2$grp == "LOCATIONID"]   # Level 3
var_woman2 <- vc2$vcov[vc2$grp == "INDIVIDID"]    # Level 2
var_resid2 <- (pi^2)/3                          # Level 1 logistic residual

total_var_full <- var_loc2 + var_woman2 + var_resid2

ICC_loc_full   <- var_loc2   / total_var_full*100
ICC_loc_full
ICC_woman_full <- var_woman2 / total_var_full*100
ICC_woman_full 
ICC_resid_full <- var_resid2 / total_var_full*100
ICC_resid_full


######### Proportional change in variance (PCV)



pcv3 <- (113 - 100) / 113*100
pcv3




#############       centering 


still_birth <- still_birth %>%
  mutate(
    age_c = age - mean(age, na.rm = TRUE)
  )



fullmodel3 <- glmer(
  stillbirth ~ age_c +
    (1 | LOCATIONID) + 
    (1 | INDIVIDID),
  data = still_birth,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5),
    calc.derivs = FALSE
  )
)

summary(fullmodel3) 




vc <- as.data.frame(VarCorr(fullmodel))
var_loc   <- vc$vcov[vc$grp == "LOCATIONID"]   # Level 3
var_woman <- vc$vcov[vc$grp == "INDIVIDID"]    # Level 2
var_resid <- (pi^2)/3                          # Level 1 logistic residual

total_var <- var_loc + var_woman + var_resid

ICC_loc   <- var_loc   / total_var
ICC_woman <- var_woman / total_var
ICC_resid <- var_resid / total_var








