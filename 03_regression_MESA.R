library(sjPlot)
library(glmnet)
library(lsmeans)
library(mediation)
library(performance)

#########################################################
# DATA LOADING / PREPROCESSING
#########################################################

df <- read.csv("../output/csv/df_concat_R_1sec_MESA.csv")
df$race <- relevel(as.factor(df$race), ref = "caucasian")
df$diabetes <- relevel(as.factor(df$diabetes), ref = "normal")
# levels(df$diabetes) <- list("0"=c("normal", "impaired"), "1"=c("treated", "untreated"))
table(df$diabetes)

# Remove useless columns
df <- subset(df, select = -c(sw_ndpac, sw_pp, sw_pp_thr_supzero, sw_coinc_spindles))

# Define covariates
covar_nosleep <- "+ age + male + race + bmi + hypertension" # + income + smoking"
covar <- paste(covar_nosleep, "+ ahi + sstats_SME + sstats_SPT")

#########################################################
# STANDARD REGRESSIONS
#########################################################

# Fasting glucose
tab_model(lm(paste("fasting_glucose ~ sw_ndpac_prop_supzero", covar), data=df), show.std=T)
tab_model(lm(paste("fasting_glucose ~ sw_ndpac_thr_supzero", covar), data=df), show.std=T)

#########################################################
# ASSOCIATION / MEDIATION WITH HRV
#########################################################

# Fasting glucose
m1 <- lm(paste("hrv_rmssd ~ sw_ndpac_prop_supzero", covar), data=df)
m2 <- lm(paste("fasting_glucose ~ hrv_rmssd", covar), data=df)
m3 <- lm(paste("fasting_glucose ~ hrv_rmssd + sw_ndpac_prop_supzero", covar), data=df)
set.seed(42)
summary(mediate(m1, m3, treat="sw_ndpac_prop_supzero", mediator="hrv_rmssd", dropobs=TRUE, sims=10000))
tab_model(m1, m2, show.std=T, show.stat=F, show.df=F, collapse.ci=F)

# Coupling strength
m1 <- lm(paste("hrv_rmssd ~ sw_ndpac_thr_supzero", covar), data=df)
m2 <- lm(paste("fasting_glucose ~ hrv_rmssd", covar), data=df)
m3 <- lm(paste("fasting_glucose ~ hrv_rmssd + sw_ndpac_thr_supzero", covar), data=df)
tab_model(m1, m2, show.std=T, show.stat=F, show.df=F, collapse.ci=F)
set.seed(42)
summary(mediate(m1, m3, treat="sw_ndpac_thr_supzero", mediator="hrv_rmssd", dropobs=TRUE, sims=10000))