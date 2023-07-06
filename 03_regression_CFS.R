library(sjPlot)
library(lsmeans)
library(lme4)
library(lmerTest)
library(mediation)
library(performance)
library(effectsize)
library(ggeffects)

#########################################################
# DATA LOADING / PREPROCESSING
#########################################################

df <- read.csv("../output/csv/df_concat_R_1sec_CFS.csv")
df$race[df$race == "hispanic"] <- "other"
df$race <- relevel(as.factor(df$race), ref = "caucasian")
df$subj <- as.factor(df$subj)
df$sstats_SPT <- df$sstats_SPT / 60

# Remove useless columns
df <- subset(df, select = -c(sw_ndpac, sw_ndpac_thr, sw_pp, sw_pp_thr_supzero, sw_coinc_spindles))
df <- subset(df, select = -c(ahi_rem, ahi_nrem))

# Define covariates
# NOTE: Many missing values for income and smoking in CFS. 
# This leads to a drastic reduction of sample size if included as covariates.
# >>> + income_cfs + smoking_cfs
covar_lmm <- "+ age + male + race + bmi + hypertension + ahi + sstats_SME + sstats_SPT + (1 | family)"
covar_lmm_nosleep <- "+ age + male + race + bmi + hypertension + (1 | family)"

#########################################################
# FASTING GLUCOSE
#########################################################

# 1) FULLY ADJUSTED MODELS

# Proportion of coupled SO
m <- lmer(paste("fasting_glucose ~ sw_ndpac_prop_supzero", covar_lmm), data=df)
tab_model(m, show.std=T)
plot_model(m, type="pred", show.data=T, terms="sw_ndpac_prop_supzero")

# Calculate marginal means (effects are marginalized, averaged over the levels of factors)
round(quantile(df$sw_ndpac_prop_supzero, c(0, 0.01, .025, .05, .5, .95, .975, 0.99, 1), na.rm=T), 3)
ggeffect(m, "sw_ndpac_prop_supzero [0.77:0.94 by=0.01]")

# Coupling strength
m <- lmer(paste("fasting_glucose ~ sw_ndpac_thr_supzero + diabetes_cfs", covar_lmm), data=df)
tab_model(m, show.std=T)
plot_model(m, type="pred", show.data=T, terms="sw_ndpac_thr_supzero")

round(quantile(df$sw_ndpac_thr_supzero, c(.01, .025, .05, .5, .95, .975, .99), na.rm=T), 3)
ggeffect(m, "sw_ndpac_thr_supzero [0.27:0.38 by=0.02]")

# 2) RANKING OF SLEEP PREDICTORS
# AHI, SPT and efficiency are included as DV and not covariates
feat_sleep <- colnames(df[,match("ai_all", names(df)):match("bp_delta_REM", names(df))])
feat_sleep <- setdiff(feat_sleep, c("diabetes_cfs"))
df_sleep <- data.frame(row.names=feat_sleep)

# mask <- complete.cases(df[feat_sleep])  # Use complete case
# df_complete <- df[mask, ]

for (c in feat_sleep) {
  if (grepl(c, covar_lmm_nosleep)) {
    next
  }
  m <- lmer(paste("fasting_glucose ~", c, covar_lmm_nosleep), data=df)
  sm <- summary(m)$coefficients
  perf <- performance(m,  metrics = "R2")
  std <- standardize_parameters(m, method="basic")
  std <- data.frame(std, row.names = "Parameter")
  df_sleep[c, "beta"] <- sm[c, "Estimate"]
  df_sleep[c, "std_beta"] <- std[c, "Std_Coefficient"]
  df_sleep[c, "df"] <- sm[c, "df"]
  df_sleep[c, "n_valid"] <- sum(!is.na(df[c]))
  df_sleep[c, "r2_cond"] <- perf$R2_conditional
  df_sleep[c, "t"] <- sm[c, "t value"]
  df_sleep[c, "p"] <- sm[c, "Pr(>|t|)"]
}

# Add extra columns
df_sleep$pcorr <- p.adjust(df_sleep$p)
m_empty <- lmer(paste("fasting_glucose ~", covar_lmm_nosleep), data=df)
perf_empty <- performance(m,  metrics = "R2")
df_sleep$r2_gain <- df_sleep$r2_cond - perf_empty$R2_conditional
df_sleep <- df_sleep[order(df_sleep$p), ]
head(df_sleep, 10)

# Save to .csv
write.csv(df_sleep, file="../output/csv/rank_sleep_CFS.csv")

# 3) MEDIATION WITH HRV

detach("package:lmerTest", unload=TRUE)

# Proportion of coupled SO
m1 <- lmer(paste("hrv_rmssd ~ sw_ndpac_prop_supzero", covar_lmm), data=df)
m2 <- lmer(paste("fasting_glucose ~ hrv_rmssd", covar_lmm), data=df)
m3 <- lmer(paste("fasting_glucose ~ hrv_rmssd + sw_ndpac_prop_supzero", covar_lmm), data=df)
tab_model(m1, m2, show.std=T, show.stat=F, show.df=F, collapse.ci=F)
set.seed(42)
summary(mediate(m1, m3, treat="sw_ndpac_prop_supzero", mediator="hrv_rmssd", dropobs=TRUE, sims=10000))

# Coupling strength
m1 <- lmer(paste("hrv_rmssd ~ sw_ndpac_thr_supzero", covar_lmm), data=df)
m2 <- lmer(paste("fasting_glucose ~ hrv_rmssd", covar_lmm), data=df)
m3 <- lmer(paste("fasting_glucose ~ hrv_rmssd + sw_ndpac_thr_supzero", covar_lmm), data=df)
tab_model(m1, m2, show.std=T, show.stat=F, show.df=F, collapse.ci=F)
set.seed(42)
summary(mediate(m1, m3, treat="sw_ndpac_thr_supzero", mediator="hrv_rmssd", dropobs=TRUE, sims=10000))

#########################################################
# HOMA
#########################################################

# 1) FULLY-ADJUSTED MODEL
# Proportion of coupled SO
tab_model(lmer(paste("homa ~ sw_ndpac_prop_supzero", covar_lmm), data=df), show.std=T)

# Coupling strength
tab_model(lmer(paste("homa ~ sw_ndpac_thr_supzero", covar_lmm), data=df), show.std=T)

# 2) RANKING OF SLEEP PREDICTORS
df_sleep <- data.frame(row.names=feat_sleep)
for (c in feat_sleep) {
  if (grepl(c, covar_lmm_nosleep)) {
    next
  }
  m <- lmer(paste("homa ~", c, covar_lmm_nosleep), data=df)
  sm <- summary(m)$coefficients
  perf <- performance(m,  metrics = "R2")
  std <- standardize_parameters(m, method="basic")
  std <- data.frame(std, row.names = "Parameter")
  df_sleep[c, "beta"] <- sm[c, "Estimate"]
  df_sleep[c, "std_beta"] <- std[c, "Std_Coefficient"]
  df_sleep[c, "df"] <- sm[c, "df"]
  df_sleep[c, "n_valid"] <- sum(!is.na(df[c]))
  df_sleep[c, "r2_cond"] <- perf$R2_conditional
  df_sleep[c, "t"] <- sm[c, "t value"]
  df_sleep[c, "p"] <- sm[c, "Pr(>|t|)"]
}

# Add extra columns
df_sleep$pcorr <- p.adjust(df_sleep$p)
m_empty <- lmer(paste("homa ~", covar_lmm_nosleep), data=df)
perf_empty <- performance(m,  metrics = "R2")
df_sleep$r2_gain <- df_sleep$r2_cond - perf_empty$R2_conditional
df_sleep <- df_sleep[order(df_sleep$p), ]
head(df_sleep, 10)

# Save to .csv
write.csv(df_sleep, file="../output/csv/rank_sleep_homa_CFS.csv")

# 3) MEDIATION
m1 <- lmer(paste("hrv_rmssd ~ sw_ndpac_prop_supzero", covar_lmm), data=df)
m2 <- lmer(paste("homa ~ hrv_rmssd", covar_lmm), data=df)
m3 <- lmer(paste("homa ~ hrv_rmssd + sw_ndpac_prop_supzero", covar_lmm), data=df)
tab_model(m1, m2, show.std=T, show.stat=F, show.df=F, collapse.ci=F)
set.seed(42)
summary(mediate(m1, m3, treat="sw_ndpac_prop_supzero", mediator="hrv_rmssd", dropobs=TRUE, sims=10000))