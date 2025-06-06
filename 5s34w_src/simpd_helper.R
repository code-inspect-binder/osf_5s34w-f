# R script to acccompany Knowles & Badh 2023: Impact of face masks on speech in Parkinson's disease: Effect of clear and loud speech styles
# This R script will be sourced in the .Rmd manuscript file.

# Load libraries ----
library(tidyverse)

# Functions ----
# Coefs with EffSizes ----
make_coefs_with_effsize <- function(m, measure = "measure"){
  tmp_coefs <- summary(m)$coefficients %>% as.data.frame()
  tmp_vars <- summary(m)$varcor %>% as.data.frame()
  totvar <- sum(tmp_vars$vcov)
  output <- tmp_coefs %>% 
    mutate(d = abs(Estimate/sqrt(totvar))) %>%
    mutate(measure = measure) %>%
    select(-df) %>%
    mutate_if(is.numeric,round,3) %>%
    rownames_to_column(var = "Contrast")
  
  output
}

contrast_terms_cleaned <- c(
  "cond_speech" = "",
  "cond_mask" = "",
  "group" = "",
  "Habit vs Altered" = "h-cl",
  "Habit vs Clear/Loud" = "h-cl",
  "Clear vs Loud" = "c-l",
  "NM vs Masks" = "nm-ms",
  "SM vs KN" = "sm-kn",
  "genderW vs M" = "gen",
  "intake_genderm" = "gen",
  "OC vs PD" = "gr"
)

# Make df of beta-values for in-text usage
make_b <- function(m, terms = contrast_terms_cleaned){
  coefs_cleaned <- summary(m)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "Contrast") %>%
    mutate(Contrast = str_replace_all(Contrast, contrast_terms_cleaned)) %>%
    mutate_if(is.numeric, round, 3)
  
  coefsB <- column_to_rownames(coefs_cleaned, "Contrast")["Estimate"]
  
  coefsB
}

# Make df of p-values for in-text usage
make_p <- function(m, terms = contrast_terms_cleaned){
  coefs_cleaned <- summary(m)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "Contrast") %>%
    mutate(Contrast = str_replace_all(Contrast, contrast_terms_cleaned)) %>%
    mutate_if(is.numeric, round, 3)
  
  coefsP <- column_to_rownames(coefs_cleaned, "Contrast")["Pr(>|t|)"]
  
  coefsP
}

# Load in cleaned data ----
df <- readr::read_csv("simpd_data_cleaned.csv")


# Models ----
# .. Contrasts ----
# (Reverse) Helmert for cond speech, cond mask
# Treatment for speech_mask
# Sum for group

# Factor variables
df <- df %>%
  mutate(
    group = factor(group), # oc, pd
    gender = factor(intake_gender), # f, m 
    cond_mask = factor(cond_mask,
                       levels = c("nm","sm","kn")),
    cond_speech = factor(cond_speech,
                         levels = c("habitual","clear","loud"))
  ) %>%
  mutate(speech_mask = paste(cond_speech,cond_mask,sep="_"))

levels(df$group) # oc, pd
levels(df$gender) # f, m
levels(df$cond_speech) # habitual, clear, loud
levels(df$cond_mask) # nm, sm, kn

# Reference level: habitual
contrasts(df$cond_speech) <- matrix(c(2/3, -1/3, -1/3, # habitual vs [clear & loud]
                                      0, 1/2, -1/2), # clear vs loud
                                    ncol = 2)

colnames(contrasts(df$cond_speech)) <- c("Habit vs Clear/Loud", "Clear vs Loud")
solve(cbind(1,contrasts(df$cond_speech))) # Habit + vs C/L -; Clear + vs Loud -

contrasts(df$cond_mask) <- matrix(c(2/3, -1/3, -1/3, # nm vs [sm & kn]
                                    0, 1/2, -1/2), # sm vs kn
                                  ncol = 2)
colnames(contrasts(df$cond_mask)) <- c("NM vs Masks", "SM vs KN")
solve(cbind(1,contrasts(df$cond_mask))) # NM + vs Masks - ; SM + vs KN -

contrasts(df$group) <- contr.sum(2) # oc = +1; pd = -1
colnames(contrasts(df$group)) <- c("OC vs PD")
solve(cbind(1,contrasts(df$group)))

contrasts(df$gender) <- contr.sum(2) # f = +1; m = -1
colnames(contrasts(df$gender)) <- c("W vs M")
solve(cbind(1,contrasts(df$gender)))


# ..Slopes ----
# Speech condition contrasts
contrasts(df$cond_speech)
model.matrix(~cond_speech, df) %>% head()
df$HvCL <- model.matrix(~cond_speech, df)[,2]
df$CvL <- model.matrix(~cond_speech, df)[,3]

# Mask condition contrasts
contrasts(df$cond_mask)
model.matrix(~cond_mask, df) %>% head()
df$NMvM <- model.matrix(~cond_mask, df)[,2]
df$SMvKN <- model.matrix(~cond_mask, df)[,3]


#..RQ1: habitual speech ----
# Fixed effects: group*mask + gender
# Random effects: participant, utterance intercepts
# Non-convergence with random slopes for all models; proceed w/ random intercepts only
df_habitual <- df %>% 
  filter(cond_speech == "habitual",
         channel == "ch1")

# ....Spectral moments ----
# Note: based on model inspection, we will log-transform spectral moments
# Moments (M1-M4): COG, specsd, skew, kurt
# Include gender as covariate

# COG (M1)
# log-transformed
m_cog_habit <- lmerTest::lmer(
  log(cog) ~ group*cond_mask + 
    gender +
    (1 + (NMvM + SMvKN) | participant)+ #converged
    (1 | utterance),
  data = df_habitual
)
qqnorm(resid(m_cog_habit)) # approx normal
qqline(resid(m_cog_habit))
plot(resid(m_cog_habit)) # ok

sjPlot::tab_model(m_cog_habit)
# Effect of group (oc > pd)
# Effect of masks (nm > m; sm > kn)
# Effect of gender (f > m)
# No interactions

# Confirming direction of effects with plots of estimated main effects
sjPlot::plot_model(m_cog_habit, type="pred", terms = c("gender"))
sjPlot::plot_model(m_cog_habit, type="pred", terms = c("group"))
sjPlot::plot_model(m_cog_habit, type="pred", terms = c("cond_mask"))

# COG sd (M2)
# Log-transform
# var name = specsd
m_cogsd_habit <- lmerTest::lmer(
  log(specsd) ~ group*cond_mask +  
    gender +
    (1 + (NMvM + SMvKN) | participant)+ #converged
    (1 | utterance),
  data = df_habitual
)
# Check assumptions
qqnorm(resid(m_cogsd_habit)) # approx normal
qqline(resid(m_cogsd_habit))
plot(resid(m_cogsd_habit)) #ok

sjPlot::tab_model(m_cogsd_habit)
# No effect of group
# Effect of masks (nm > m; sm > kn)
# Effect of gender (f > m)
# no interactions



# Skew (M3)
# log-transform
m_skew_habit <- lmerTest::lmer(
  log(skew) ~ group*cond_mask +  
    gender +
    (1 + (NMvM + SMvKN) | participant)+ #converged
    (1 | utterance),
  control = lme4::lmerControl(optimizer ="Nelder_Mead"), # for convergence
  data = df_habitual
)
# Check assumptions
qqnorm(resid(m_skew_habit)) # approx normal
qqline(resid(m_skew_habit))
plot(resid(m_skew_habit)) # ok

sjPlot::tab_model(m_skew_habit)
# Effect of group (oc < pd)
# Effect of masks (nm < m; sm < kn)
# no effect of gender
# no interactions

# Kurtosis (M4)
# log-transform
m_kurt_habit <- lmerTest::lmer(
  log(kurt) ~ group*cond_mask +  
    gender +
    (1 + (NMvM + SMvKN) | participant)+ #converged
    (1 | utterance),
  data = df_habitual
)

# Check assumptions
qqnorm(resid(m_kurt_habit)) # approx normal
qqline(resid(m_kurt_habit))
plot(resid(m_kurt_habit)) # ok

sjPlot::tab_model(m_kurt_habit)
# Small effect of group (oc < pd)
# Effects of masks (nm < m; sm < kn)
# No effect of gender
# No interactions


#..RQ2: all speech conditions, no masks ----

df_nm <- df %>%
  filter(cond_mask=="nm",
         channel == "ch1")


# ....Spectral Moments ----
# All log-transformed
m_cog_nm <- lmerTest::lmer(
  log(cog) ~ group*cond_speech + gender +
    (1 + (HvCL + CvL) | participant)+ # converged
    (1 | utterance),
  data = df_nm
)
qqnorm(resid(m_cog_nm)) # normal
qqline(resid(m_cog_nm))
plot(resid(m_cog_nm)) # ok

sjPlot::tab_model(m_cog_nm)
# Effect of group (oc > pd)
# Effect of speech (habit < CL; clear < loud)
# Effect of gender (f > m)
# no interactions

m_cogsd_nm <- lmerTest::lmer(
  log(specsd) ~ group*cond_speech + gender +
    (1 + (HvCL + CvL) | participant)+ # converged
    (1 | utterance),
  data = df_nm
)

qqnorm(resid(m_cogsd_nm)) # approx normal
qqline(resid(m_cogsd_nm))
plot(resid(m_cogsd_nm)) # ok

sjPlot::tab_model(m_cogsd_nm)
# only gender is significant (f > m)

m_skew_nm <- lmerTest::lmer(
  log(skew) ~ group*cond_speech + gender +
    (1 + (HvCL + CvL) | participant)+ # converged
    (1 | utterance),
  data = df_nm
)

qqnorm(resid(m_skew_nm)) # approx normal
qqline(resid(m_skew_nm))
plot(resid(m_skew_nm)) # ok

sjPlot::tab_model(m_skew_nm)
# No effect of group
# Effect of speech (habit > CL; clear > loud)
# No effect of gender
# No interactions

m_kurt_nm <- lmerTest::lmer(
  log(kurt) ~ group*cond_speech + gender +
    (1 + (HvCL + CvL) | participant)+ # converged
    (1 | utterance),
  data = df_nm
)

qqnorm(resid(m_kurt_nm)) # approx normal
qqline(resid(m_kurt_nm))
plot(resid(m_kurt_nm)) # ok

sjPlot::tab_model(m_kurt_nm)
# No effect of group
# Effect of speech (habit > CL; clear > loud)
# No effect of gender
# No interactions

# Compare all
sjPlot::tab_model(m_cog_nm, m_cogsd_nm, m_skew_nm, m_kurt_nm)


# .. RQ3: all speech, masks ----
# all data (df)
# of interest: spectral tilt, mid 1-3khz, intensity (from knowles2022)

# ....Intensity ----
# Use corrected intensity (accounts for calibration factor)
m_int <- lmerTest::lmer(
  int_corrected ~ group*cond_speech*cond_mask + gender +
    (1 + (HvCL + CvL + NMvM + SMvKN) | participant)+ # converged
    (1 | utterance),
  data=subset(df,channel=="ch1")
)

# Check assumptions
qqnorm(resid(m_int)) # approx normal
qqline(resid(m_int))
plot(resid(m_int)) # ok

sjPlot::tab_model(m_int)
# All main effects significant EXCEPT group
# No effect of group
# Effect of speech (habit < CL; clear < loud)
# Effect of mask (nm > m; sm < kn; unexpected direction of surgical mask)
# Effect of gender (f < m)
# Several interactions (mostly small)

# Explore effects with forest plot
sjPlot::plot_model(m_int, show.values = TRUE, show.p = TRUE)

# ....Tilt ----
m_tilt <- lmerTest::lmer(
  tilt ~ group*cond_speech*cond_mask + gender +
    (1 + (HvCL + CvL + NMvM + SMvKN) | participant)+ # converged
    (1 | utterance),
  data=subset(df,channel=="ch1")
)

# Check assumptions
qqnorm(resid(m_tilt)) # approx normal
qqline(resid(m_tilt))
plot(resid(m_tilt)) # ok

sjPlot::tab_model(m_tilt)
# All main effects significant
# Effect of group (oc > pd)
# Effect of speech (habit < CL; clear < loud)
# Effect of mask (nm > m; sm > kn)
# Effect of gender (f > m)
# Very few interactions; all very small

# Explore 3-way interaction (not significant)
sjPlot::plot_model(m_tilt, type="pred",
                   terms = c("group",
                             "cond_mask",
                             "cond_speech"))

# ....Mid ----
m_mid <- lmerTest::lmer(
  mid_1_3k ~ group*cond_speech*cond_mask + gender +
    (1 + (HvCL + CvL + NMvM + SMvKN) | participant)+ # converged
    (1 | utterance),
  data=subset(df,channel=="ch1")
)

# Check assumptions
qqnorm(resid(m_mid)) # approx normal
qqline(resid(m_mid))
plot(resid(m_mid)) # ok

sjPlot::tab_model(m_mid)
# All main effects significant EXCEPT gender
# Effect of group (oc > pd)
# Effect of speech (habit < CL; clear < loud)
# Effect of mask (nm > m; sm > kn)
# Effect of gender (f > m)
# Some interactions; all very small

# Explore 3-way interaction (very small effect size for one comparison)
sjPlot::plot_model(m_mid, type="pred",
                   terms = c("group",
                             "cond_mask",
                             "cond_speech"))


# Coefs ----
# ..RQ1 ----
coefs_cog_h <- make_coefs_with_effsize(m_cog_habit, "M1-COG")
coefs_cogsd_h <- make_coefs_with_effsize(m_cogsd_habit, "M2-COGsd")
coefs_skew_h <- make_coefs_with_effsize(m_skew_habit, "M3-skew")
coefs_kurt_h <- make_coefs_with_effsize(m_kurt_habit, "M4-kurt")

coefs_h_all <- rbind(coefs_cog_h, coefs_cogsd_h,
                     coefs_skew_h, coefs_kurt_h)

coefs_h_all <- coefs_h_all %>%
  select(Contrast, measure, everything()) %>%
  rename("Measure" = "measure",
         "t" = "t value",
         "p" = "Pr(>|t|)") %>%
  mutate("Effect size" = case_when(
    as.numeric(p) > 0.05 ~ "ns",
    d < 0.2 ~ "very small",
    d >= 0.2 & d <0.5 ~ "small",
    d >= 0.5 & d <0.8 ~ "medium",
    d >= 0.8 ~ "large"
  )) %>%
  rename("Effect size parameter" = "d") %>%
  mutate(Contrast = str_remove_all(Contrast,c("cond_speech"))) %>%
  mutate(Contrast = str_remove_all(Contrast,c("cond_mask"))) %>%
  mutate(Contrast = str_remove_all(Contrast,c("group"))) %>%
  mutate(Contrast = str_replace_all(Contrast,"genderW vs M","Women vs. Men")) %>%
  mutate(p = ifelse(p < 0.001,"<0.001",p))


# beta, p dfs
habit_cogB <- make_b(m_cog_habit)
habit_cogsdB <- make_b(m_cogsd_habit)
habit_skewB <- make_b(m_skew_habit)
habit_kurtB <- make_b(m_kurt_habit)

habit_cogP <- make_p(m_cog_habit)
habit_cogsdP <- make_p(m_cogsd_habit)
habit_skewP <- make_p(m_skew_habit)
habit_kurtP <- make_p(m_kurt_habit)


# ..RQ2 ----
coefs_cog_nm <- make_coefs_with_effsize(m_cog_nm, "M1-COG")
coefs_cogsd_nm <- make_coefs_with_effsize(m_cogsd_nm, "M2-COGsd")
coefs_skew_nm <- make_coefs_with_effsize(m_skew_nm, "M3-skew")
coefs_kurt_nm <- make_coefs_with_effsize(m_kurt_nm, "M4-kurt")

coefs_nm_all <- rbind(coefs_cog_nm, coefs_cogsd_nm,
                      coefs_skew_nm, coefs_kurt_nm)

coefs_nm_all <- coefs_nm_all %>%
  select(Contrast, measure, everything()) %>%
  rename("Measure" = "measure",
         "t" = "t value",
         "p" = "Pr(>|t|)") %>%
  mutate("Effect size" = case_when(
    as.numeric(p) > 0.05 ~ "ns",
    d < 0.2 ~ "very small",
    d >= 0.2 & d <0.5 ~ "small",
    d >= 0.5 & d <0.8 ~ "medium",
    d >= 0.8 ~ "large"
  )) %>%
  rename("Effect size parameter" = "d") %>%
  mutate(Contrast = str_remove_all(Contrast,c("cond_speech"))) %>%
  mutate(Contrast = str_remove_all(Contrast,c("cond_mask"))) %>%
  mutate(Contrast = str_remove_all(Contrast,c("group"))) %>%
  mutate(Contrast = str_replace_all(Contrast,"genderW vs M","Women vs. Men")) %>%
  mutate(p = ifelse(p < 0.001,"<0.001",p))

# Make b, p dfs
nm_cogB <- make_b(m_cog_nm)
nm_cogsdB <- make_b(m_cogsd_nm)
nm_skewB <- make_b(m_skew_nm)
nm_kurtB <- make_b(m_kurt_nm)

nm_cogP <- make_p(m_cog_nm)
nm_cogsdP <- make_p(m_cogsd_nm)
nm_skewP <- make_p(m_skew_nm)
nm_kurtP <- make_p(m_kurt_nm)

# ..RQ3 ----
coefs_int <- make_coefs_with_effsize(m_int, "intensity")
coefs_tilt <- make_coefs_with_effsize(m_tilt, "tilt")
coefs_mid <- make_coefs_with_effsize(m_mid, "mid 1-3khz")

coefs_all <- rbind(coefs_int, coefs_tilt, coefs_mid)

coefs_all <- coefs_all %>%
  select(Contrast, measure, everything()) %>%
  rename("Measure" = "measure",
         "t" = "t value",
         "p" = "Pr(>|t|)") %>%
  mutate("Effect size" = case_when(
    as.numeric(p) > 0.05 ~ "ns",
    d < 0.2 ~ "very small",
    d >= 0.2 & d <0.5 ~ "small",
    d >= 0.5 & d <0.8 ~ "medium",
    d >= 0.8 ~ "large"
  )) %>%
  rename("Effect size parameter" = "d") %>%
  mutate(Contrast = str_remove_all(Contrast,c("cond_speech"))) %>%
  mutate(Contrast = str_remove_all(Contrast,c("cond_mask"))) %>%
  mutate(Contrast = str_remove_all(Contrast,c("group"))) %>%
  mutate(Contrast = str_replace_all(Contrast,"genderW vs M","Women vs. Men")) %>%
  mutate(p = ifelse(p < 0.001,"<0.001",p))

# Make b, p dfs
intB <- make_b(m_int)
tiltB <- make_b(m_tilt)
midB <- make_b(m_mid)

intP <- make_p(m_int)
tiltP <- make_p(m_tilt)
midP <- make_p(m_mid)

# Ready to be sourced!

