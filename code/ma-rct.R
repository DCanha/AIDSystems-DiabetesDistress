# Load necessary libraries
library(readxl)  # For read_excel
library(dplyr)  # For data manipulation
library(esc)    # For effect size calculation
library(metafor)  # For meta-analysis
library(meta)   # For meta-analysis and forest plots

# Set a seed for reproducibility
set.seed(5555)

# All data values are removed directly from the articles (from the search results)

############################################################
#### Adult population ####
############################################################

# Convert median to mean for McAuley 2020
set.seed(5555)
meansd_exp_t0 <- bc.mean.sd(q1.val = 12.5, med.val = 20, q3.val = 38.8, n = 61)
set.seed(5555)
meansd_control_t0 <- bc.mean.sd(q1.val = 11.9, med.val = 17.5, q3.val = 34.4, n= 59)
set.seed(5555)
meansd_exp_t1 <- bc.mean.sd(q1.val = 10.2, med.val = 16.7, q3.val = 27.4, n = 61)
set.seed(5555)
meansd_control_t1 <- bc.mean.sd(q1.val = 9.5, med.val = 21.2, q3.val = 36.2, n= 59)

# Prepare data adults with p-values
# Considering Kudva 2021 at 3 months for consistency

data_adults <- data.frame(
  study = c("Kudva 2021", "Mc Auley 2020", "Mc Auley 2022"),
  type_rct = c("parallel", "parallel", "cross-over"),
  pi = c(0.04, 0.10, 0.46),
  n_exp = c(112, 61, 30),
  # for crossover design, we put n_control = 0 for calculations
  n_control = c(56, 59, 0),
  m_exp_0 = c(1.8, meansd_exp_t0$est.mean, 4.30),
  m_control_0 = c(1.8, meansd_control_t0$est.mean, 4.60),
  sd_exp_0 = c(0.6, meansd_exp_t0$est.sd, 2.90),
  sd_control_0 = c(0.7, meansd_control_t0$est.sd, 3.20),
  # for crossover design, we put NA for calculations (there is only two values; not 4)
  m_exp_1 = c(1.7, meansd_exp_t1$est.mean, NA),
  m_control_1 = c(1.9, meansd_control_t1$est.mean, NA),
  sd_exp_1 = c(0.60, meansd_exp_t1$est.sd, NA),
  sd_control_1 = c(0.70, meansd_control_t1$est.sd, NA),
  md_between = c(-0.17, NA, -0.3),
  #in months
  study_duration = c(3, 6, 4),
  risk_bias = c("low", "low", "some")
)

# Calculate correlation coefficients
data_adults <- escalc(measure="UCOR", pi = pi, ni = n_exp+n_control, data = data_adults)
data_adults <- data_adults[,-c(3,18)]  # Remove unnecessary columns
colnames(data_adults)[16] <- "ri"  # Rename column for clarity

# Calculate raw mean differences
data_adults$md_between <- ifelse(is.na(data_adults$md_between), 
                                 (data_adults$m_exp_1 - data_adults$m_exp_0) - 
                                   (data_adults$m_control_1 - data_adults$m_control_0),
                                 data_adults$md_between)

# Calculate pooled SDs

data_adults$sd_exp_pooled <- ifelse(data_adults$type_rct == "parallel",
                                    sqrt(data_adults$sd_exp_0^2 + data_adults$sd_exp_1^2 - 
                                           2 * data_adults$ri * data_adults$sd_exp_0 * data_adults$sd_exp_1),
                                    #for crossover is already sd exp and sd control
                                    data_adults$sd_exp_0)

data_adults$sd_control_pooled <- ifelse(data_adults$type_rct == "parallel",
                                    sqrt(data_adults$sd_control_0^2 + data_adults$sd_control_1^2 - 
                                           2 * data_adults$ri * data_adults$sd_control_0 * data_adults$sd_control_1),
                                    #for crossover is already sd exp and sd control
                                    data_adults$sd_control_0)

data_adults$sd_change_pooled <- ifelse(data_adults$type_rct == "parallel",
                                       sqrt(((data_adults$n_exp - 1) * data_adults$sd_exp_pooled^2 + 
                                               (data_adults$n_control - 1) * data_adults$sd_control_pooled^2) / 
                                              (data_adults$n_exp + data_adults$n_control - 2)),
                                       #for crossover is different
                                       sqrt((data_adults$sd_exp_pooled^2 + 
                                               data_adults$sd_control_pooled^2) / 2))

# Calculate SMDs
data_adults$smd_between <- data_adults$md_between / data_adults$sd_change_pooled

# Calculate SEs
data_adults$se_between <- ifelse(data_adults$type_rct == "parallel",
                                 sqrt((data_adults$n_exp + data_adults$n_control) / (data_adults$n_exp * data_adults$n_control) + 
                                        data_adults$smd_between^2 / (2 * (data_adults$n_exp + data_adults$n_control))),
                                 #for crossover is different
                                 sqrt(1 / data_adults$n_exp + data_adults$smd_between^2 / 
                                        (2 * data_adults$n_exp)) * sqrt(2 * (1 - data_adults$ri)))

# Meta-analysis using a random-effects model
random_adults_RCT <- metagen(TE = smd_between, seTE = se_between, studlab = study, data = data_adults, 
                         sm = "SMD", fixed = FALSE, random = TRUE, 
                         method.tau = "REML", hakn = TRUE, 
                         title = "Closed-loop RCT in adults")
summary(random_adults_RCT)

# Forest plot visualization
forest.meta(random_adults_RCT, sortvar = risk_bias, prediction = TRUE, print.tau2 = FALSE, 
            leftcols = c("studlab"), leftlabs = c("Author"))

# Funnel plot visualization
funnel.meta(random_adults_RCT, studlab = TRUE, pos.studlab = 3, col = "blue")


############################################################
#### Pediatric population ####
############################################################

# Prepare data children with p-values
# Considering Hood 2022 at 3 months for consistency

data_peds <- data.frame(
  study = c("Abraham 2021", "Hood 2022", "Cobry 2021"),
  type_rct = c("parallel", "parallel", "parallel"),
  pi = c(0.14, 0.99, 0.54),
  n_exp = c(58, 37, 78),
  n_control = c(53, 34, 22),
  m_exp_0 = c(30.3, 51.93, 35),
  m_control_0 = c(28.9, 47.65, 37),
  sd_exp_0 = c(19.7, 19.92, 16),
  sd_control_0 = c(18.9, 23.61, 14),
  m_exp_1 = c(27.9, 52.71, 34),
  m_control_1 = c(31.1, 48.56, 37),
  sd_exp_1 = c(21.8, 21.17, 16),
  sd_control_1 = c(20.7, 20.09, 17),
  md_between = c(-4.5, NA, -1.5),
  #in months
  study_duration = c(6, 3, 3.75),
  risk_bias = c("low", "some", "low")
)

# Calculate correlation coefficients
data_peds <- escalc(measure="UCOR", pi = pi, ni = n_exp+n_control, data = data_peds)
data_peds <- data_peds[,-c(3,18)]  # Remove unnecessary columns
colnames(data_peds)[16] <- "ri"  # Rename column for clarity

# Calculate raw mean differences
data_peds$md_between <- ifelse(is.na(data_peds$md_between), 
                                 (data_peds$m_exp_1 - data_peds$m_exp_0) - 
                                   (data_peds$m_control_1 - data_peds$m_control_0),
                                 data_peds$md_between)

# Calculate pooled SDs

data_peds$sd_exp_pooled <- ifelse(data_peds$type_rct == "parallel",
                                    sqrt(data_peds$sd_exp_0^2 + data_peds$sd_exp_1^2 - 
                                           2 * data_peds$ri * data_peds$sd_exp_0 * data_peds$sd_exp_1),
                                    #for crossover is already sd exp and sd control
                                    data_peds$sd_exp_0)

data_peds$sd_control_pooled <- ifelse(data_peds$type_rct == "parallel",
                                        sqrt(data_peds$sd_control_0^2 + data_peds$sd_control_1^2 - 
                                               2 * data_peds$ri * data_peds$sd_control_0 * data_peds$sd_control_1),
                                        #for crossover is already sd exp and sd control
                                        data_peds$sd_control_0)

data_peds$sd_change_pooled <- ifelse(data_peds$type_rct == "parallel",
                                       sqrt(((data_peds$n_exp - 1) * data_peds$sd_exp_pooled^2 + 
                                               (data_peds$n_control - 1) * data_peds$sd_control_pooled^2) / 
                                              (data_peds$n_exp + data_peds$n_control - 2)),
                                       #for crossover is different
                                       sqrt((data_peds$sd_exp_pooled^2 + 
                                               data_peds$sd_control_pooled^2) / 2))

# Calculate SMDs
data_peds$smd_between <- data_peds$md_between / data_peds$sd_change_pooled

# Calculate SEs
data_peds$se_between <- ifelse(data_peds$type_rct == "parallel",
                                 sqrt((data_peds$n_exp + data_peds$n_control) / (data_peds$n_exp * data_peds$n_control) + 
                                        data_peds$smd_between^2 / (2 * (data_peds$n_exp + data_peds$n_control))),
                                 #for crossover is different
                                 sqrt(1 / data_peds$n_exp + data_peds$smd_between^2 / 
                                        (2 * data_peds$n_exp)) * sqrt(2 * (1 - data_peds$ri)))

# Meta-analysis using a random-effects model
random_peds_RCT <- metagen(TE = smd_between, seTE = se_between, studlab = study, data = data_peds, 
                             sm = "SMD", fixed = FALSE, random = TRUE, 
                             method.tau = "REML", hakn = TRUE, 
                             title = "Closed-loop RCT in adults")
summary(random_peds_RCT)

# Forest plot visualization
forest.meta(random_peds_RCT, sortvar = risk_bias, prediction = TRUE, print.tau2 = FALSE, 
            leftcols = c("studlab"), leftlabs = c("Author"))

# Funnel plot visualization
funnel.meta(random_peds_RCT, studlab = TRUE, pos.studlab = 3, col = "blue")


############################################################
#### Caregiver population ####
############################################################

# compute mean and sd for hood 2022 

#sample sizes are the same in the 3 subscales so:
mean_exp_t0 <- c(1.06, 1.01, 1.48)
mean_control_t0 <- c(1.01, 0.88, 1.64)
mean_exp_t1 <- c(0.92, 1.08, 1.36)
mean_control_t1 <- c(1.13, 1.01, 1.80)
mean_exp_t2 <- c(0.97, 1.22, 1.53)
mean_control_t2 <- c(1.22, 1.09, 1.79)

sd_exp_t0 <- c(0.94, 0.87, 0.87)
sd_control_t0 <- c(0.90, 0.76, 1.11)
sd_exp_t1 <- c(0.85, 0.81, 0.81)
sd_control_t1 <- c(1.01, 0.90, 1.00)
sd_exp_t2 <- c(0.79, 0.96, 0.86)
sd_control_t2 <- c(1.21, 1.02, 1.18)

# Compute pooled mean
pmean_exp_t0 <- mean(mean_exp_t0)
pmean_control_t0 <- mean(mean_control_t0)
pmean_exp_t1 <- mean(mean_exp_t1)
pmean_control_t1 <- mean(mean_control_t1)
pmean_exp_t2 <- mean(mean_exp_t2)
pmean_control_t2 <- mean(mean_control_t2)

# Compute pooled variance
pvar_exp_t0 <- mean(sd_exp_t0^2)
pvar_control_t0 <- mean(sd_control_t0^2)
pvar_exp_t1 <- mean(sd_exp_t1^2)
pvar_control_t1 <- mean(sd_control_t1^2)
pvar_exp_t2 <- mean(sd_exp_t2^2)
pvar_control_t2 <- mean(sd_control_t2^2)

# Compute pooled standard deviation
psd_exp_t0 <- sqrt(pvar_exp_t0)
psd_control_t0 <- sqrt(pvar_control_t0)
psd_exp_t1 <- sqrt(pvar_exp_t1)
psd_control_t1 <- sqrt(pvar_control_t1)
psd_exp_t2 <- sqrt(pvar_exp_t2)
psd_control_t2 <- sqrt(pvar_control_t2)

# Prepare data caregivers with p-values
# Considering Hood 2022 at 3 months (t1) for consistency

data_care <- data.frame(
  study = c("Kudva 2021", "Hood 2022", "Cobry 2021"),
  type_rct = c("parallel", "parallel", "parallel"),
  pi = c(0.29, 0.39, 0.08),
  n_exp = c(31, 48, 78),
  n_control = c(17, 46, 22),
  m_exp_0 = c(1.2, pmean_exp_t0, 43),
  m_control_0 = c(1, pmean_control_t0, 46),
  sd_exp_0 = c(0.7, psd_exp_t0, 16),
  sd_control_0 = c(0.5, psd_control_t0, 15),
  m_exp_1 = c(1, pmean_exp_t1, 36),
  m_control_1 = c(1.2, pmean_control_t1, 45),
  sd_exp_1 = c(0.8, psd_exp_t1, 12),
  sd_control_1 = c(0.7, psd_control_t1, 17),
  md_between = c(-0.34, NA, -7.4),
  #in months
  study_duration = c(6, 3, 3.75),
  risk_bias = c("low", "some", "low")
)

# Calculate correlation coefficients
data_care <- escalc(measure="UCOR", pi = pi, ni = n_exp+n_control, data = data_care)
data_care <- data_care[,-c(3,18)]  # Remove unnecessary columns
colnames(data_care)[16] <- "ri"  # Rename column for clarity

# Calculate raw mean differences
data_care$md_between <- ifelse(is.na(data_care$md_between), 
                               (data_care$m_exp_1 - data_care$m_exp_0) - 
                                 (data_care$m_control_1 - data_care$m_control_0),
                               data_care$md_between)

# Calculate pooled SDs

data_care$sd_exp_pooled <- ifelse(data_care$type_rct == "parallel",
                                  sqrt(data_care$sd_exp_0^2 + data_care$sd_exp_1^2 - 
                                         2 * data_care$ri * data_care$sd_exp_0 * data_care$sd_exp_1),
                                  #for crossover is already sd exp and sd control
                                  data_care$sd_exp_0)

data_care$sd_control_pooled <- ifelse(data_care$type_rct == "parallel",
                                      sqrt(data_care$sd_control_0^2 + data_care$sd_control_1^2 - 
                                             2 * data_care$ri * data_care$sd_control_0 * data_care$sd_control_1),
                                      #for crossover is already sd exp and sd control
                                      data_care$sd_control_0)

data_care$sd_change_pooled <- ifelse(data_care$type_rct == "parallel",
                                     sqrt(((data_care$n_exp - 1) * data_care$sd_exp_pooled^2 + 
                                             (data_care$n_control - 1) * data_care$sd_control_pooled^2) / 
                                            (data_care$n_exp + data_care$n_control - 2)),
                                     #for crossover is different
                                     sqrt((data_care$sd_exp_pooled^2 + 
                                             data_care$sd_control_pooled^2) / 2))

# Calculate SMDs
data_care$smd_between <- data_care$md_between / data_care$sd_change_pooled

# Calculate SEs
data_care$se_between <- ifelse(data_care$type_rct == "parallel",
                               sqrt((data_care$n_exp + data_care$n_control) / (data_care$n_exp * data_care$n_control) + 
                                      data_care$smd_between^2 / (2 * (data_care$n_exp + data_care$n_control))),
                               #for crossover is different
                               sqrt(1 / data_care$n_exp + data_care$smd_between^2 / 
                                      (2 * data_care$n_exp)) * sqrt(2 * (1 - data_care$ri)))

# Meta-analysis using a random-effects model
random_care_RCT <- metagen(TE = smd_between, seTE = se_between, studlab = study, data = data_care, 
                           sm = "SMD", fixed = FALSE, random = TRUE, 
                           method.tau = "REML", hakn = TRUE, 
                           title = "Closed-loop RCT in adults")
summary(random_care_RCT)

# Forest plot visualization
forest.meta(random_care_RCT, sortvar = risk_bias, prediction = TRUE, print.tau2 = FALSE, 
            leftcols = c("studlab"), leftlabs = c("Author"))

# Funnel plot visualization
funnel.meta(random_care_RCT, studlab = TRUE, pos.studlab = 3, col = "blue")
