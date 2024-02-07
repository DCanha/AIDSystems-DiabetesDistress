# Load necessary libraries
library(metafor)
library(dplyr)
library(meta)
library(esc)
library(gsl)
library(estmeansd)

# Set a seed for reproducibility
set.seed(5555)

# All data values are removed directly from the articles (from the search results)

############################################################
#### Adult population ####
############################################################

# Adult data preparation

# We have t-test value for one study; we will use it to calculate r
dataT <- data.frame(
  study = "Beato-Víbora 2020",
  ti = 3.258,
  ni = 51,
  m_pre = 40,
  m_post = 34,
  sd_pre = 19,
  sd_post = 16,
  md_within = -5.1,
  study_duration = 3,
  quality = 6
)

# Calculate correlation coefficients for t-test data
dataTCorr <- escalc(measure="UCOR", ti = ti, ni = ni, data = dataT)
dataTCorr <- dataTCorr[,-c(2,12)]  # Remove unnecessary columns
colnames(dataTCorr)[10] <- "ri"  # Rename column for clarity

# Prepare data with p-values 
dataP <- data.frame(
  study = c("Schneider-Utaka 2023", "Polonsky 2022", "Bisio 2022", "Reznik 2023", "Akiyama 2023"),
  pi = c(0.03, 0.0001, 0.046, 0.001, 0.43),
  ni = c(35, 115, 15, 202, 22),
  m_pre = c(1.54, 1.64, 1.5, 47.3, 40.8),
  m_post = c(1.43, 1.48, 1.35, 39.9, 38.6),
  sd_pre = c(0.35, 0.51, 0.36, 18.9, 21.3),
  sd_post = c(0.32, 0.4, 0.25, 20.8, 23.4),
  md_within = c(NA, -0.16, NA, -6.9, NA),
  study_duration = c(3.7, 3, 1, 6, 3),
  quality = c(3, 2, 5, 2, 4)
)

dataPCorr <- escalc(measure="UCOR", pi = pi, ni = ni, data = dataP)
dataPCorr <- dataPCorr[,-c(2,12)]  # Remove unnecessary columns
colnames(dataPCorr)[10] <- "ri"  # Rename column for clarity

# Combine adult data
dataAdults <- rbind(dataTCorr, dataPCorr)

# Calculate raw mean differences and SMD
dataAdults$md_within <- ifelse(is.na(dataAdults$md_within), dataAdults$m_post - dataAdults$m_pre,
                               dataAdults$md_within)
dataAdults$smd_within <- dataAdults$md_within / dataAdults$sd_pre
dataAdults$se_within <- sqrt(((2 * (1 - dataAdults$ri)) / dataAdults$ni) + 
                               (dataAdults$smd_within^2 / (2 * dataAdults$ni)))

# Meta-analysis using a random-effects model
random_adults <- metagen(TE = smd_within, seTE = se_within, studlab = study, data = dataAdults, 
                        sm = "SMD", fixed = FALSE, random = TRUE, 
                        method.tau = "REML", hakn = TRUE, 
                        title = "Closed-loop before-after in adults")
summary(random_adults)

# Forest plot visualization
forest.meta(random_adults, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, 
            leftcols = c("studlab"), leftlabs = c("Author"))


# Funnel plot visualization
funnel.meta(random_adults, studlab = TRUE, pos.studlab = 3, col = "blue")

# Categorize studies by quality
dataAdults$Group <- ifelse(dataAdults$quality < 6, "Good/Fair", "Poor")

# Perform meta-analysis considering subgroup analysis for study quality
adults_quality <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                          data = dataAdults[order(dataAdults$Group == "Good/Fair", decreasing = TRUE), ], 
                          sm = "SMD", fixed = FALSE, random = TRUE, method.tau = "REML", hakn = TRUE, 
                          title = "Closed-loop before-after in adults", 
                          subgroup = Group, test.subgroup = FALSE)

# Summary of the meta-analysis with subgroup analysis
summary(adults_quality)

# Forest plot visualizing the subgroup analysis
forest.meta(adults_quality, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)


############################################################
#### Pediatric population ####
############################################################

# Convert median to mean for Cobry 2022
set.seed(5555)
meansd_pre <- bc.mean.sd(q1.val = 17, med.val = 23, q3.val = 32, n = 49)
set.seed(5555)
meansd_post <- bc.mean.sd(q1.val = 17, med.val = 21, q3.val = 31, n = 49)

# Prepare child data

# We have t-test value for one study; we will use it to calculate r
dataChildt <- data.frame(
  study = "Cobry 2020",
  ti = -0.23,
  ni = 23,
  m_pre = 22.39,
  m_post = 23.04,
  sd_pre = 15.16,
  sd_post = 17.57,
  md_within = NA,
  #in months
  study_duration = 3,
  quality = 8
)

# Calculate correlation coefficients for t-test data
dataChildT <- escalc(measure="UCOR", ti = ti, ni=ni, data=dataChildt)
dataChildT <- dataChildT[,-c(2,12)]  # Remove unnecessary columns
colnames(dataChildT)[10] <- "ri"  # Rename column for clarity

dataChildp <- data.frame(
  study = c("Bisio 2021", "Cobry 2022", "Gianini 2022", "Reznik 2023", "Berget 2020", "Cobry 2021"),
  pi = c(NA, 0.153, 0.001, NA, NA, 0.1),
  ni = c(13, 49, 24, 55, 92, 22),
  m_pre = c(17.26, meansd_pre$est.mean, 19.32, 34.3, 35.8, 37),
  m_post = c(18.15, meansd_post$est.mean, 8.65, 34.6, 36.4, 37),
  sd_pre = c(13.17, meansd_pre$est.sd, 12.25, 19.4, 2.4 * sqrt(92), 17),
  sd_post = c(16.22, meansd_post$est.sd, 8.26, 23.5, 3 * sqrt(92), 19),
  md_within = c(NA, NA, NA, 1.4, NA, NA),
  study_duration = c(1, 3.7, 4, 6, 3, 2.8),
  quality = c(7, 2, 6, 2, 0, 5)
)

# Calculate correlation coefficients
dataChildP <- escalc(measure="UCOR", pi = pi, ni = ni, data = dataChildp)
dataChildP <- dataChildP[,-c(2,12)]  # Remove unnecessary columns
colnames(dataChildP)[10] <- "ri"  # Rename column for clarity

# Combine child data
dataChild <- rbind(dataChildT, dataChildP)

# Adjust correlation coefficients based on study similarities
#Bisio is similar to Cobry 2020 in terms of ni and m_pre m_post
dataChild$ri[2] <- dataChild$ri[1]
#Reznik is similar to Cobry 2022 in the same vars; but with opposite signs
dataChild$ri[5] <- -dataChild$ri[3]
#Berget is similar to Reznik
dataChild$ri[6] <- dataChild$ri[5]

# Calculate raw mean differences and standardized mean differences
dataChild$md_within <- ifelse(is.na(dataChild$md_within), dataChild$m_post - dataChild$m_pre,
                              dataChild$md_within)
dataChild$smd_within <- dataChild$md_within / dataChild$sd_pre
dataChild$se_within <- sqrt(((2 * (1 - dataChild$ri)) / dataChild$ni) + 
                              (dataChild$smd_within^2 / (2 * dataChild$ni)))

# Perform meta-analysis
random_child <- metagen(TE = smd_within, seTE = se_within, studlab = study, data = dataChild, 
                        sm = "SMD", fixed = FALSE, random = TRUE, 
                        method.tau = "REML", hakn = TRUE, 
                        title = "Closed-loop before-after in children")
summary(random_child)

# Visualize with a forest plot
forest.meta(random_child, sortvar = dataChild$quality, prediction = TRUE, print.tau2 = FALSE, 
            leftcols = c("studlab"), leftlabs = c("Author"))


# Funnel plot visualization
funnel.meta(random_child, studlab = TRUE, pos.studlab = 3, col = "blue")

# High heterogeneity 
# Perform sensitivity analysis to evaluate the impact of each study

sensitivity_results <- lapply(1:length(dataChild$smd_within), function(i) {
  excluded_study <- dataChild$study[i]
  ma_temp <- metagen(TE = dataChild$smd_within[-i], seTE = dataChild$se_within[-i], 
                     studlab = dataChild$study[-i], method.tau = "REML", hakn = TRUE)
  list(excluded = excluded_study, ma_temp = ma_temp)
})

# Format the results into a readable table
sensitivity_table <- data.frame(
  Excluded_Study = sapply(sensitivity_results, function(x) x$excluded),
  SMD_95CI = sapply(sensitivity_results, function(x) paste(round(x$ma_temp$TE.random, 2), "[", round(x$ma_temp$lower.random, 2), "-", round(x$ma_temp$upper.random, 2), "]")),
  I2_pvalue = sapply(sensitivity_results, function(x) {
    pval <- round(x$ma_temp$pval.Q, 2)
    p_text <- if (pval < 0.01) "p<0.01" else paste("p=", format.pval(pval, digits = 2))
    paste(round(x$ma_temp$I2 * 100, 0), "%", " (", p_text, ")")
  })
)

# View the sensitivity analysis table
print(sensitivity_table)

# Extract required data
TE = sapply(sensitivity_results, function(x) x$ma_temp$TE.random)
lower = sapply(sensitivity_results, function(x) x$ma_temp$lower.random)
upper = sapply(sensitivity_results, function(x) x$ma_temp$upper.random)
labels = sapply(sensitivity_results, function(x) paste(x$excluded))

# Create a new meta-analysis object without  performing a new analysis - just for visualization
sensitivity_child <- metagen(
  TE = c(random_child$TE.random, TE),
  lower = c(random_child$lower.random, lower),
  upper = c(random_child$upper.random, upper),
  studlab = c("None", labels),
  common = FALSE,
  random = FALSE,
  sm = "Random effects model",
  method.tau = "ML"
)

# Plot
forest(sensitivity_child,
       weight.study = "same",
       xlim = c(min(lower) - 0.1, max(upper) + 0.3), # Adjust as needed
       main = "Sensitivity Analysis: Leave-One-Out",
       rightlabs = c("SMD", "95% CI"),
       leftcols = c("studlab"),
       leftlabs = c("Excluded Study")
)

# Categorize studies by quality
dataChild$Group <- ifelse(dataChild$quality < 6, "Good/Fair", "Poor")

# Perform meta-analysis considering subgroup analysis for study quality
child_quality <- metagen(TE = smd_within, seTE = se_within, studlab = study, data = dataChild[order(dataChild$Group == "Good/Fair", decreasing = TRUE), ], sm = "SMD", fixed = FALSE, random = TRUE, method.tau = "REML", hakn = TRUE, title = "Closed-loop before-after in children", subgroup = Group, test.subgroup = FALSE)

# Summary of the meta-analysis with subgroup analysis
summary(child_quality)

# Forest plot visualizing the subgroup analysis
forest.meta(child_quality, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)

############################################################
#### Caregiver population ####
############################################################

# Convert median to mean for Cobry 2022
set.seed(5555)
meansd_pre <- bc.mean.sd(q1.val = 36, med.val = 44, q3.val = 56,n = 49)
set.seed(5555)
meansd_post <- bc.mean.sd(q1.val = 26, med.val = 35, q3.val = 43, n= 49)

# Prepare caregiver data
# Removing Cobry 2020 because the values do not make sense
dataCar <- data.frame(
  study = c("Bisio 2021", "Cobry 2022", "Berget 2020", "Cobry 2021"),
  pi = c(0.032, 0.001, -0.52, NA),
  ni = c(13, 49, 89, 22),
  m_pre = c(51.03, meansd_pre$est.mean, 44.2, 45),
  m_post = c(36.85, meansd_post$est.mean, 46.2, 40),
  sd_pre = c(17.19, meansd_pre$est.sd, 2 * sqrt(89), 17),
  sd_post = c(19.11, meansd_post$est.sd, 2.3 * sqrt(89), 17),
  md_within = rep(NA, 4),
  study_duration = c(1, 3.7, 3, 2.8),
  quality = c(7, 2, 0, 5)
)

# Calculate correlation coefficients
dataCare <- escalc(measure="COR", pi = pi, ni = ni, data = dataCar)
dataCare <- dataCare[,-c(2,12)]  # Remove unnecessary columns
colnames(dataCare)[10] <- "ri"  # Rename column for clarity

# Adjust correlation coefficients for Cobry 2021
dataCare$ri[4] <- mean(dataCare$ri, na.rm = TRUE)

# Calculate raw mean differences and standardized mean differences
dataCare$md_within <- ifelse(is.na(dataCare$md_within), dataCare$m_post - dataCare$m_pre, 
                             dataCare$md_within)
dataCare$smd_within <- dataCare$md_within / dataCare$sd_pre
dataCare$se_within <- sqrt(((2 * (1 - dataCare$ri)) / dataCare$ni) + 
                             (dataCare$smd_within^2 / (2 * dataCare$ni)))

# Perform meta-analysis
random_care <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                       data = dataCare, sm = "SMD", fixed = FALSE, random = TRUE, 
                       method.tau = "REML", hakn = TRUE, 
                       title = "Closed-loop before-after in caregivers")
summary(random_care)

# Visualize with a forest plot
forest.meta(random_care, sortvar = dataCare$quality, prediction = TRUE, 
            print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"))


# Funnel plot visualization
funnel.meta(random_care,  studlab = TRUE, pos.studlab = 3, col = "blue")

# High heterogeneity 
# Perform sensitivity analysis to evaluate the impact of each study

sensitivity_results <- lapply(1:length(dataCare$smd_within), function(i) {
  excluded_study <- dataCare$study[i]
  ma_temp <- metagen(TE = dataCare$smd_within[-i], seTE = dataCare$se_within[-i], 
                     studlab = dataCare$study[-i], method.tau = "REML", hakn = TRUE)
  list(excluded = excluded_study, ma_temp = ma_temp)
})

# Format the results into a readable table
sensitivity_table <- data.frame(
  Excluded_Study = sapply(sensitivity_results, function(x) x$excluded),
  SMD_95CI = sapply(sensitivity_results, function(x) paste(round(x$ma_temp$TE.random, 2), "[", round(x$ma_temp$lower.random, 2), "-", round(x$ma_temp$upper.random, 2), "]")),
  I2_pvalue = sapply(sensitivity_results, function(x) {
    pval <- round(x$ma_temp$pval.Q, 2)
    p_text <- if (pval < 0.01) "p<0.01" else paste("p=", format.pval(pval, digits = 2))
    paste(round(x$ma_temp$I2 * 100, 0), "%", " (", p_text, ")")
  })
)

# View the sensitivity analysis table
print(sensitivity_table)

# Extract required data
TE = sapply(sensitivity_results, function(x) x$ma_temp$TE.random)
lower = sapply(sensitivity_results, function(x) x$ma_temp$lower.random)
upper = sapply(sensitivity_results, function(x) x$ma_temp$upper.random)
labels = sapply(sensitivity_results, function(x) paste(x$excluded))

# Create a new meta-analysis object without  performing a new analysis - just for visualization
sensitivity_care <- metagen(
  TE = c(random_care$TE.random, TE),
  lower = c(random_care$lower.random, lower),
  upper = c(random_care$upper.random, upper),
  studlab = c("None", labels),
  common = FALSE,
  random = FALSE,
  sm = "Random effects model",
  method.tau = "ML"
)

# Plot
forest(sensitivity_care,
       weight.study = "same",
       xlim = c(min(lower) - 0.1, max(upper) + 0.3), # Adjust as needed
       main = "Sensitivity Analysis: Leave-One-Out",
       rightlabs = c("SMD", "95% CI"),
       leftcols = c("studlab"),
       leftlabs = c("Excluded Study")
)

# Categorize studies by quality
dataCare$Group <- ifelse(dataCare$quality < 6, "Good/Fair", "Poor")

# Perform meta-analysis considering subgroup analysis for study quality
care_quality <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                        data = dataCare[order(dataCare$Group == "Good/Fair", decreasing = TRUE), ],
                        sm = "SMD", fixed = FALSE, random = TRUE, method.tau = "REML", 
                        hakn = TRUE, title = "Closed-loop before-after in caregivers", 
                        subgroup = Group, test.subgroup = FALSE)

# Summary of the meta-analysis with subgroup analysis
summary(care_quality)

# Forest plot visualizing the subgroup analysis
forest.meta(care_quality, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, 
            leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)
