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


# Convert median to mean for van Bon 2024
set.seed(5555)
meansd_pre_van <- bc.mean.sd(q1.val = 18.8, med.val = 30.0, q3.val = 41.3, n = 71)
set.seed(5555)
meansd_post_van <- bc.mean.sd(q1.val = 3.8, med.val = 10.0, q3.val = 21.3, n = 71)

# Prepare data with p-values 
dataP <- data.frame(
  study = c("Schneider-Utaka 2023", "Polonsky 2022", "Bisio 2022", "Reznik 2024", "Akiyama 2023", 
            "Beato-Víbora 2021", "van Bon 2024", "Beato-Víbora 2024"),
  # when p-value < 0.00(0)1, we put 0.00(0)1
  pi = c(0.03, 0.0001, 0.046, 0.001, 0.43, 0.873, 0.0001, 0.001),
  ni = c(35, 115, 15, 202, 22, 49, 71, 150),
  m_pre = c(1.54, 1.64, 1.5, 47.3, 40.8, 2.207, meansd_pre_van$est.mean, 2.7),
  m_post = c(1.43, 1.48, 1.35, 39.9, 38.6, 2.1825, meansd_post_van$est.mean, 2.2),
  sd_pre = c(0.35, 0.51, 0.36, 18.9, 21.3, 1.1095, meansd_pre_van$est.sd, 1.3),
  sd_post = c(0.32, 0.4, 0.25, 20.8, 23.4, 0.97623, meansd_post_van$est.sd, 1.1),
  md_within = c(NA, -0.16, NA, -6.9, NA, NA, NA, NA),
  study_duration = c(3.7, 3, 1, 6, 3, 3, 12, 3),
  quality = c(3, 2, 5, 2, 4, 9, 2, 2)
)

dataPCorr <- escalc(measure="UCOR", pi = pi, ni = ni, data = dataP)
dataPCorr <- dataPCorr[,-c(2,12)]  # Remove unnecessary columns
colnames(dataPCorr)[10] <- "ri"  # Rename column for clarity

# Combine adult data
adults_ba <- rbind(dataTCorr, dataPCorr)

# Calculate raw mean differences and SMD
adults_ba$md_within <- ifelse(is.na(adults_ba$md_within), adults_ba$m_post - adults_ba$m_pre,
                               adults_ba$md_within)
adults_ba$smd_within <- adults_ba$md_within / adults_ba$sd_pre
adults_ba$se_within <- sqrt(((2 * (1 - adults_ba$ri)) / adults_ba$ni) + 
                               (adults_ba$smd_within^2 / (2 * adults_ba$ni)))

# Meta-analysis using a random-effects model
random_adults <- metagen(TE = smd_within, seTE = se_within, studlab = study, data = adults_ba, 
                        sm = "SMD", fixed = FALSE, random = TRUE, 
                        method.tau = "REML", hakn = TRUE, 
                        title = "Closed-loop before-after in adults")
summary(random_adults)

# Forest plot visualization
forest(random_adults, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, 
            leftcols = c("studlab"), leftlabs = c("Author"))


# Funnel plot visualization
funnel(random_adults, studlab = TRUE, pos.studlab = 3, col = "blue")

# High heterogeneity - outlier detection
# https://github.com/MathiasHarrer/dmetar/blob/master/R/find.outliers.R
find.outliers(random_adults) # van Bon 2024 is an outlier --> see what happens when removed


adults_ba$remove <- ifelse(adults_ba$study == "van Bon 2024", "yes", "no")

# Perform meta-analysis considering subgroup analysis for study quality
adults_out <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                     data = adults_ba, 
                     sm = "SMD", fixed = FALSE, random = TRUE, method.tau = "REML", hakn = TRUE, 
                     title = "Closed-loop before-after in adults", subgroup = remove, test.subgroup = FALSE)

# Forest plot visualizing the subgroup analysis
forest(adults_out, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)

# van Bon increases the heterogeneity from 0% to 68% --> remove!

dataAdult_out <- adults_ba[adults_ba$study!= "van Bon 2024",]

# Perform meta-analysis
random_adult_out <- metagen(TE = smd_within, seTE = se_within, studlab = study, data = dataAdult_out, 
                            sm = "SMD", fixed = FALSE, random = TRUE, 
                            method.tau = "REML", hakn = TRUE,
                            title = "Closed-loop before-after in adults")

forest(random_adult_out, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)

# Funnel plot visualization - ok
funnel(random_adult_out, studlab = TRUE, pos.studlab = 3, col = "blue")

# Categorize studies by quality
dataAdult_out$Group <- ifelse(dataAdult_out$quality < 6, "Good/Fair", "Poor")

# Perform meta-analysis considering subgroup analysis for study quality
adults_quality <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                          data = dataAdult_out[order(dataAdult_out$Group == "Good/Fair", decreasing = TRUE), ], 
                          sm = "SMD", fixed = TRUE, random = TRUE, method.tau = "REML", hakn = TRUE, 
                          title = "Closed-loop before-after in adults", 
                          subgroup = Group, test.subgroup = TRUE)

# Summary of the meta-analysis with subgroup analysis
summary(adults_quality)

# Forest plot visualizing the subgroup analysis
forest(adults_quality, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)


############################################################
#### Pediatric population ####
############################################################

# Considering Marks 2024 & Cobry 2024 & Berget 2020 at 3 months for consistency

# Convert median to mean for Cobry 2022
set.seed(5555)
meansd_pre <- bc.mean.sd(q1.val = 17, med.val = 23, q3.val = 32, n = 49)
set.seed(5555)
meansd_post <- bc.mean.sd(q1.val = 17, med.val = 21, q3.val = 31, n = 49)
# Convert median to mean for Marks 2024
set.seed(5555)
meansd_pre_Marks <- bc.mean.sd(q1.val = 27.1, med.val = 44.3, q3.val = 58.6, n = 13)
set.seed(5555)
meansd_post_Marks <- bc.mean.sd(q1.val = 21.4, med.val = 24.3, q3.val = 47.1, n = 13)

# Prepare child data

# We have t-test value for one study; we will use it to calculate r
dataChildt <- data.frame(
  study = "Cobry 2020",
  Group = "Teenagers",
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
dataChildT <- dataChildT[,-c(3,13)]  # Remove unnecessary columns
colnames(dataChildT)[11] <- "ri"  # Rename column for clarity

dataChildp <- data.frame(
  study = c("Bisio 2021", "Cobry 2022", "Gianini 2022", "Reznik 2024", "Berget 2020", "Cobry 2021", 
            "Marks 2024", "Hood 2023", "Hood 2023", "Cobry 2024"),
  Group = c("Children", "Children", "Teenagers", "Teenagers", "Both", "Children", "Both", "Children", "Teenagers", "Both"),
  # We put negative signal when effect is the opposite
  pi = c(NA, 0.153, 0.001, -0.5873, -0.91, 0.1, 0.006, 0.0010, 0.0451, -0.122),
  ni = c(13, 49, 24, 55, 92, 22, 13, 75, 38, 33),
  m_pre = c(17.26, meansd_pre$est.mean, 19.32, 33.2, 35.8, 37, meansd_pre_Marks$est.mean, 27.4, 30.5, 57.1),
  m_post = c(18.15, meansd_post$est.mean, 8.65, 34.6, 36.4, 37, meansd_post_Marks$est.mean, 24.2, 27.1, 63.1),
  # SD = SE*sqrt(n)
  sd_pre = c(13.17, meansd_pre$est.sd, 12.25, 18.6, 2.4 * sqrt(92), 17, meansd_pre_Marks$est.sd, 9.8, 11.4, 3.6*sqrt(33)),
  sd_post = c(16.22, meansd_post$est.sd, 8.26, 23.5, 3 * sqrt(92), 19, meansd_post_Marks$est.sd, 9.3, 10.8, 3.6*sqrt(33)),
  md_within = c(NA, NA, NA, 1.4, NA, NA, NA, -3.2, -3.4, NA),
  study_duration = c(1, 3.7, 4, 6, 3, 2.8, 3, 3, 3, 3),
  quality = c(7, 2, 6, 2, 0, 5, 2, 2, 2, 3)
)

# Calculate correlation coefficients
dataChildP <- escalc(measure="UCOR", pi = pi, ni = ni, data = dataChildp)
dataChildP <- dataChildP[,-c(3,13)]  # Remove unnecessary columns
colnames(dataChildP)[11] <- "ri"  # Rename column for clarity

# Combine child data
dataChild <- rbind(dataChildT, dataChildP)

# Calculate raw mean differences and standardized mean differences
dataChild$md_within <- ifelse(is.na(dataChild$md_within), dataChild$m_post - dataChild$m_pre,
                              dataChild$md_within)
dataChild$smd_within <- dataChild$md_within / dataChild$sd_pre

# Adjust correlation coefficients based on study similarities
#Bisio2021 is similar to Cobry 2020 in terms of ni and smd
dataChild$ri[2] <- dataChild$ri[1]


dataChild$se_within <- sqrt(((2 * (1 - dataChild$ri)) / dataChild$ni) + 
                              (dataChild$smd_within^2 / (2 * dataChild$ni)))

# Perform meta-analysis
random_child <- metagen(TE = smd_within, seTE = se_within, studlab = study, data = dataChild, 
                        sm = "SMD", fixed = FALSE, random = TRUE, 
                        method.tau = "REML", hakn = TRUE, 
                        title = "Closed-loop before-after in children")
summary(random_child)

# Visualize with a forest plot
forest(random_child, sortvar = dataChild$quality, prediction = TRUE, print.tau2 = FALSE, 
            leftcols = c("studlab"), leftlabs = c("Author"))


# Moderate heterogeneity - outlier detection
find.outliers(random_child) # Gianini is an outlier --> see what happens when removed


dataChild$remove <- ifelse(dataChild$study == "Gianini 2022", "yes", "no")

# Perform meta-analysis considering subgroup analysis for study quality
child_out <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                         data = dataChild, 
                         sm = "SMD", fixed = FALSE, random = TRUE, method.tau = "REML", hakn = TRUE, 
                         title = "Closed-loop before-after in children", subgroup = remove, test.subgroup = FALSE)

# Forest plot visualizing the subgroup analysis
forest(child_out, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)

# Gianini increases the heterogeneity from 6.5% to 51.7% --> remove!

dataChild_out <- dataChild[dataChild$study!= "Gianini 2022",]

# Perform meta-analysis
random_child_out <- metagen(TE = smd_within, seTE = se_within, studlab = study, data = dataChild_out, 
                            sm = "SMD", fixed = FALSE, random = TRUE, 
                            method.tau = "REML", hakn = TRUE,
                            title = "Closed-loop before-after in children")

forest(random_child_out, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)

# Funnel plot visualization - ok
funnel(random_child_out, studlab = TRUE, pos.studlab = 3, col = "blue")

metabias(random_child_out, method.bias = "egger") # p-value 0.20

# Categorize studies by quality
dataChild_out$Quality <- ifelse(dataChild_out$quality < 6, "Good/Fair", "Poor")

# Perform meta-analysis considering subgroup analysis for study quality
child_quality <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                         data = dataChild_out[order(dataChild_out$Quality == "Good/Fair", decreasing = TRUE), ], 
                         sm = "SMD", fixed = TRUE, random = TRUE, method.tau = "REML", hakn = TRUE, 
                         title = "Closed-loop before-after in children", subgroup = Quality, test.subgroup = TRUE)

# Summary of the meta-analysis with subgroup analysis
summary(child_quality)

# Forest plot visualizing the subgroup analysis
forest(child_quality, sortvar = Quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)
# Quality = Poor should be replaced by Common effects model - only 2 studies

# Perform meta-analysis considering subgroup analysis for age group
child_age <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                         data = dataChild_out[dataChild_out$Group!="Both",],
                         sm = "SMD", fixed = FALSE, random = TRUE, method.tau = "REML", hakn = TRUE, 
                         title = "Closed-loop before-after in children", subgroup = Group, test.subgroup = TRUE)

# Summary of the meta-analysis with subgroup analysis
summary(child_age)

# Forest plot visualizing the subgroup analysis
forest(child_age, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)



############################################################
#### Caregiver population  ####
############################################################

# We have t-test value for one study; we will use it to calculate r
dataCarT <- data.frame(
  study = "Cobry 2020",
  Group = "Teenagers",
  ti = 0.86,
  ni = 24,
  m_pre = 6.79,
  m_post = 6.21,
  sd_pre = 4.09,
  sd_post = 6.61,
  md_within = NA,
  #in months
  study_duration = 3,
  quality = 8
)

# Calculate correlation coefficients for t-test data
dataCarT <- escalc(measure="UCOR", ti = ti, ni=ni, data=dataCarT)
dataCarT <- dataCarT[,-c(3,13)]  # Remove unnecessary columns
colnames(dataCarT)[11] <- "ri"  # Rename column for clarity

# Convert median to mean for Cobry 2022
set.seed(5555)
meansd_pre <- bc.mean.sd(q1.val = 36, med.val = 44, q3.val = 56,n = 49)
set.seed(5555)
meansd_post <- bc.mean.sd(q1.val = 26, med.val = 35, q3.val = 43, n= 49)
# Convert median to mean for Marks 2024 - use 3 months for consistency
set.seed(5555)
meansd_pre_Marks <- bc.mean.sd(q1.val = 65.0, med.val = 73.3, q3.val = 76.3 ,n = 15)
set.seed(5555)
meansd_post_Marks <- bc.mean.sd(q1.val = 30.4, med.val = 45.0, q3.val = 67.7, n= 15)

# Prepare caregiver data
# Removing Cobry 2020 because the values do not make sense
dataCar <- data.frame(
  study = c("Bisio 2021", "Cobry 2022", "Berget 2020", "Cobry 2021", "Marks 2024", "Hood 2023", "Hood 2023", "Pulkkinen 2024",
            "MacLeish 2024", "Cobry 2024"),
  Group = c("Children", "Children", "Both", "Children", "Both", "Children", "Teenagers", "Children",
            "Children", "Both"),
  pi = c(0.032, 0.001, -0.52, NA, 0.001, 0.0001, 0.0014, 0.006, 0.0001, -0.233),
  ni = c(13, 49, 89, 22, 15, 82, 42, 35, 70, 39),
  m_pre = c(51.03, meansd_pre$est.mean, 44.2, 45, meansd_pre_Marks$est.mean, 47.1, 45.0, 37.5, 46.5, 68.2),
  m_post = c(36.85, meansd_post$est.mean, 46.2, 40, meansd_post_Marks$est.mean, 40.7, 38.0, 27.5, 37.0, 70.5),
  sd_pre = c(17.19, meansd_pre$est.sd, 2 * sqrt(89), 17, meansd_pre_Marks$est.sd, 14.3, 14.2, 18.2, 13.9, 2.6*sqrt(39)),
  sd_post = c(19.11, meansd_post$est.sd, 2.3 * sqrt(89), 17, meansd_post_Marks$est.sd, 10.8, 13.8, 14.8, 10.3, 2.6*sqrt(39)),
  md_within = c(NA, NA, NA, NA, NA, -6.3, -7.0, NA, -9.4, NA),
  study_duration = c(1, 3.7, 3, 2.8, 3, 3, 3, 3, 3, 3),
  quality = c(7, 2, 0, 5, 2, 2, 2, 1, 3, 3)
)



# Calculate correlation coefficients
dataCare <- escalc(measure="COR", pi = pi, ni = ni, data = dataCar)
dataCare <- dataCare[,-c(3,13)]  # Remove unnecessary columns
colnames(dataCare)[11] <- "ri"  # Rename column for clarity

# Combine caregiver data
dataCare <- rbind(dataCarT, dataCare)

# Calculate raw mean differences and standardized mean differences
dataCare$md_within <- ifelse(is.na(dataCare$md_within), dataCare$m_post - dataCare$m_pre, 
                             dataCare$md_within)

dataCare$smd_within <- dataCare$md_within / dataCare$sd_pre

# Adjust correlation coefficients for Cobry 2021
# similar to Cobry 2020 in children in terms of n and smd
dataCare$ri[5] <- dataCare$ri[1]

dataCare$se_within <- sqrt(((2 * (1 - dataCare$ri)) / dataCare$ni) + 
                             (dataCare$smd_within^2 / (2 * dataCare$ni)))

# Perform meta-analysis
random_care <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                       data = dataCare, sm = "SMD", fixed = FALSE, random = TRUE, 
                       method.tau = "REML", hakn = TRUE, 
                       title = "Closed-loop before-after in caregivers")
summary(random_care)

# Visualize with a forest plot
forest(random_care, sortvar = dataCare$quality, prediction = TRUE, 
            print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"))

# High heterogeneity - outliers detection
find.outliers(random_care) # Marks 2024 is an outlier (random) --> see what happens when removed

dataCare$remove <- ifelse(dataCare$study == "Marks 2024", "yes", "no")

# Perform meta-analysis considering subgroup analysis for study quality
care_out <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                     data = dataCare, 
                     sm = "SMD", fixed = FALSE, random = TRUE, method.tau = "REML", hakn = TRUE, 
                     title = "Closed-loop before-after in children", subgroup = remove, test.subgroup = FALSE)

# Forest plot visualizing
forest(care_out, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)

# Marks does not significantly influence the results or decreases heterogeneity - do not remove!

# Funnel plot visualization
funnel(random_care,  studlab = TRUE, pos.studlab = 3, col = "blue")

metabias(random_care, method.bias  = "egger", k.min = 9) # p-value 0.70


# Categorize studies by quality
dataCare$Quality <- ifelse(dataCare$quality < 6, "Good/Fair", "Poor")

# Perform meta-analysis considering subgroup analysis for study quality
care_quality <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                        data = dataCare[order(dataCare$Quality == "Good/Fair", decreasing = TRUE), ],
                        sm = "SMD", fixed = TRUE, random = TRUE, method.tau = "REML", 
                        hakn = TRUE, title = "Closed-loop before-after in caregivers", 
                        subgroup = Quality, test.subgroup = TRUE)

# Summary of the meta-analysis with subgroup analysis
summary(care_quality)

# Forest plot visualizing the subgroup analysis
forest(care_quality, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, 
            leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)


# Perform meta-analysis considering subgroup analysis for age group
care_ba_age <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                        data = dataCare[dataCare$Group!= "Both",], 
                        sm = "SMD", fixed = TRUE, random = TRUE, method.tau = "REML", hakn = TRUE,
                        subgroup = Group, test.subgroup = TRUE)

# Summary of the meta-analysis with subgroup analysis
summary(care_ba_age)

# Forest plot visualizing the subgroup analysis
forest(care_ba_age, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)


# Perform meta-analysis removing group = Both

# Categorize studies by quality
dataCare$group <- ifelse(dataCare$Group == "Both", "Both", "Clear")

care_ba_both <- metagen(TE = smd_within, seTE = se_within, studlab = study, 
                       data = dataCare,
                       sm = "SMD", fixed = FALSE, random = TRUE, method.tau = "REML", hakn = TRUE,
                       subgroup = group, test.subgroup = TRUE)

# Summary of the meta-analysis 
summary(care_ba_both)

# Forest plot 
forest(care_ba_both, sortvar = quality, prediction = TRUE, print.tau2 = FALSE, leftcols = c("studlab"), leftlabs = c("Author"), print.Q.subgroup = FALSE)
