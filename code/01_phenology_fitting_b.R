library(chillR)
library(tidyverse)
library(patchwork)

# For reproducibility, this part will set the seed for pseudo random number generator. In theory, I will run this code chunk
# only in the final version of the analysis. In the meanwhile, results will differ every time. Use system time to avoid
# anchoring the seed selection
set.seed(52243)

# This script will run the PhenoFlex approach using the phenology data collected in our experiment at campus Klein-Altendorf.
# Since the data for pears only contains 32 different treatments, I will only analyze the data for apples (66 seasons).
# I will implement two versions of the analysis. The first one will use all experimental seasons minus 2 seasons that 
# resulted in no full bloom (extreme warm treatments). In the second version I will remove 6 low chill treatments that will
# rarely be observed in normal conditions (?), with some of them producing controversial blooming.

# Load the data from the folder
data <- read.csv("data/final_weather_data_S1_S2_apple_hourly.csv")

# Generate a new column (Year_2) to simulate the year and comply with the format of PhenoFlex functions
data["Year_2"] <- data$Treatment + data$Year 

# Since this experiment was conducted during two consecutive seasons, the next step will fix a small continuity issue
# generated during the season 2
data[data$Treatment >= 34, "Year_2"] <- data[data$Treatment >= 34, "Year_2"] - 1

# For further compatibility, I will now select the columns needed and will drop "Year" (the original one)
data <- data[c("YEARMODA", "Year_2", "Month", "Day", "Hour", "JDay", "Temp")]

# To replace the missing "Year" column, I will now change the name of the column
colnames(data)[which(colnames(data) == "Year_2")] <- "Year"


# Import the phenology data from the repository
pheno <- read.csv("data/final_bio_data_S1_S2_apple.csv")

# Generate two data sets according to the version of the analysis. I will also exclude repeated treatments 22, 23, 25, and 26
pheno_v1 <- pheno[!(pheno$Treatment %in% c(36, 3, 23, 24, 17, 18, 61)), ]
pheno_v2 <- pheno[!(pheno$Treatment %in% c(36, 3, 23, 24, 17, 18, 61, # Excluding these treatments from both datasets
                                           13, 42, 46, 58, 62)), ] # These treatments are the marginal seasons (n = 5)

# Add the year column to match the data in the weather data frame
pheno_v1["Treatment"] <- pheno_v1$Treatment + 2019
pheno_v2["Treatment"] <- pheno_v2$Treatment + 2019

# Select only the relevant columns for further analysis
pheno_v1 <- pheno_v1[c("Treatment", "pheno")]
pheno_v2 <- pheno_v2[c("Treatment", "pheno")]

# Rename the columns for further compatibility
colnames(pheno_v1) <- c("Year", "pheno")
colnames(pheno_v2) <- c("Year", "pheno")


# Define a vector of calibration and validation seasons. V1 includes the marginal seasons
calibration_seasons <- sort(sample(pheno_v2$Year, 40, replace = FALSE))
calibration_seasons_v1 <- sort(c(sample(calibration_seasons, 35, replace = FALSE),
                                 pheno_v1$Year[which(!(pheno_v1$Year %in% pheno_v2$Year))]))
calibration_seasons_v2 <- calibration_seasons

# Common validation seasons
validation_seasons <- sort(pheno_v2[!(pheno_v2$Year %in% calibration_seasons), "Year"])

# Define the list of seasons (weather data)
season_list_v1 <- genSeasonList(data, mrange = c(9, 7), years = calibration_seasons_v1)
season_list_v2 <- genSeasonList(data, mrange = c(9, 7), years = calibration_seasons_v2)



# Model run number 1 ####
# Set the initial parameters (wide ranges)
#          yc,  zc,  s1, Tu,     E0,      E1,     A0,          A1,   Tf, Tc, Tb, slope
lower <- c(20, 100, 0.1,  0, 3000.0,  9000.0, 6000.0,       5.e13,    0,  0,  0,  0.05)
par   <- c(40, 190, 0.5, 25, 3372.8,  9900.3, 6319.5, 5.939917e13,    4, 36,  4,  1.60)
upper <- c(80, 500, 1.0, 30, 4000.0, 10000.0, 7000.0,       6.e13,   10, 40, 10, 50.00)

# Run the fitter
pheno_fit_v1_r1 <- phenologyFitter(par.guess = par,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower,
                                   upper = upper,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1000,
                                                  nb.stop.improvement = 5))

# Same for version 2
pheno_fit_v2_r1 <- phenologyFitter(par.guess = par,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower,
                                   upper = upper,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1000,
                                                  nb.stop.improvement = 5))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r1 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r1 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r1[["Predicted"]] <- pheno_fit_v1_r1$pbloomJDays
out_df_v2_r1[["Predicted"]] <- pheno_fit_v2_r1$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r1[["Error"]] <- pheno_fit_v1_r1$bloomJDays - pheno_fit_v1_r1$pbloomJDays
out_df_v2_r1[["Error"]] <- pheno_fit_v2_r1$bloomJDays - pheno_fit_v2_r1$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r1 <- RMSEP(out_df_v1_r1$Predicted, out_df_v1_r1$pheno)
RMSEP_calib_v2_r1 <- RMSEP(out_df_v2_r1$Predicted, out_df_v2_r1$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r1 <- RPIQ(out_df_v1_r1$Predicted, out_df_v1_r1$pheno)
RPIQ_calib_v2_r1 <- RPIQ(out_df_v2_r1$Predicted, out_df_v2_r1$pheno)

# Estimate the mean error
mean(out_df_v1_r1$Error)
mean(out_df_v2_r1$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r1$Error))
mean(abs(out_df_v2_r1$Error))

# Compute the AICc
aic_fit_v1_r1 <- (2 * length(par)) + (length(calibration_seasons_v1) * 
                                        log(pheno_fit_v1_r1$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par)^2 + 2 * length(par)) /
                                                                                                                    length(calibration_seasons_v1) - length(par) - 1)
aic_fit_v2_r1 <- (2 * length(par)) + (length(calibration_seasons_v2) * 
                                        log(pheno_fit_v2_r1$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par)^2 + 2 * length(par)) /
                                                                                                                    length(calibration_seasons_v2) - length(par) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r1, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r1, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")



# Model run number 2 ####
# Set the parameters using the results from the previous run
# Version 1 (pheno_fit_v1_r1$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v1_r2 <- c(          45,          290,            0,           18,       3000.0,       9000.0,       6000.0,       5.2e13,            0,           29,            0,         0.05)
par_v1_r2   <- c(5.898501e+01, 3.260215e+02, 1.606328e-01, 2.251842e+01, 3.309868e+03, 9.900106e+03, 6.295673e+03, 5.939926e+13, 6.143360e+00, 3.986939e+01, 8.630051e+00, 1.459628e+01)
upper_v1_r2 <- c(          65,          380,          0.5,           30,       4000.0,      10000.0,       7000.0,       6.2e13,           10,           44,           12,        30.00)

# Run the fitter
pheno_fit_v1_r2 <- phenologyFitter(par.guess = par_v1_r2,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r2,
                                   upper = upper_v1_r2,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1000,
                                                  nb.stop.improvement = 5))

# Same for version 2 (pheno_fit_v2_r1$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r2 <- c(          25,          340,            0,           17,       3000.0,       9000.0,       6000.0,       5.2e13,            0,           30,            0,         1.00)
par_v2_r2   <- c(3.492167e+01, 3.872746e+02, 2.715607e-01, 2.176161e+01, 3.371032e+03, 9.900773e+03, 6.283549e+03, 5.939925e+13, 5.119899e+00, 3.999171e+01, 6.298089e+00, 3.999151e+01)
upper_v2_r2 <- c(          40,          420,          0.5,           32,       4000.0,      10000.0,       7000.0,       6.2e13,           10,           42,           10,        50.00)

# Run the fitter
pheno_fit_v2_r2 <- phenologyFitter(par.guess = par_v2_r2,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r2,
                                   upper = upper_v2_r2,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1000,
                                                  nb.stop.improvement = 5))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r2 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r2 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r2[["Predicted"]] <- pheno_fit_v1_r2$pbloomJDays
out_df_v2_r2[["Predicted"]] <- pheno_fit_v2_r2$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r2[["Error"]] <- pheno_fit_v1_r2$bloomJDays - pheno_fit_v1_r2$pbloomJDays
out_df_v2_r2[["Error"]] <- pheno_fit_v2_r2$bloomJDays - pheno_fit_v2_r2$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r2 <- RMSEP(out_df_v1_r2$Predicted, out_df_v1_r2$pheno)
RMSEP_calib_v2_r2 <- RMSEP(out_df_v2_r2$Predicted, out_df_v2_r2$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r2 <- RPIQ(out_df_v1_r2$Predicted, out_df_v1_r2$pheno)
RPIQ_calib_v2_r2 <- RPIQ(out_df_v2_r2$Predicted, out_df_v2_r2$pheno)

# Estimate the mean error
mean(out_df_v1_r2$Error)
mean(out_df_v2_r2$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r2$Error))
mean(abs(out_df_v2_r2$Error))

# Compute the AICc
aic_fit_v1_r2 <- (2 * length(par_v1_r2)) + (length(calibration_seasons_v1) * 
                                              log(pheno_fit_v1_r2$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par_v1_r2)^2 + 2 * length(par_v1_r2)) /
                                                                                                                          length(calibration_seasons_v1) - length(par_v1_r2) - 1)
aic_fit_v2_r2 <- (2 * length(par_v2_r2)) + (length(calibration_seasons_v2) * 
                                              log(pheno_fit_v2_r2$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par_v2_r2)^2 + 2 * length(par_v2_r2)) /
                                                                                                                          length(calibration_seasons_v2) - length(par_v2_r2) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r2, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r2, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")



# Model run number 3 ####
# Set the parameters using the results from the previous run
# Version 1 (pheno_fit_v1_r2$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v1_r3 <- c(          48,          290,            0,           15,       3000.0,       9000.0,       6000.0,       5.2e13,            0,           24,            0,         0.05)
par_v1_r3   <- c(5.837468e+01, 3.322524e+02, 3.899379e-01, 2.204617e+01, 3.310447e+03, 9.900260e+03, 6.344819e+03, 5.939986e+13, 7.235736e+00, 3.467449e+01, 8.883098e+00, 1.817915e+01)
upper_v1_r3 <- c(          68,          370,         0.65,           32,       4000.0,      10500.0,       7000.0,       6.2e13,           11,           44,           11,        35.00)

# Run the fitter
pheno_fit_v1_r3 <- phenologyFitter(par.guess = par_v1_r3,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r3,
                                   upper = upper_v1_r3,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1000,
                                                  nb.stop.improvement = 5))

# Same for version 2 (pheno_fit_v2_r2$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r3 <- c(          20,          350,            0,           18,       3000.0,       9000.0,       6000.0,       5.2e13,            0,           30,            0,         1.00)
par_v2_r3   <- c(3.379713e+01, 4.062356e+02, 2.361878e-01, 2.137393e+01, 3.371039e+03, 9.901300e+03, 6.214008e+03, 5.939925e+13, 1.454457e+00, 4.196362e+01, 6.094298e+00, 1.070739e+01)
upper_v2_r3 <- c(          50,          450,         0.55,           35,       4000.0,      10500.0,       7000.0,       6.2e13,           10,           46,           10,        35.00)

# Run the fitter
pheno_fit_v2_r3 <- phenologyFitter(par.guess = par_v2_r3,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r3,
                                   upper = upper_v2_r3,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1500,
                                                  nb.stop.improvement = 10))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r3 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r3 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r3[["Predicted"]] <- pheno_fit_v1_r3$pbloomJDays
out_df_v2_r3[["Predicted"]] <- pheno_fit_v2_r3$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r3[["Error"]] <- pheno_fit_v1_r3$bloomJDays - pheno_fit_v1_r3$pbloomJDays
out_df_v2_r3[["Error"]] <- pheno_fit_v2_r3$bloomJDays - pheno_fit_v2_r3$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r3 <- RMSEP(out_df_v1_r3$Predicted, out_df_v1_r3$pheno)
RMSEP_calib_v2_r3 <- RMSEP(out_df_v2_r3$Predicted, out_df_v2_r3$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r3 <- RPIQ(out_df_v1_r3$Predicted, out_df_v1_r3$pheno)
RPIQ_calib_v2_r3 <- RPIQ(out_df_v2_r3$Predicted, out_df_v2_r3$pheno)

# Estimate the mean error
mean(out_df_v1_r3$Error)
mean(out_df_v2_r3$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r3$Error))
mean(abs(out_df_v2_r3$Error))

# Compute the AICc
aic_fit_v1_r3 <- (2 * length(par_v1_r3)) + (length(calibration_seasons_v1) * 
                                              log(pheno_fit_v1_r3$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par_v1_r3)^2 + 2 * length(par_v1_r3)) /
                                                                                                                          length(calibration_seasons_v1) - length(par_v1_r3) - 1)
aic_fit_v2_r3 <- (2 * length(par_v2_r3)) + (length(calibration_seasons_v2) * 
                                              log(pheno_fit_v2_r3$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par_v2_r3)^2 + 2 * length(par_v2_r3)) /
                                                                                                                          length(calibration_seasons_v2) - length(par_v2_r3) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r3, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r3, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")



# Model run number 4 ####
# Set the parameters using the results from the previous run
# Version 1 (pheno_fit_v1_r3$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v1_r4 <- c(          45,          290,            0,           15,       3000.0,       9000.0,       5500.0,       5.2e13,            0,           20,            0,            1)
par_v1_r4   <- c(5.982188e+01, 3.391820e+02, 5.362487e-01, 2.298870e+01, 3.310334e+03, 9.901517e+03, 6.346304e+03, 5.939996e+13, 7.197365e+00, 3.069973e+01, 7.613608e+00, 7.334750e+00)
upper_v1_r4 <- c(          70,          370,         0.75,           32,       4000.0,      10300.0,       7000.0,       6.2e13,           12,           40,           12,        20.00)

# Run the fitter
pheno_fit_v1_r4 <- phenologyFitter(par.guess = par_v1_r4,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r4,
                                   upper = upper_v1_r4,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1000,
                                                  nb.stop.improvement = 5))

# Same for version 2 (pheno_fit_v2_r3$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r4 <- c(          28,          360,            0,           17,       3000.0,       9300.0,       5800.0,       5.3e13,            0,           35,            0,           10)
par_v2_r4   <- c(3.412481e+01, 3.962168e+02, 1.958515e-01, 2.161528e+01, 3.371039e+03, 9.901299e+03, 6.214040e+03, 5.939950e+13, 1.454457e+00, 4.599390e+01, 5.920062e+00, 1.112166e+01)
upper_v2_r4 <- c(          44,          440,         0.45,           30,       4000.0,      10500.0,       6800.0,       6.5e13,           10,           48,           12,        40.00)

# Run the fitter
pheno_fit_v2_r4 <- phenologyFitter(par.guess = par_v2_r4,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r4,
                                   upper = upper_v2_r4,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1000,
                                                  nb.stop.improvement = 5))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r4 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r4 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r4[["Predicted"]] <- pheno_fit_v1_r4$pbloomJDays
out_df_v2_r4[["Predicted"]] <- pheno_fit_v2_r4$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r4[["Error"]] <- pheno_fit_v1_r4$bloomJDays - pheno_fit_v1_r4$pbloomJDays
out_df_v2_r4[["Error"]] <- pheno_fit_v2_r4$bloomJDays - pheno_fit_v2_r4$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r4 <- RMSEP(out_df_v1_r4$Predicted, out_df_v1_r4$pheno)
RMSEP_calib_v2_r4 <- RMSEP(out_df_v2_r4$Predicted, out_df_v2_r4$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r4 <- RPIQ(out_df_v1_r4$Predicted, out_df_v1_r4$pheno)
RPIQ_calib_v2_r4 <- RPIQ(out_df_v2_r4$Predicted, out_df_v2_r4$pheno)

# Estimate the mean error
mean(out_df_v1_r4$Error)
mean(out_df_v2_r4$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r4$Error))
mean(abs(out_df_v2_r4$Error))

# Compute the AICc
aic_fit_v1_r4 <- (2 * length(par_v1_r4)) + (length(calibration_seasons_v1) * 
                                              log(pheno_fit_v1_r4$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par_v1_r4)^2 + 2 * length(par_v1_r4)) /
                                                                                                                          length(calibration_seasons_v1) - length(par_v1_r4) - 1)
aic_fit_v2_r4 <- (2 * length(par_v2_r4)) + (length(calibration_seasons_v2) * 
                                              log(pheno_fit_v2_r4$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par_v2_r4)^2 + 2 * length(par_v2_r4)) /
                                                                                                                          length(calibration_seasons_v2) - length(par_v2_r4) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r4, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r4, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")



# Model run number 5 ####
# Set the parameters using the results from the previous run
# Version 1 (pheno_fit_v1_r4$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v1_r5 <- c(          50,          330,            0,           18,       3000.0,       9600.0,       6000.0,       5.5e13,            2,           30,            2,          0.5)
par_v1_r5   <- c(6.002517e+01, 3.538325e+02, 4.668646e-01, 2.279862e+01, 3.310305e+03, 9.901535e+03, 6.344317e+03, 5.939992e+13, 7.145970e+00, 3.601947e+01, 7.611601e+00, 7.627144e+00)
upper_v1_r5 <- c(          75,          370,         0.55,           28,       4000.0,      10300.0,       7000.0,       6.5e13,           10,           42,           12,        15.00)

# Run the fitter
pheno_fit_v1_r5 <- phenologyFitter(par.guess = par_v1_r5,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r5,
                                   upper = upper_v1_r5,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1500,
                                                  nb.stop.improvement = 10))

# Same for version 2 (pheno_fit_v2_r4$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r5 <- c(          27,          360,            0,           18,       3000.0,       9500.0,       5800.0,       5.5e13,            0,           35,            1,         1.00)
par_v2_r5   <- c(3.398745e+01, 3.832050e+02, 2.058296e-01, 2.245362e+01, 3.371035e+03, 9.901300e+03, 6.213937e+03, 5.939939e+13, 1.454543e+00, 4.555801e+01, 5.917282e+00, 1.081221e+01)
upper_v2_r5 <- c(          37,          400,         0.55,           32,       4000.0,      10300.0,       6800.0,       6.3e13,            8,           48,           10,        25.00)

# Run the fitter
pheno_fit_v2_r5 <- phenologyFitter(par.guess = par_v2_r5,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r5,
                                   upper = upper_v2_r5,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1000,
                                                  nb.stop.improvement = 5))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r5 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r5 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r5[["Predicted"]] <- pheno_fit_v1_r5$pbloomJDays
out_df_v2_r5[["Predicted"]] <- pheno_fit_v2_r5$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r5[["Error"]] <- pheno_fit_v1_r5$bloomJDays - pheno_fit_v1_r5$pbloomJDays
out_df_v2_r5[["Error"]] <- pheno_fit_v2_r5$bloomJDays - pheno_fit_v2_r5$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r5 <- RMSEP(out_df_v1_r5$Predicted, out_df_v1_r5$pheno)
RMSEP_calib_v2_r5 <- RMSEP(out_df_v2_r5$Predicted, out_df_v2_r5$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r5 <- RPIQ(out_df_v1_r5$Predicted, out_df_v1_r5$pheno)
RPIQ_calib_v2_r5 <- RPIQ(out_df_v2_r5$Predicted, out_df_v2_r5$pheno)

# Estimate the mean error
mean(out_df_v1_r5$Error)
mean(out_df_v2_r5$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r5$Error))
mean(abs(out_df_v2_r5$Error))

# Compute the AICc
aic_fit_v1_r5 <- (2 * length(par_v1_r5)) + (length(calibration_seasons_v1) * 
                                              log(pheno_fit_v1_r5$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par_v1_r5)^2 + 2 * length(par_v1_r5)) /
                                                                                                                          length(calibration_seasons_v1) - length(par_v1_r5) - 1)
aic_fit_v2_r5 <- (2 * length(par_v2_r5)) + (length(calibration_seasons_v2) * 
                                              log(pheno_fit_v2_r5$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par_v2_r5)^2 + 2 * length(par_v2_r5)) /
                                                                                                                          length(calibration_seasons_v2) - length(par_v2_r5) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r5, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r5, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")



# Model run number 6 ####
# Set the parameters using the results from the previous run
# Version 1 (pheno_fit_v1_r5$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v1_r6 <- c(          40,          290,            0,           18,       2900.0,       9000.0,       5700.0,       5.3e13,            0,           20,            0,         0.01)
par_v1_r6   <- c(6.003778e+01, 3.430605e+02, 3.238688e-01, 2.507193e+01, 3.310284e+03, 9.901609e+03, 6.342868e+03, 5.940006e+13, 5.636596e+00, 3.360488e+01, 6.974146e+00, 2.238448e+00)
upper_v1_r6 <- c(          80,          390,         0.55,           32,       3900.0,      10500.0,       6900.0,       6.8e13,           12,           48,           12,        20.00)


# Run the fitter
pheno_fit_v1_r6 <- phenologyFitter(par.guess = par_v1_r6,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r6,
                                   upper = upper_v1_r6,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 2000,
                                                  nb.stop.improvement = 100))

# Same for version 2 (pheno_fit_v2_r5$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r6 <- c(          25,          330,            0,           18,       3000.0,       9200.0,       5700.0,       5.4e13,         0.05,           40,            0,          0.5)
par_v2_r6   <- c(3.467235e+01, 3.613910e+02, 1.481830e-01, 2.242230e+01, 3.371040e+03, 9.901317e+03, 6.213769e+03, 5.939943e+13, 1.453357e+00, 4.798370e+01, 5.917437e+00, 1.100240e+01)
upper_v2_r6 <- c(          40,          390,         0.35,           30,       4000.0,      10800.0,       6900.0,       6.6e13,           10,           50,           10,        20.00)

# Run the fitter
pheno_fit_v2_r6 <- phenologyFitter(par.guess = par_v2_r6,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r6,
                                   upper = upper_v2_r6,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 2000,
                                                  nb.stop.improvement = 100))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r6 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r6 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r6[["Predicted"]] <- pheno_fit_v1_r6$pbloomJDays
out_df_v2_r6[["Predicted"]] <- pheno_fit_v2_r6$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r6[["Error"]] <- pheno_fit_v1_r6$bloomJDays - pheno_fit_v1_r6$pbloomJDays
out_df_v2_r6[["Error"]] <- pheno_fit_v2_r6$bloomJDays - pheno_fit_v2_r6$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r6 <- RMSEP(out_df_v1_r6$Predicted, out_df_v1_r6$pheno)
RMSEP_calib_v2_r6 <- RMSEP(out_df_v2_r6$Predicted, out_df_v2_r6$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r6 <- RPIQ(out_df_v1_r6$Predicted, out_df_v1_r6$pheno)
RPIQ_calib_v2_r6 <- RPIQ(out_df_v2_r6$Predicted, out_df_v2_r6$pheno)

# Estimate the mean error
mean(out_df_v1_r6$Error)
mean(out_df_v2_r6$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r6$Error))
mean(abs(out_df_v2_r6$Error))

# Compute the AICc
aic_fit_v1_r6 <- (2 * length(par_v1_r6)) + (length(calibration_seasons_v1) * 
                                              log(pheno_fit_v1_r6$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par_v1_r6)^2 + 2 * length(par_v1_r6)) /
                                                                                                                          length(calibration_seasons_v1) - length(par_v1_r6) - 1)
aic_fit_v2_r6 <- (2 * length(par_v2_r6)) + (length(calibration_seasons_v2) * 
                                              log(pheno_fit_v2_r6$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par_v2_r6)^2 + 2 * length(par_v2_r6)) /
                                                                                                                          length(calibration_seasons_v2) - length(par_v2_r6) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r6, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r6, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")



# Model run number 7 ####
# Set the parameters using the results from the previous run
# Version 1 (pheno_fit_v1_r6$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v1_r7 <- c(          55,          325,            0,           20,       3000.0,       9000.0,       5500.0,       5.5e13,            0,           30,            2,            0)
par_v1_r7   <- c(6.075154e+01, 3.341590e+02, 5.318899e-01, 2.495407e+01, 3.310303e+03, 9.901636e+03, 6.342682e+03, 5.940032e+13, 6.418330e+00, 3.388484e+01, 6.956173e+00, 1.573674e+00)
upper_v1_r7 <- c(          65,          345,         0.65,           30,       4000.0,      10500.0,       6600.0,       6.5e13,           10,           40,           10,        10.00)

# Run the fitter
pheno_fit_v1_r7 <- phenologyFitter(par.guess = par_v1_r7,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r7,
                                   upper = upper_v1_r7,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 3000,
                                                  nb.stop.improvement = 50))

# Same for version 2 (pheno_fit_v2_r6$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r7 <- c(          30,          365,            0,           20,       3000.0,       9500.0,       5800.0,       5.5e13,            0,           35,            0,          1.0)
par_v2_r7   <- c(3.416959e+01, 3.796795e+02, 1.588349e-01, 2.637063e+01, 3.371036e+03, 9.901315e+03, 6.214195e+03, 5.939915e+13, 1.268310e+00, 4.119790e+01, 3.340906e+00, 1.120830e+01)
upper_v2_r7 <- c(          40,          385,        0.225,           30,       4000.0,      10500.0,       6800.0,       6.4e13,            6,           45,            8,        20.00)

# Run the fitter
pheno_fit_v2_r7 <- phenologyFitter(par.guess = par_v2_r7,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r7,
                                   upper = upper_v2_r7,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 3000,
                                                  nb.stop.improvement = 50))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r7 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r7 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r7[["Predicted"]] <- pheno_fit_v1_r7$pbloomJDays
out_df_v2_r7[["Predicted"]] <- pheno_fit_v2_r7$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r7[["Error"]] <- pheno_fit_v1_r7$bloomJDays - pheno_fit_v1_r7$pbloomJDays
out_df_v2_r7[["Error"]] <- pheno_fit_v2_r7$bloomJDays - pheno_fit_v2_r7$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r7 <- RMSEP(out_df_v1_r7$Predicted, out_df_v1_r7$pheno)
RMSEP_calib_v2_r7 <- RMSEP(out_df_v2_r7$Predicted, out_df_v2_r7$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r7 <- RPIQ(out_df_v1_r7$Predicted, out_df_v1_r7$pheno)
RPIQ_calib_v2_r7 <- RPIQ(out_df_v2_r7$Predicted, out_df_v2_r7$pheno)

# Estimate the mean error
mean(out_df_v1_r7$Error)
mean(out_df_v2_r7$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r7$Error))
mean(abs(out_df_v2_r7$Error))

# Compute the AICc
aic_fit_v1_r7 <- (2 * length(par_v1_r7)) + (length(calibration_seasons_v1) * 
                                              log(pheno_fit_v1_r7$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par_v1_r7)^2 + 2 * length(par_v1_r7)) /
                                                                                                                          length(calibration_seasons_v1) - length(par_v1_r7) - 1)
aic_fit_v2_r7 <- (2 * length(par_v2_r7)) + (length(calibration_seasons_v2) * 
                                              log(pheno_fit_v2_r7$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par_v2_r7)^2 + 2 * length(par_v2_r7)) /
                                                                                                                          length(calibration_seasons_v2) - length(par_v2_r7) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r7, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r7, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")


# Model run number 8 ####
# Set the parameters using the results from the previous run
# Version 1 (pheno_fit_v1_r7$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v1_r8 <- c(          55,          315,          0.3,           22,       2800.0,       9500.0,       5900.0,       5.6e13,            2,           26,            2,            0)
par_v1_r8   <- c(6.241611e+01, 3.363827e+02, 5.244054e-01, 2.711304e+01, 3.310309e+03, 9.901660e+03, 6.395421e+03, 5.940064e+13, 6.474790e+00, 3.023036e+01, 5.553607e+00, 1.386503e+00)
upper_v1_r8 <- c(          70,          355,          0.7,           32,       3800.0,      10200.0,       6700.0,       6.4e13,            9,           36,            8,        15.00)

# Run the fitter
pheno_fit_v1_r8 <- phenologyFitter(par.guess = par_v1_r8,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r8,
                                   upper = upper_v1_r8,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 4000,
                                                  nb.stop.improvement = 150))

# Same for version 2 (pheno_fit_v2_r7$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r8 <- c(          28,          345,         0.01,           20,       2800.0,       9600.0,       5800.0,       5.6e13,            0,           38,            0,            7)
par_v2_r8   <- c(3.471534e+01, 3.747534e+02, 1.284166e-01, 2.636967e+01, 3.371029e+03, 9.901260e+03, 6.214426e+03, 5.939895e+13, 1.330276e+00, 4.472405e+01, 2.962285e+00, 1.795351e+01)
upper_v2_r8 <- c(          40,          405,          0.3,           32,       3700.0,      10300.0,       6800.0,       6.3e13,            8,           52,            8,        30.00)

# Run the fitter
pheno_fit_v2_r8 <- phenologyFitter(par.guess = par_v2_r8,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r8,
                                   upper = upper_v2_r8,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 4000,
                                                  nb.stop.improvement = 150))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r8 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r8 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r8[["Predicted"]] <- pheno_fit_v1_r8$pbloomJDays
out_df_v2_r8[["Predicted"]] <- pheno_fit_v2_r8$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r8[["Error"]] <- pheno_fit_v1_r8$bloomJDays - pheno_fit_v1_r8$pbloomJDays
out_df_v2_r8[["Error"]] <- pheno_fit_v2_r8$bloomJDays - pheno_fit_v2_r8$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r8 <- RMSEP(out_df_v1_r8$Predicted, out_df_v1_r8$pheno)
RMSEP_calib_v2_r8 <- RMSEP(out_df_v2_r8$Predicted, out_df_v2_r8$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r8 <- RPIQ(out_df_v1_r8$Predicted, out_df_v1_r8$pheno)
RPIQ_calib_v2_r8 <- RPIQ(out_df_v2_r8$Predicted, out_df_v2_r8$pheno)

# Estimate the mean error
mean(out_df_v1_r8$Error)
mean(out_df_v2_r8$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r8$Error))
mean(abs(out_df_v2_r8$Error))

# Compute the AICc
aic_fit_v1_r8 <- (2 * length(par_v1_r8)) + (length(calibration_seasons_v1) * 
                                              log(pheno_fit_v1_r8$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par_v1_r8)^2 + 2 * length(par_v1_r8)) /
                                                                                                                          length(calibration_seasons_v1) - length(par_v1_r8) - 1)
aic_fit_v2_r8 <- (2 * length(par_v2_r8)) + (length(calibration_seasons_v2) * 
                                              log(pheno_fit_v2_r8$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par_v2_r8)^2 + 2 * length(par_v2_r8)) /
                                                                                                                          length(calibration_seasons_v2) - length(par_v2_r8) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r8, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r8, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")


# Model run number 9 ####
# Set the parameters using the results from the previous run
# Version 1 (pheno_fit_v1_r8$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v1_r9 <- c(          55,          300,          0.5,           20,       3000.0,       9500.0,       5500.0,       5.0e13,            1,           20,            2,            0)
par_v1_r9   <- c(6.308038e+01, 3.150055e+02, 6.960001e-01, 2.789891e+01, 3.310344e+03, 9.901652e+03, 6.396146e+03, 5.940060e+13, 6.469237e+00, 2.789990e+01, 5.700527e+00, 1.392844e+00)
upper_v1_r9 <- c(          75,          325,         0.85,           35,       4000.0,      10500.0,       6900.0,       7.0e13,           11,           37,           10,        25.00)

# Run the fitter
pheno_fit_v1_r9 <- phenologyFitter(par.guess = par_v1_r9,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r9,
                                   upper = upper_v1_r9,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 5000,
                                                  nb.stop.improvement = 300))

# Same for version 2 (pheno_fit_v2_r8$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r9 <- c(          23,          325,            0,           16,       3000.0,       9500.0,       5700.0,       5.6e13,            0,           40,            0,            1)
par_v2_r9   <- c(3.398599e+01, 3.567737e+02, 1.355020e-01, 2.684285e+01, 3.371005e+03, 9.901244e+03, 6.214071e+03, 5.939792e+13, 1.373655e+00, 5.133324e+01, 3.632503e+00, 8.446999e+00)
upper_v2_r9 <- c(          43,          375,          0.5,           36,       3800.0,      10500.0,       6700.0,       6.3e13,           10,           55,           10,        20.00)

# Run the fitter
pheno_fit_v2_r9 <- phenologyFitter(par.guess = par_v2_r9,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r9,
                                   upper = upper_v2_r9,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 5000,
                                                  nb.stop.improvement = 200))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r9 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r9 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r9[["Predicted"]] <- pheno_fit_v1_r9$pbloomJDays
out_df_v2_r9[["Predicted"]] <- pheno_fit_v2_r9$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r9[["Error"]] <- pheno_fit_v1_r9$bloomJDays - pheno_fit_v1_r9$pbloomJDays
out_df_v2_r9[["Error"]] <- pheno_fit_v2_r9$bloomJDays - pheno_fit_v2_r9$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r9 <- RMSEP(out_df_v1_r9$Predicted, out_df_v1_r9$pheno)
RMSEP_calib_v2_r9 <- RMSEP(out_df_v2_r9$Predicted, out_df_v2_r9$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r9 <- RPIQ(out_df_v1_r9$Predicted, out_df_v1_r9$pheno)
RPIQ_calib_v2_r9 <- RPIQ(out_df_v2_r9$Predicted, out_df_v2_r9$pheno)

# Estimate the mean error
mean(out_df_v1_r9$Error)
mean(out_df_v2_r9$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r9$Error))
mean(abs(out_df_v2_r9$Error))

# Compute the AICc
aic_fit_v1_r9 <- (2 * length(par_v1_r9)) + (length(calibration_seasons_v1) * 
                                              log(pheno_fit_v1_r9$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par_v1_r9)^2 + 2 * length(par_v1_r9)) /
                                                                                                                          length(calibration_seasons_v1) - length(par_v1_r9) - 1)
aic_fit_v2_r9 <- (2 * length(par_v2_r9)) + (length(calibration_seasons_v2) * 
                                              log(pheno_fit_v2_r9$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par_v2_r9)^2 + 2 * length(par_v2_r9)) /
                                                                                                                          length(calibration_seasons_v2) - length(par_v2_r9) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r9, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r9, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")



# Model run number 10 ####
# Set the parameters using the results from the previous run
# Version 1 (pheno_fit_v1_r9$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v1_r10 <- c(          58,          305,         0.15,           20,       2800.0,       9500.0,       5800.0,       5.7e13,            0,           20,            0,            0)
par_v1_r10   <- c(6.327477e+01, 3.190987e+02, 8.449871e-01, 2.788848e+01, 3.310345e+03, 9.901638e+03, 6.396174e+03, 5.940024e+13, 6.474602e+00, 2.790773e+01, 5.591152e+00, 1.393939e+00)
upper_v1_r10 <- c(          68,          335,         1.05,           35,       3700.0,      10400.0,       6600.0,       6.2e13,           12,           34,           10,        15.00)

# Run the fitter
pheno_fit_v1_r10 <- phenologyFitter(par.guess = par_v1_r10,
                                    modelfn = PhenoFlex_GDHwrapper,
                                    bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                    SeasonList = season_list_v1,
                                    lower = lower_v1_r10,
                                    upper = upper_v1_r10,
                                    control = list(smooth = FALSE,
                                                   verbose = FALSE,
                                                   maxit = 10000,
                                                   nb.stop.improvement = 500))

# Same for version 2 (pheno_fit_v2_r9$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r10 <- c(          28,          360,         0.05,           20,       2900.0,       9500.0,       5800.0,       5.7e13,            0,           45,            0,            0)
par_v2_r10   <- c(3.370071e+01, 3.709721e+02, 1.332810e-01, 2.495372e+01, 3.371002e+03, 9.901248e+03, 6.214569e+03, 5.939814e+13, 1.739936e+00, 5.334214e+01, 4.014245e+00, 3.168656e+00)
upper_v2_r10 <- c(          38,          380,        0.325,           30,       3700.0,      10400.0,       6600.0,       6.2e13,           10,           58,           10,        15.00)

# Run the fitter
pheno_fit_v2_r10 <- phenologyFitter(par.guess = par_v2_r10,
                                    modelfn = PhenoFlex_GDHwrapper,
                                    bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                    SeasonList = season_list_v2,
                                    lower = lower_v2_r10,
                                    upper = upper_v2_r10,
                                    control = list(smooth = FALSE,
                                                   verbose = FALSE,
                                                   maxit = 10000,
                                                   nb.stop.improvement = 500))

# Some intermediate results
# Generate a data set to collect the outputs of the fitting for the calibration data 
out_df_v1_r10 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2_r10 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1_r10[["Predicted"]] <- pheno_fit_v1_r10$pbloomJDays
out_df_v2_r10[["Predicted"]] <- pheno_fit_v2_r10$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1_r10[["Error"]] <- pheno_fit_v1_r10$bloomJDays - pheno_fit_v1_r10$pbloomJDays
out_df_v2_r10[["Error"]] <- pheno_fit_v2_r10$bloomJDays - pheno_fit_v2_r10$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1_r10 <- RMSEP(out_df_v1_r10$Predicted, out_df_v1_r10$pheno)
RMSEP_calib_v2_r10 <- RMSEP(out_df_v2_r10$Predicted, out_df_v2_r10$pheno)

# Compute the RPIQ (Ratio of Performance to InterQuartile range) to account for non-normally distributed samples
RPIQ_calib_v1_r10 <- RPIQ(out_df_v1_r10$Predicted, out_df_v1_r10$pheno)
RPIQ_calib_v2_r10 <- RPIQ(out_df_v2_r10$Predicted, out_df_v2_r10$pheno)

# Estimate the mean error
mean(out_df_v1_r10$Error)
mean(out_df_v2_r10$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1_r10$Error))
mean(abs(out_df_v2_r10$Error))

# Compute the AICc
aic_fit_v1_r10 <- (2 * length(par_v1_r10)) + (length(calibration_seasons_v1) * 
                                                log(pheno_fit_v1_r10$model_fit$value / length(calibration_seasons_v1))) + ((2 * length(par_v1_r10)^2 + 2 * length(par_v1_r10)) /
                                                                                                                             length(calibration_seasons_v1) - length(par_v1_r10) - 1)
aic_fit_v2_r10 <- (2 * length(par_v2_r10)) + (length(calibration_seasons_v2) * 
                                                log(pheno_fit_v2_r10$model_fit$value / length(calibration_seasons_v2))) + ((2 * length(par_v2_r10)^2 + 2 * length(par_v2_r10)) /
                                                                                                                             length(calibration_seasons_v2) - length(par_v2_r10) - 1)

# Plot the observed versus predicted values
ggplot(out_df_v1_r10, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r10, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

# Finishing the fitting replications

# After 10 extensive fitting procedures (see maxit and nb.stop.improvement parameters in fitting 10), we observed the
# smallest RMSEP after run 9. We therefore, will proceed using the set of parameters fitted by the model in that run




# Validation ####
# Compute the bloom dates for the validation datasets
# Generate a validation data set with phenology data
valid_df_v1 <- pheno_v1[pheno_v1$Year %in% validation_seasons, ]
valid_df_v2 <- pheno_v2[pheno_v2$Year %in% validation_seasons, ]

# Add a third validation which is predicting the bloom dates in the five marginal seasons using the parameters fitted
# when calibrating the model WHITOUT these seasons. This is to test how capable the version excluding the marginal seasons is
# to predict phenology in marginal environments
valid_df_v2_marg_seasons <- pheno_v1[which(!(pheno_v1$Year %in% pheno_v2$Year)), ]

# Generate a list of seasons with weather data for the validation procedure
valid_season_list <- genSeasonList(data, mrange = c(9, 7), years = validation_seasons)
valid_season_list_marg_seasons <- genSeasonList(data, mrange = c(9, 7), years = valid_df_v2_marg_seasons$Year)

# Estimate the bloom dates with PhenoFlexGDHwrapper
for (i in 1 : nrow(valid_df_v1)) {
  
  valid_df_v1[i, "Predicted"] <- PhenoFlex_GDHwrapper(valid_season_list[[i]], pheno_fit_v1_r9$par)
}

# The same for the second version
for (i in 1 : nrow(valid_df_v2)) {
  
  valid_df_v2[i, "Predicted"] <- PhenoFlex_GDHwrapper(valid_season_list[[i]], pheno_fit_v2_r9$par)
}

# Same for the marginal seasons validation
for (i in 1 : nrow(valid_df_v2_marg_seasons)) {
  
  valid_df_v2_marg_seasons[i, "Predicted"] <- PhenoFlex_GDHwrapper(valid_season_list_marg_seasons[[i]],
                                                                   pheno_fit_v2_r9$par)
}

# Compute the error (observed - predicted)
valid_df_v1[["Error"]] <- valid_df_v1$pheno - valid_df_v1$Predicted
valid_df_v2[["Error"]] <- valid_df_v2$pheno - valid_df_v2$Predicted
valid_df_v2_marg_seasons[["Error"]] <- valid_df_v2_marg_seasons$pheno - valid_df_v2_marg_seasons$Predicted

# Estimate the RMSEP
RMSEP_valid_v1 <- RMSEP(valid_df_v1$Predicted, valid_df_v1$pheno, na.rm = TRUE)
RMSEP_valid_v2 <- RMSEP(valid_df_v2$Predicted, valid_df_v2$pheno, na.rm = TRUE)
RMSEP_valid_v2_marg_seasons <- RMSEP(valid_df_v2_marg_seasons$Predicted, valid_df_v2_marg_seasons$pheno, na.rm = TRUE)

# Compute de RPIQ for the validation set
RPIQ_valid_v1 <- RPIQ(valid_df_v1$Predicted, valid_df_v1$pheno)
RPIQ_valid_v2 <- RPIQ(valid_df_v2$Predicted, valid_df_v2$pheno)
RPIQ_valid_v2_marg_seasons <- RPIQ(valid_df_v2_marg_seasons$Predicted, valid_df_v2_marg_seasons$pheno)

# Plot the calibrated and validated 
ggplot(out_df_v1_r10, aes(pheno, Predicted)) +
  geom_point() +
  geom_point(data = valid_df_v1, aes(pheno, Predicted), 
             color = "red") + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")


ggplot(out_df_v2_r10, aes(pheno, Predicted)) +
  geom_point(shape = 1, size = 2) +
  geom_point(data = valid_df_v2, aes(pheno, Predicted, color = "'Normal' seasons"),
             size = 2.5) +
  geom_point(data = valid_df_v2_marg_seasons, aes(pheno, Predicted, color = "Marginal seasons"),
             size = 2.5) +
  scale_color_manual(values = c("firebrick", "cadetblue")) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed",
       title = "Framework excluding the marginal seasons in calibration",
       color = "Validation approach") +
  theme_bw()

ggsave("figures/validation_approach_for_excl_marginal_seasons.png")




# Bootstrap ####
# Implement the bootstrap of residuals for estimating errors. In this process, I implemented the idea based on my own
# understanding after looking at the PhenoFlex vignette. This means that I would need to check the main concept with the
# authors of the original paper. In brief, I did bootstrap the residuals of the calibration fitting to estimate the
# parameters boot.R times. Then, I computed the bloom dates in the validation data set boot.R times for each season.
# Since the parameters differ a little, the bloom dates also differ, making possible the estimation of a standard deviation 
# across boot.R values.

# Bootstrap the calibration fitting to obtain boot.R sets of parameters
fit_res_boot_v1 <- bootstrap.phenologyFit(pheno_fit_v1_r9,
                                          boot.R = 10,
                                          lower = lower,
                                          upper = upper,
                                          control = list(smooth = FALSE,
                                                         verbose = TRUE,
                                                         maxit = 1000,
                                                         nb.stop.improvement = 10))

# Same as above using version 2
fit_res_boot_v2 <- bootstrap.phenologyFit(pheno_fit_v2_r9,
                                          boot.R = 10,
                                          lower = lower,
                                          upper = upper,
                                          control = list(smooth = FALSE,
                                                         verbose = TRUE,
                                                         maxit = 1000,
                                                         nb.stop.improvement = 10))

# Take some quick look at the outputs
summary(fit_res_boot_v1)
summary(fit_res_boot_v2)

# Plot the bootstrap element
plot(fit_res_boot_v1)
plot(fit_res_boot_v2)

# Implement a for loop to compute the bloom dates based on the different sets of parameters from the bootstrap element.
# This loop will add the date and boot.R column to the validation data frame
for (k in 1 : length(fit_res_boot_v1$res)) {
  
  par_k <- fit_res_boot_v1$res[[k]]$par
  
  name <- paste0("Pred_Boot_", k)
  
  for (i in 1 : nrow(valid_df_v1)) {
    
    valid_df_v1[i, name] <- PhenoFlex_GDHwrapper(valid_season_list[[i]], par_k)
  }
}

# Same for version 2
for (k in 1 : length(fit_res_boot_v2$res)) {
  
  par_k <- fit_res_boot_v2$res[[k]]$par
  
  name <- paste0("Pred_Boot_", k)
  
  for (i in 1 : nrow(valid_df_v2)) {
    
    valid_df_v2[i, name] <- PhenoFlex_GDHwrapper(valid_season_list[[i]], par_k)
  }
}

# Compute the sd across boot.R bloom dates
valid_df_v1 <- valid_df_v1 %>% pivot_longer(starts_with("Pred_Boot_")) %>% group_by(Year, pheno, Predicted, Error) %>% 
  
  summarise(SD_boot = sd(value))

# Same for version
valid_df_v2 <- valid_df_v2 %>% pivot_longer(starts_with("Pred_Boot_")) %>% group_by(Year, pheno, Predicted, Error) %>% 
  
  summarise(SD_boot = sd(value))

# Generate a data set that contains all data for easy-faceting
out_df <- bind_rows("PhenoFlex[all]" = out_df_v1_r9,
                    "PhenoFlex[excluded]" = out_df_v2_r9, .id = "Version")
valid_df <- bind_rows("PhenoFlex[all]" = valid_df_v1,
                      "PhenoFlex[excluded]" = valid_df_v2, .id = "Version")

# Create a data set that computes de RSMEP for each facet
RMSEP_text <- data.frame(pheno = 48,
                         Predicted = c(153, 148, 143, 138, 133, 153, 148, 143, 138, 133, 133),
                         Version = c(rep("PhenoFlex[all]", 5), rep("PhenoFlex[excluded]", 6)),
                         Dataset = c("Calibration", "Validation", "Calibration", "Validation",
                                     "Calibration", "Validation", "Calibration", "Validation",
                                     "Calibration", "Validation", "Validation"))

# Plot all results including the error for the validation dots
ggplot() +
  geom_abline(intercept = 0, slope = 1, alpha = 0.35) +
  geom_point(data = filter(out_df, !(Year %in% pheno_v1$Year[which(!(pheno_v1$Year %in% pheno_v2$Year))])),
             aes(pheno, Predicted, shape = "Calibration"), color = "slategrey") +
  geom_point(data = filter(out_df, Year %in% pheno_v1$Year[which(!(pheno_v1$Year %in% pheno_v2$Year))]),
             aes(pheno, Predicted, fill = "Marginal seasons"), shape = 2, color = "firebrick") +
  geom_pointrange(data = valid_df,
                  aes(pheno, Predicted, ymin = Predicted - SD_boot, ymax = Predicted + SD_boot, color = "Validation"), 
                  size = 0.3, fatten = 0.1) +
  geom_text(data = RMSEP_text,  aes(pheno, Predicted),
            label = c(bquote("RMSE"["calib"]*" : "*.(round(RMSEP_calib_v1_r9, 1))),
                      bquote("RMSE"["valid"]*" : "*.(round(RMSEP_valid_v1, 1))),
                      bquote("RPIQ"["calib"]*"   : "*.(round(RPIQ_calib_v1_r9, 1))),
                      bquote("RPIQ"["valid"]*"   : "*.(round(RPIQ_valid_v1, 1))),
                      bquote("AICc"["calib"]*"    : "*.(round(aic_fit_v1_r9, 1))),
                      bquote("RMSE"["calib"]*" : "*.(round(RMSEP_calib_v2_r9, 1))),
                      bquote("RMSE"["valid"]*" : "*.(round(RMSEP_valid_v2, 1))),
                      bquote("RPIQ"["calib"]*"   : "*.(round(RPIQ_calib_v2_r9, 1))),
                      bquote("RPIQ"["valid"]*"   : "*.(round(RPIQ_valid_v2, 1))),
                      bquote("AICc"["calib"]*"    : "*.(sprintf("%0.1f", round(aic_fit_v2_r9, 1)))),
                      expression("")),
            hjust = 0, size = 2.5, fontface = "italic") +
  scale_x_continuous(breaks = seq(50, 125, 25),
                     labels = function (x) format(dormancyR::JDay_to_date(x, 2021), "%b %d")) +
  scale_y_continuous(breaks = seq(60, 150, 30),
                     labels = function (x) format(dormancyR::JDay_to_date(x, 2021), "%b %d")) +
  scale_shape_manual(values = 1) +
  scale_color_manual(values = "deepskyblue4") +
  labs(x = "Observed bloom date",
       y = "Predicted bloom date",
       color = NULL,
       shape = NULL,
       fill = NULL) +
  facet_grid(. ~ factor(Version, levels = c("PhenoFlex[all]", "PhenoFlex[excluded]")),
             labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.spacing = unit(-0.75, "cm"),
        legend.position = c(0.83, 0.09),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 8))

# Save the final plot to folder
ggsave("figures/model_performance_final_d_new_ana.png", width = 12, height = 10, units = "cm", dpi = 600)



# Estimate chill and heat responses from the calibration fitting
# First, load some helper functions to organize the data.
source("code/00_helper_functions.R")

# Create a data set with theoretical temperatures and heat and chill responses
temp_response_v1 <- data.frame(Temp = seq(-5, 60, 0.1),
                               Chill_res = gen_bell(pheno_fit_v1_r9$par, seq(-5, 60, 0.1)),
                               Heat_res = GDH_response(pheno_fit_v1_r9$par, seq(-5, 60, 0.1)))

temp_response_v2 <- data.frame(Temp = seq(-5, 60, 0.1),
                               Chill_res = gen_bell(pheno_fit_v2_r9$par, seq(-5, 60, 0.1)),
                               Heat_res = GDH_response(pheno_fit_v2_r9$par, seq(-5, 60, 0.1)))

# Pivot longer to generate a panel plot
temp_response_v1 <- pivot_longer(temp_response_v1, -Temp, names_to = "Var", values_to = "Response")
temp_response_v2 <- pivot_longer(temp_response_v2, -Temp, names_to = "Var", values_to = "Response")

# Generate a single data set
temp_response <- bind_rows("PhenoFlex[all]" = temp_response_v1,
                           "PhenoFlex[excluded]" = temp_response_v2,
                           .id = "version")

# Implement the plot. Generate two plots and then merge them to overcome the issue produced by
# facet_grid (free scales working but not as intended in this particular case)

# Chill plot
chill_response_plot <- ggplot(filter(temp_response, Var == "Chill_res"), aes(Temp, Response)) +
  geom_line(size = 1, color = "blue4") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(limits = c(-5, 25),
                     labels = function (x) paste0(x, " °C")) +
  scale_color_manual(values = c("blue", "red")) +
  labs(y = "Arbitrary units") +
  facet_grid(factor(version, levels = c("PhenoFlex[all]",
                                        "PhenoFlex[excluded]")) ~
               factor(Var, labels = c("Chill response"))) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        strip.text.y = element_blank(),
        strip.background = element_blank())

# Heat plot
heat_response_plot <- ggplot(filter(temp_response, Var == "Heat_res"), aes(Temp, Response, color = Var)) +
  geom_line(size = 1, color = "firebrick") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(labels = function (x) paste0(x, " °C")) +
  scale_color_manual(values = c("blue", "red")) +
  labs(y = "Arbitrary units") +
  facet_grid(factor(version, levels = c("PhenoFlex[all]",
                                        "PhenoFlex[excluded]")) ~ 
               factor(Var, labels = c("Heat~response")), labeller = label_parsed) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 8),
        strip.background = element_blank())

# Use patchwork syntax and functionality to merge the plots
(chill_response_plot + heat_response_plot) + plot_annotation(caption = "Temperature") &
  theme(plot.caption = element_text(hjust = 0.5, vjust = 1, size = 11))

# Save the final plot to folder
ggsave("figures/temp_responses_final_new_ana.png", width = 12, height = 10, units = "cm", dpi = 600)


# Compute some metrics for model validation
# Median error
med_residuals_v1 <- median(valid_df_v1$Error)
med_residuals_v2 <- median(valid_df_v2$Error)

# Mean absolute error
mean_abs_error_v1 <- mean(abs(valid_df_v1$Error))
mean_abs_error_v2 <- mean(abs(valid_df_v2$Error))

# IQR
IQR(valid_df_v1$Error)
IQR(valid_df_v2$Error)

# Create a small data frame for adding the metrics to the text
metrics_text <- data.frame(Version = c("PhenoFlex[all]", "PhenoFlex[all]",
                                       "PhenoFlex[excluded]", "PhenoFlex[excluded]",
                                       "PhenoFlex[excluded]"),
                           y = c(9, 8.5, 9, 8.5, 8))

# Plot the residuals to test for model quality
ggplot(valid_df, aes(Version, Error, fill = Version)) +
  geom_hline(yintercept = 0, alpha = 0.45, linetype = 2) +
  geom_boxplot(width = 0.25, size = 0.2, outlier.size = 0.5) +
  scale_x_discrete(limits = c("PhenoFlex[all]",
                              "PhenoFlex[excluded]")) +
  geom_text(data = metrics_text, aes(Version, y),
            label = c(bquote("Median error"*"      : "*.(round(med_residuals_v1, 1))),
                      bquote("Mean abs. error"*" : "*.(round(mean_abs_error_v1, 1))),
                      bquote("Median error"*"      : "*.(round(med_residuals_v2, 1))),
                      bquote("Mean abs. error"*" : "*.(round(mean_abs_error_v2, 1))),
                      expression("")),
            size = 1.7, hjust = 0, nudge_x = -0.35) +
  scale_fill_manual(values = c("cadetblue", "firebrick"),
                    breaks = c("PhenoFlex[all]",
                               "PhenoFlex[excluded]"),
                    labels = c(bquote("PhenoFlex"["all"]),
                               bquote("PhenoFlex"["excluded"]))) +
  labs(x = NULL,
       y = "Residuals (days)",
       fill = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.63, 0.075),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6))

ggsave("figures/validation_errors_b_new_ana.png", dpi = 600, width = 5, height = 6, units = "cm")  

