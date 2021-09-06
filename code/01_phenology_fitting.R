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
pheno_v1 <- pheno[!(pheno$Treatment %in% c(36, 3, 23, 24, 17, 18)), ]
pheno_v2 <- pheno[!(pheno$Treatment %in% c(36, 3, 61, 13, 46, 62, 42, 58, 23, 24, 17, 18)), ]

# Add the year column to match the data in the weather data frame
pheno_v1["Treatment"] <- pheno_v1$Treatment + 2019
pheno_v2["Treatment"] <- pheno_v2$Treatment + 2019

# Select only the relevant columns for further analysis
pheno_v1 <- pheno_v1[c("Treatment", "pheno")]
pheno_v2 <- pheno_v2[c("Treatment", "pheno")]

# Rename the columns for further compatibility
colnames(pheno_v1) <- c("Year", "pheno")
colnames(pheno_v2) <- c("Year", "pheno")


# Define a vector of calibration and validation seasons (20% of the seasons)
calibration_seasons <- sort(sample(pheno_v2$Year, round(nrow(pheno_v2) * 0.75), replace = FALSE))
calibration_seasons_v1 <- sort(c(calibration_seasons, pheno_v1$Year[which(!(pheno_v1$Year %in% pheno_v2$Year))]))
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
lower_v1_r2 <- c(          50,          300,            0,           20,       3000.0,       9000.0,       6000.0,        5.e13,            0,           26,            0,         0.05)
par_v1_r2   <- c(7.973532e+01, 3.515716e+02, 1.456674e-01, 2.550465e+01, 3.268874e+03, 9.897931e+03, 6.067611e+03, 5.939862e+13, 5.361414e+00, 3.067708e+01, 5.759079e+00, 8.604118e-01)
upper_v1_r2 <- c(          85,          400,          0.5,           30,       4000.0,      10000.0,       7000.0,        6.e13,           10,           36,           10,        30.00)

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
lower_v2_r2 <- c(          20,          300,            0,           18,       3000.0,       9000.0,       6000.0,        5.e13,            0,           30,            0,           10)
par_v2_r2   <- c(3.492167e+01, 3.872746e+02, 2.715607e-01, 2.176161e+01, 3.371032e+03, 9.900773e+03, 6.283549e+03, 5.939925e+13, 5.119899e+00, 3.999171e+01, 6.298089e+00, 3.999151e+01)
upper_v2_r2 <- c(          40,          400,          0.5,           30,       4000.0,      10000.0,       7000.0,        6.e13,           10,           45,           10,        55.00)

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
lower_v1_r3 <- c(          60,          320,            0,           20,       3000.0,       9000.0,       6000.0,        5.e13,            0,           25,            0,         0.05)
par_v1_r3   <- c(8.169113e+01, 3.402895e+02, 1.715215e-01, 2.572267e+01, 3.269949e+03, 9.897904e+03, 6.078208e+03, 5.939879e+13, 5.371916e+00, 2.928168e+01, 5.759231e+00, 1.290212e+00)
upper_v1_r3 <- c(          90,          370,          0.25,          30,       4000.0,      10000.0,       7000.0,        6.e13,           10,           35,           10,        30.00)

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
lower_v2_r3 <- c(          28,          330,            0,           20,       3000.0,       9000.0,       6000.0,        5.e13,            0,           30,            0,            5)
par_v2_r3   <- c(3.577411e+01, 3.528342e+02, 1.514894e-01, 2.207405e+01, 3.371040e+03, 9.900867e+03, 6.283443e+03, 5.939932e+13, 2.278396e+00, 4.498529e+01, 6.481487e+00, 4.256795e+01)
upper_v2_r3 <- c(          42,          380,         0.25,           30,       4000.0,      10000.0,       7000.0,        6.e13,           10,           45,           10,        45.00)

# Run the fitter
pheno_fit_v2_r3 <- phenologyFitter(par.guess = par_v2_r3,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r3,
                                   upper = upper_v2_r3,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 1000,
                                                  nb.stop.improvement = 5))

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
lower_v1_r4 <- c(          60,          340,            0,           22,       3000.0,       9000.0,       5500.0,        5.e13,            0,           24,            0,            1)
par_v1_r4   <- c(8.127247e+01, 3.621326e+02, 2.027210e-01, 2.786980e+01, 3.269950e+03, 9.897904e+03, 6.078365e+03, 5.939873e+13, 5.371916e+00, 2.795902e+01, 4.313384e+00, 1.290212e+00)
upper_v1_r4 <- c(          95,          380,          0.3,           32,       4000.0,      10500.0,       7000.0,       6.2e13,           10,           36,           10,        20.00)

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
lower_v2_r4 <- c(          25,          350,            0,           20,       3000.0,       9000.0,       5500.0,       5.0e13,            0,           35,            0,            5)
par_v2_r4   <- c(3.493029e+01, 3.748525e+02, 1.407212e-01, 2.330585e+01, 3.371041e+03, 9.899681e+03, 6.240408e+03, 5.939942e+13, 1.375911e+00, 4.498077e+01, 4.881018e+00, 3.311156e+01)
upper_v2_r4 <- c(          45,          390,         0.25,           30,       4000.0,      11000.0,       7000.0,       6.5e13,           10,           46,           10,        40.00)

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
lower_v1_r5 <- c(          70,          350,            0,           20,       3000.0,       9000.0,       5500.0,       5.5e13,            0,           20,            0,            1)
par_v1_r5   <- c(8.126363e+01, 3.744358e+02, 2.414810e-01, 2.788736e+01, 3.269947e+03, 9.897909e+03, 6.078255e+03, 5.939866e+13, 5.371916e+00, 2.806959e+01, 4.040326e+00, 1.290576e+00)
upper_v1_r5 <- c(          90,          390,         0.35,           35,       4000.0,      11000.0,       7000.0,       6.5e13,           10,           36,           10,        20.00)

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
lower_v2_r5 <- c(          27,          340,            0,           20,       3000.0,       9000.0,       5500.0,       5.5e13,            0,           35,            0,            5)
par_v2_r5   <- c(3.559138e+01, 3.556880e+02, 1.306868e-01, 2.340172e+01, 3.371042e+03, 9.899680e+03, 6.240416e+03, 5.939915e+13, 1.342226e+00, 4.566080e+01, 4.881018e+00, 1.184611e+01)
upper_v2_r5 <- c(          43,          370,         0.25,           30,       4000.0,      10500.0,       6500.0,       6.2e13,           10,           47,           10,        20.00)

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
lower_v1_r6 <- c(          75,          360,            0,           22,       3000.0,       9000.0,       5500.0,       5.5e13,            0,           23,            0,            1)
par_v1_r6   <- c(8.172607e+01, 3.808086e+02, 2.119380e-01, 2.772426e+01, 3.269964e+03, 9.897902e+03, 6.076986e+03, 5.939900e+13, 5.370861e+00, 2.773408e+01, 3.727458e+00, 1.295076e+00)
upper_v1_r6 <- c(          88,          390,        0.325,           32,       4000.0,      10500.0,       6500.0,       6.5e13,           10,           33,           10,        15.00)


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
                                                  nb.stop.improvement = 20))

# Same for version 2 (pheno_fit_v2_r5$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r6 <- c(          27,          325,            0,           20,       3000.0,       9500.0,       5500.0,       5.5e13,            0,           35,            0,            5)
par_v2_r6   <- c(3.715861e+01, 3.410498e+02, 8.301206e-02, 2.543932e+01, 3.371036e+03, 9.899691e+03, 6.240187e+03, 5.939902e+13, 1.322409e+00, 4.666922e+01, 3.368638e+00, 1.742363e+01)
upper_v2_r6 <- c(          47,          365,         0.15,           30,       4000.0,      10500.0,       6800.0,       6.2e13,           10,           47,           10,        25.00)

# Run the fitter
pheno_fit_v2_r6 <- phenologyFitter(par.guess = par_v2_r6,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r6,
                                   upper = upper_v2_r6,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 3000,
                                                  nb.stop.improvement = 50))

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
lower_v1_r7 <- c(          80,          375,            0,           24,       3000.0,       9000.0,       5500.0,       5.5e13,            0,           23,            0,            0)
par_v1_r7   <- c(8.406747e+01, 3.868252e+02, 2.802052e-01, 2.803642e+01, 3.269967e+03, 9.897790e+03, 6.123488e+03, 5.939797e+13, 5.543387e+00, 2.821836e+01, 3.435302e+00, 1.697243e+00)
upper_v1_r7 <- c(          88,          395,        0.325,           34,       4000.0,      10500.0,       6600.0,       6.5e13,           10,           33,           10,        10.00)

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
lower_v2_r7 <- c(          34,          335,            0,           20,       3000.0,       9400.0,       5800.0,       5.5e13,            0,           35,            0,            5)
par_v2_r7   <- c(3.715861e+01, 3.410498e+02, 8.301206e-02, 2.543932e+01, 3.371036e+03, 9.899691e+03, 6.240187e+03, 5.939902e+13, 1.322409e+00, 4.666922e+01, 3.368638e+00, 1.742363e+01)
upper_v2_r7 <- c(          40,          355,        0.125,           30,       4000.0,      10500.0,       6800.0,       6.4e13,           10,           48,           10,        25.00)

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
lower_v1_r8 <- c(          70,          380,            0,           24,       3000.0,       9300.0,       5600.0,       5.7e13,            0,           24,            0,            0)
par_v1_r8   <- c(8.069311e+01, 3.947530e+02, 3.011341e-01, 2.870726e+01, 3.270072e+03, 9.897215e+03, 6.066225e+03, 5.939720e+13, 5.830716e+00, 2.871877e+01, 3.013004e+00, 1.513857e+00)
upper_v1_r8 <- c(          90,          410,        0.375,           34,       4000.0,      10500.0,       6400.0,       6.2e13,           10,           32,           10,         8.00)

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
lower_v2_r8 <- c(          34,          335,            0,           22,       3000.0,       9500.0,       5800.0,       5.6e13,            0,           40,            0,            5)
par_v2_r8   <- c(3.703158e+01, 3.422739e+02, 8.370522e-02, 2.543982e+01, 3.371008e+03, 9.899714e+03, 6.239101e+03, 5.939914e+13, 1.323499e+00, 4.771689e+01, 3.387846e+00, 1.767365e+01)
upper_v2_r8 <- c(          40,          355,        0.125,           28,       4000.0,      10400.0,       6800.0,       6.3e13,           10,           48,           10,        25.00)

# Run the fitter
pheno_fit_v2_r8 <- phenologyFitter(par.guess = par_v2_r8,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r8,
                                   upper = upper_v2_r8,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 3000,
                                                  nb.stop.improvement = 50))

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
lower_v1_r9 <- c(          70,          380,            0,           22,       2800.0,       9400.0,       5500.0,       5.6e13,            0,           23,            0,            0)
par_v1_r9   <- c(8.110983e+01, 4.091862e+02, 3.495934e-01, 2.808232e+01, 3.270080e+03, 9.897752e+03, 6.067306e+03, 5.939717e+13, 5.457257e+00, 2.811405e+01, 3.068110e+00, 1.356180e+00)
upper_v1_r9 <- c(          90,          430,          0.4,           36,       4200.0,      10600.0,       6500.0,       6.3e13,           10,           36,           10,        10.00)

# Run the fitter
pheno_fit_v1_r9 <- phenologyFitter(par.guess = par_v1_r9,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r9,
                                   upper = upper_v1_r9,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 4000,
                                                  nb.stop.improvement = 200))

# Same for version 2 (pheno_fit_v2_r8$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r9 <- c(          33,          335,            0,           22,       3100.0,       9500.0,       5900.0,       5.6e13,            0,           40,            0,            5)
par_v2_r9   <- c(3.600597e+01, 3.433006e+02, 9.626594e-02, 2.617425e+01, 3.371018e+03, 9.899718e+03, 6.239022e+03, 5.939945e+13, 1.373597e+00, 4.696788e+01, 3.387809e+00, 1.553505e+01)
upper_v2_r9 <- c(          39,          355,        0.125,           30,       3800.0,      10500.0,       6500.0,       6.3e13,           10,           48,           10,        22.00)

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
lower_v1_r10 <- c(          78,          415,          0.2,           25,       3000.0,       9600.0,       5800.0,       5.7e13,            0,           25,            0,            0)
par_v1_r10   <- c(8.082702e+01, 4.291137e+02, 3.759529e-01, 2.818768e+01, 3.270077e+03, 9.897764e+03, 6.067065e+03, 5.939721e+13, 5.456977e+00, 2.832324e+01, 2.578843e+00, 1.355878e+00)
upper_v1_r10 <- c(          82,          435,          0.4,           31,       3500.0,      10200.0,       6200.0,       6.2e13,           10,           31,           10,        10.00)

# Run the fitter
pheno_fit_v1_r10 <- phenologyFitter(par.guess = par_v1_r10,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, "pheno"],
                                   SeasonList = season_list_v1,
                                   lower = lower_v1_r10,
                                   upper = upper_v1_r10,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 6000,
                                                  nb.stop.improvement = 300))

# Same for version 2 (pheno_fit_v2_r9$model_fit$par)
#                          yc,           zc,           s1,           Tu,           E0,           E1,           A0,           A1,           Tf,           Tc,           Tb,        slope
lower_v2_r10 <- c(          33,          330,            0,           22,       3100.0,       9500.0,       5900.0,       5.6e13,            0,           42,            0,            8)
par_v2_r10   <- c(3.575275e+01, 3.362357e+02, 9.665035e-02, 2.653338e+01, 3.371009e+03, 9.899714e+03, 6.238888e+03, 5.939968e+13, 1.357235e+00, 4.489110e+01, 3.387779e+00, 1.340575e+01)
upper_v2_r10 <- c(          39,          340,        0.125,           30,       3800.0,      10500.0,       6500.0,       6.3e13,            5,           46,            6,        16.00)

# Run the fitter
pheno_fit_v2_r10 <- phenologyFitter(par.guess = par_v2_r10,
                                   modelfn = PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                   SeasonList = season_list_v2,
                                   lower = lower_v2_r10,
                                   upper = upper_v2_r10,
                                   control = list(smooth = FALSE,
                                                  verbose = FALSE,
                                                  maxit = 7000,
                                                  nb.stop.improvement = 400))

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



# Compute the bloom dates for the validation datasets
# Generate a validation data set with phenology data
valid_df_v1 <- pheno_v1[pheno_v1$Year %in% validation_seasons, ]
valid_df_v2 <- pheno_v2[pheno_v2$Year %in% validation_seasons, ]

# Generate a list of seasons with weather data for the validation procedure
valid_season_list <- genSeasonList(data, mrange = c(9, 7), years = validation_seasons)

# Estimate the bloom dates with PhenoFlexGDHwrapper
for (i in 1 : nrow(valid_df_v1)) {
  
  valid_df_v1[i, "Predicted"] <- PhenoFlex_GDHwrapper(valid_season_list[[i]], pheno_fit_v1_r10$par)
}

# The same for the second version
for (i in 1 : nrow(valid_df_v2)) {
  
  valid_df_v2[i, "Predicted"] <- PhenoFlex_GDHwrapper(valid_season_list[[i]], pheno_fit_v2_r9$par)
}

# Compute the error (observed - predicted)
valid_df_v1[["Error"]] <- valid_df_v1$pheno - valid_df_v1$Predicted
valid_df_v2[["Error"]] <- valid_df_v2$pheno - valid_df_v2$Predicted

# Estimate the RMSEP
RMSEP_valid_v1 <- RMSEP(valid_df_v1$Predicted, valid_df_v1$pheno, na.rm = TRUE)
RMSEP_valid_v2 <- RMSEP(valid_df_v2$Predicted, valid_df_v2$pheno, na.rm = TRUE)

# Compute de RPIQ for the validation set
RPIQ_valid_v1 <- RPIQ(valid_df_v1$Predicted, valid_df_v1$pheno)
RPIQ_valid_v2 <- RPIQ(valid_df_v2$Predicted, valid_df_v2$pheno)

# Plot the calibrated and validated 
ggplot(out_df_v1_r10, aes(pheno, Predicted)) +
  geom_point() +
  geom_point(data = valid_df_v1, aes(pheno, Predicted), 
             color = "red") + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2_r9, aes(pheno, Predicted)) +
  geom_point() +
  geom_point(data = valid_df_v2, aes(pheno, Predicted), 
             color = "red") + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

# Implement the bootstrap of residuals for estimating errors. In this process, I implemented the idea based on my own
# understanding after looking at the PhenoFlex vignette. This means that I would need to check the main concept with the
# authors of the original paper. In brief, I did bootstrap the residuals of the calibration fitting to estimate the
# parameters boot.R times. Then, I computed the bloom dates in the validation data set boot.R times for each season.
# Since the parameters differ a little, the bloom dates also differ, making possible the estimation of a standard deviation 
# across boot.R values.

# Bootstrap the calibration fitting to obtain boot.R sets of parameters
fit_res_boot_v1 <- bootstrap.phenologyFit(pheno_fit_v1_r10,
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
out_df <- bind_rows("Version 1" = out_df_v1_r10, "Version 2" = out_df_v2_r10, .id = "Version")
valid_df <- bind_rows("Version 1" = valid_df_v1, "Version 2" = valid_df_v2, .id = "Version")

# Create a data set that computes de RSMEP for each facet
RMSEP_text <- data.frame(pheno = 48,
                         Predicted = c(153, 148, 143, 138, 133, 153, 148, 143, 138, 133, 133),
                         Version = c("Version 1", "Version 1", "Version 1", "Version 1", "Version 1",
                                     "Version 2", "Version 2", "Version 2", "Version 2", "Version 2", "Version 2"),
                         Dataset = c("Calibration", "Validation", "Calibration", "Validation",
                                     "Calibration", "Validation", "Calibration", "Validation",
                                     "Calibration", "Validation", "Validation"))

# Plot all results including the error for the validation dots
ggplot() +
  geom_abline(intercept = 0, slope = 1, alpha = 0.35) +
  geom_point(data = out_df, aes(pheno, Predicted, shape = "Calibration"), color = "slategrey") +
  geom_point(data = filter(out_df, Year %in% pheno_v1$Year[which(!(pheno_v1$Year %in% pheno_v2$Year))]),
             aes(pheno, Predicted, fill = "Limiting seasons"), shape = 1, color = "firebrick4") +
  geom_pointrange(data = valid_df,
                  aes(pheno, Predicted, ymin = Predicted - SD_boot, ymax = Predicted + SD_boot, color = "Validation"), 
                  size = 0.3, fatten = 0.1) +
  geom_text(data = RMSEP_text,  aes(pheno, Predicted),
            label = c(bquote("RMSE"["calib"]*" : "*.(round(RMSEP_calib_v1_r10, 1))),
                      bquote("RMSE"["valid"]*" : "*.(round(RMSEP_valid_v1, 1))),
                      bquote("RPIQ"["calib"]*"   : "*.(round(RPIQ_calib_v1_r10, 1))),
                      bquote("RPIQ"["valid"]*"   : "*.(round(RPIQ_valid_v1, 1))),
                      bquote("AICc"["calib"]*"    : "*.(round(aic_fit_v1_r10, 1))),
                      bquote("RMSE"["calib"]*" : "*.(round(RMSEP_calib_v2_r10, 1))),
                      bquote("RMSE"["valid"]*" : "*.(round(RMSEP_valid_v2, 1))),
                      bquote("RPIQ"["calib"]*"   : "*.(round(RPIQ_calib_v2_r9, 1))),
                      bquote("RPIQ"["valid"]*"   : "*.(round(RPIQ_valid_v2, 1))),
                      bquote("AICc"["calib"]*"    : "*.(round(aic_fit_v2_r10, 1))),
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
  facet_grid(. ~ Version) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.spacing = unit(-0.75, "cm"),
        legend.position = c(0.85, 0.09),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 8))

# Save the final plot to folder
ggsave("figures/model_performance_final_b.png", width = 12, height = 10, units = "cm", dpi = 600)



# Estimate chill and heat responses from the calibration fitting
# First, load some helper functions to organize the data.
source("code/00_helper_functions.R")

# Create a data set with theoretical temperatures and heat and chill responses
temp_response_v1 <- data.frame(Temp = seq(-3, 48, 0.1),
                               Chill_res = gen_bell(pheno_fit_v1_r10$par, seq(-3, 48, 0.1)),
                               Heat_res = GDH_response(pheno_fit_v1_r10$par, seq(-3, 48, 0.1)))

temp_response_v2 <- data.frame(Temp = seq(-3, 48, 0.1),
                               Chill_res = gen_bell(pheno_fit_v2_r9$par, seq(-3, 48, 0.1)),
                               Heat_res = GDH_response(pheno_fit_v2_r9$par, seq(-3, 48, 0.1)))
  
# Pivot longer to generate a panel plot
temp_response_v1 <- pivot_longer(temp_response_v1, -Temp, names_to = "Var", values_to = "Response")
temp_response_v2 <- pivot_longer(temp_response_v2, -Temp, names_to = "Var", values_to = "Response")

# Generate a single data set
temp_response <- bind_rows("version 1" = temp_response_v1,
                           "version 2" = temp_response_v2,
                           .id = "version")

# Implement the plot. Generate two plots and then merge them to overcome the issue produced by
# facet_grid (free scales working but not as intended in this particular case)

# Chill plot
chill_response_plot <- ggplot(filter(temp_response, Var == "Chill_res"), aes(Temp, Response)) +
  geom_line(size = 1, color = "blue4") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(limits = c(-3, 27),
                     labels = function (x) paste0(x, "°C")) +
  scale_color_manual(values = c("blue", "red")) +
  labs(y = "Arbitrary units") +
  facet_grid(version ~ factor(Var, labels = c("Chill response"))) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        strip.text.y = element_blank(),
        strip.background = element_blank())

# Heat plot
heat_response_plot <- ggplot(filter(temp_response, Var == "Heat_res"), aes(Temp, Response, color = Var)) +
  geom_line(size = 1, color = "firebrick") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(labels = function (x) paste0(x, "°C")) +
  scale_color_manual(values = c("blue", "red")) +
  labs(y = "Arbitrary units") +
  facet_grid(factor(version, labels = c("Version 1", "Version 2")) ~ factor(Var, labels = c("Heat response"))) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 8),
        strip.background = element_blank())

# Use patchwork syntax and functionality to merge the plots
(chill_response_plot + heat_response_plot) + plot_annotation(caption = "Temperature") &
  theme(plot.caption = element_text(hjust = 0.5, vjust = 1, size = 11))

# Save the final plot to folder
ggsave("figures/temp_responses_final.png", width = 12, height = 10, units = "cm", dpi = 600)


