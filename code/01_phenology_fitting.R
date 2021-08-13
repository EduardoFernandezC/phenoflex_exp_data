library(chillR)
library(tidyverse)
library(patchwork)

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

# Generate two data sets according to the version of the analysis
pheno_v1 <- pheno[!(pheno$Treatment %in% c(36, 3)), ]
pheno_v2 <- pheno[!(pheno$Treatment %in% c(36, 3, 61, 13, 46, 62, 42, 58)), ]

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
calibration_seasons_v1 <- sort(sample(pheno_v1$Year, round(nrow(pheno_v1) * 0.75), replace = FALSE))
validation_seasons_v1 <- sort(pheno_v1[!(pheno_v1$Year %in% calibration_seasons_v1), "Year"])

calibration_seasons_v2 <- sort(sample(pheno_v2$Year, round(nrow(pheno_v2) * 0.75), replace = FALSE))
validation_seasons_v2 <- sort(pheno_v2[!(pheno_v2$Year %in% calibration_seasons_v2), "Year"])


# Set the initial parameters (wide ranges)
par <-   c(40, 190, 0.5, 25, 3372.8,  9900.3, 6319.5, 5.939917e13,  4, 36,  4,  1.60)
upper <- c(80, 500, 1.0, 30, 4000.0, 10000.0, 7000.0,       6.e13, 10, 40, 10, 20.00)
lower <- c(20, 100, 0.1, 0 , 3000.0,  9000.0, 6000.0,       5.e13,  0,  0,  0,  0.05)


# Define the list of seasons
season_list_v1 <- genSeasonList(data, mrange = c(9, 7), years = calibration_seasons_v1)
season_list_v2 <- genSeasonList(data, mrange = c(9, 7), years = calibration_seasons_v2)

# Run the fitter
pheno_fit_v1 <- phenologyFitter(par.guess = par,
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
pheno_fit_v2 <- phenologyFitter(par.guess = par,
                                modelfn = PhenoFlex_GDHwrapper,
                                bloomJDays = pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, "pheno"],
                                SeasonList = season_list_v2,
                                lower = lower,
                                upper = upper,
                                control = list(smooth = FALSE,
                                               verbose = FALSE,
                                               maxit = 1000,
                                               nb.stop.improvement = 5))

# Generate a dataset to collect the outputs of the fitting for the calibration data 
out_df_v1 <- pheno_v1[pheno_v1$Year %in% calibration_seasons_v1, ]
out_df_v2 <- pheno_v2[pheno_v2$Year %in% calibration_seasons_v2, ]

# Add the predicted values based on the fitting procedure
out_df_v1[["Predicted"]] <- pheno_fit_v1$pbloomJDays
out_df_v2[["Predicted"]] <- pheno_fit_v2$pbloomJDays

# Compute the error (observed - predicted)
out_df_v1[["Error"]] <- pheno_fit_v1$bloomJDays - pheno_fit_v1$pbloomJDays
out_df_v2[["Error"]] <- pheno_fit_v2$bloomJDays - pheno_fit_v2$pbloomJDays

# Compute the RMSEP
RMSEP_calib_v1 <- RMSEP(out_df_v1$Predicted, out_df_v1$pheno)
RMSEP_calib_v2 <- RMSEP(out_df_v2$Predicted, out_df_v2$pheno)

# Estimate the mean error
mean(out_df_v1$Error)
mean(out_df_v2$Error)

# Estimate the mean absolute error
mean(abs(out_df_v1$Error))
mean(abs(out_df_v2$Error))

# Plot the observed versus predicted values
ggplot(out_df_v1, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2, aes(pheno, Predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")


# Compute the bloom dates for the validation datasets
# Generate a validation data set with phenology data
valid_df_v1 <- pheno_v1[pheno_v1$Year %in% validation_seasons_v1, ]
valid_df_v2 <- pheno_v2[pheno_v2$Year %in% validation_seasons_v2, ]

# Generate a list of seasons with weather data for the validation procedure
valid_season_list_v1 <- genSeasonList(data, mrange = c(9, 7), years = validation_seasons_v1)
valid_season_list_v2 <- genSeasonList(data, mrange = c(9, 7), years = validation_seasons_v2)

# Estimate the bloom dates with PhenoFlexGDHwrapper
for (i in 1 : nrow(valid_df_v1)) {
  
  valid_df_v1[i, "Predicted"] <- PhenoFlex_GDHwrapper(valid_season_list_v1[[i]], pheno_fit_v1$par)
}

# The same for the second version
for (i in 1 : nrow(valid_df_v2)) {
  
  valid_df_v2[i, "Predicted"] <- PhenoFlex_GDHwrapper(valid_season_list_v2[[i]], pheno_fit_v2$par)
}

# Compute the error (observed - predicted)
valid_df_v1[["Error"]] <- valid_df_v1$pheno - valid_df_v1$Predicted
valid_df_v2[["Error"]] <- valid_df_v2$pheno - valid_df_v2$Predicted

# Estimate the RMSEP
RMSEP_valid_v1 <- RMSEP(valid_df_v1$pheno, valid_df_v1$Predicted)
RMSEP_valid_v2 <- RMSEP(valid_df_v2$pheno, valid_df_v2$Predicted)

# Plot the calibrated and validated 
ggplot(out_df_v1, aes(pheno, Predicted)) +
  geom_point() +
  geom_point(data = valid_df_v1, aes(pheno, Predicted), 
             color = "red") + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed")

ggplot(out_df_v2, aes(pheno, Predicted)) +
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
fit_res_boot_v1 <- bootstrap.phenologyFit(pheno_fit_v1,
                                          boot.R = 10,
                                          lower = lower,
                                          upper = upper,
                                          control = list(smooth = FALSE,
                                                         verbose = TRUE,
                                                         maxit = 10,
                                                         nb.stop.improvement = 5))

# Same as above using version 2
fit_res_boot_v2 <- bootstrap.phenologyFit(pheno_fit_v2,
                                          boot.R = 10,
                                          lower = lower,
                                          upper = upper,
                                          control = list(smooth = FALSE,
                                                         verbose = TRUE,
                                                         maxit = 10,
                                                         nb.stop.improvement = 5))

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
  
  valid_df_v1[i, name] <- PhenoFlex_GDHwrapper(valid_season_list_v1[[i]], par_k)
  }
}

# Same for version 2
for (k in 1 : length(fit_res_boot_v2$res)) {
  
  par_k <- fit_res_boot_v2$res[[k]]$par
  
  name <- paste0("Pred_Boot_", k)
  
  for (i in 1 : nrow(valid_df_v2)) {
    
    valid_df_v2[i, name] <- PhenoFlex_GDHwrapper(valid_season_list_v2[[i]], par_k)
  }
}

# Compute the sd across boot.R bloom dates
valid_df_v1 <- valid_df_v1 %>% pivot_longer(starts_with("Pred_Boot_")) %>% group_by(Year, pheno, Predicted, Error) %>% 
  
  summarise(SD_boot = sd(value))

# Same for version
valid_df_v2 <- valid_df_v2 %>% pivot_longer(starts_with("Pred_Boot_")) %>% group_by(Year, pheno, Predicted, Error) %>% 
  
  summarise(SD_boot = sd(value))

# Generate a data set that contains all data for easy-faceting
out_df <- bind_rows("Version 1" = out_df_v1, "Version 2" = out_df_v2, .id = "Version")
valid_df <- bind_rows("Version 1" = valid_df_v1, "Version 2" = valid_df_v2, .id = "Version")

# Create a data set that computes de RSMEP for each facet
RMSEP_text <- data.frame(pheno = 48,
                         Predicted = c(153, 148, 153, 148, 148),
                         Version = c("Version 1", "Version 1", "Version 2", "Version 2", "Version 2"),
                         Dataset = c("Calibration", "Validation", "Calibration", "Validation", "Validation"))

# Plot all results including the error for the validation dots
ggplot() +
  geom_abline(intercept = 0, slope = 1, alpha = 0.35) +
  geom_point(data = out_df, aes(pheno, Predicted, shape = "Calibration"), color = "slategrey") +
  geom_pointrange(data = valid_df,
                  aes(pheno, Predicted, ymin = Predicted - SD_boot, ymax = Predicted + SD_boot, color = "Validation"), 
                  size = 0.3, fatten = 0.1) +
  geom_text(data = RMSEP_text,  aes(pheno, Predicted),
            label = c(bquote("RMSEP"["calib"]*": "*.(round(RMSEP_calib_v1, 2))),
                      bquote("RMSEP"["valid"]*": "*.(round(RMSEP_valid_v1, 2))),
                      bquote("RMSEP"["calib"]*": "*.(round(RMSEP_calib_v2, 2))),
                      bquote("RMSEP"["valid"]*": "*.(round(RMSEP_valid_v2, 2))),
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
       shape = NULL) +
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
ggsave("figures/model_performance.png", width = 12, height = 10, units = "cm", dpi = 600)



# Estimate chill and heat responses from the calibration fitting
# First, load some helper functions to organize the data.
source("code/00_helper_functions.R")

# Create a data set with theoretical temperatures and heat and chill responses
temp_response_v1 <- data.frame(Temp = seq(-10, 35, 0.1),
                               Chill_res = gen_bell(pheno_fit_v1$par, seq(-10, 35, 0.1)),
                               Heat_res = GDH_response(pheno_fit_v1$par, seq(-10, 35, 0.1)))

temp_response_v2 <- data.frame(Temp = seq(-5, 40, 0.1),
                               Chill_res = gen_bell(pheno_fit_v2$par, seq(-5, 40, 0.1)),
                               Heat_res = GDH_response(pheno_fit_v2$par, seq(-5, 40, 0.1)))
  
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
  scale_x_continuous(labels = function (x) paste0(x, "°C")) +
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
  theme(plot.caption = element_text(hjust = 0.5, size = 11))

# Save the final plot to folder
ggsave("figures/temp_responses.png", width = 12, height = 10, units = "cm", dpi = 600)


