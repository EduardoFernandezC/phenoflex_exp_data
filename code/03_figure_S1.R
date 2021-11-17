library(chillR)
library(tidyverse)

# Load the environment generated during the study
load("fitting_09_nov_final.RData")

# Generate a validation dataframe for supplementary purpose
valid_df_v2_marg_seasons <- pheno_v1[which(!(pheno_v1$Year %in% pheno_v2$Year)), ]

# Generate a list of seasons with weather data for the validation procedure
valid_season_list_marg_seasons <- genSeasonList(data, mrange = c(9, 7), years = valid_df_v2_marg_seasons$Year)

# Estimate the predicted bloom date
for (i in 1 : nrow(valid_df_v2_marg_seasons)) {
  
  valid_df_v2_marg_seasons[i, "Predicted"] <- PhenoFlex_GDHwrapper(valid_season_list_marg_seasons[[i]],
                                                                   pheno_fit_v2_r9$par)
}

# Compute the error (observed - predicted)
valid_df_v2_marg_seasons[["Error"]] <- valid_df_v2_marg_seasons$pheno - valid_df_v2_marg_seasons$Predicted

# Estimate the RMSEP
RMSEP_valid_v2_marg_seasons <- RMSEP(valid_df_v2_marg_seasons$Predicted, valid_df_v2_marg_seasons$pheno, na.rm = TRUE)

# Compute de RPIQ for the validation set
RPIQ_valid_v2_marg_seasons <- RPIQ(valid_df_v2_marg_seasons$Predicted, valid_df_v2_marg_seasons$pheno)


# Add the erro for bootstrap ####
# Validate k times the marginal seasons
for (k in 1 : length(fit_res_boot_v2$res)) {
  
  par_k <- fit_res_boot_v1$res[[k]]$par
  
  name <- paste0("Pred_Boot_", k)
  
  for (i in 1 : nrow(valid_df_v2_marg_seasons)) {
    
    valid_df_v2_marg_seasons[i, name] <- PhenoFlex_GDHwrapper(valid_season_list_marg_seasons[[i]], par_k)
  }
}

# Compute the sd across boot.R bloom dates
valid_df_v2_marg_seasons <- valid_df_v2_marg_seasons %>% 
  pivot_longer(starts_with("Pred_Boot_")) %>% group_by(Year, pheno, Predicted, Error) %>% 
  summarise(SD_boot = sd(value))

# Text dataframe for labels
RMSEP_text_sup <- data.frame(pheno = 121,
                             Predicted = c(53, 48, 43, 38, 33, 33))

# Plot the results nicely
ggplot() +
  geom_abline(intercept = 0, slope = 1, alpha = 0.35) +
  geom_point(data = filter(out_df, !(Year %in% pheno_v1$Year[which(!(pheno_v1$Year %in% pheno_v2$Year))]) &
                             Version == "PhenoFlex[excl.~marginal~seasons]"),
             aes(pheno, Predicted, shape = "Calibration"), color = "slategrey") +
  geom_pointrange(data = valid_df_v2_marg_seasons,
                  aes(pheno, Predicted, ymin = Predicted - SD_boot, ymax = Predicted + SD_boot, color = "Validation"), 
                  size = 0.3, fatten = 0.1) +
  geom_text(data = RMSEP_text_sup,  aes(pheno, Predicted),
            label = c(bquote("RMSE"["calib"]*" : "*.(round(RMSEP_calib_v2_r9, 1))),
                      bquote("RMSE"["valid"]*" : "*.(round(RMSEP_valid_v2_marg_seasons, 1))),
                      bquote("RPIQ"["calib"]*"   : "*.(round(RPIQ_calib_v2_r9, 1))),
                      bquote("RPIQ"["valid"]*"   : "*.(round(RPIQ_valid_v2_marg_seasons, 1))),
                      bquote("AICc"["calib"]*"    : "*.(sprintf("%0.1f", round(aic_fit_v2_r9, 1)))),
                      expression("")),
            hjust = 0, size = 2.5, fontface = "italic") +
  scale_x_continuous(breaks = seq(50, 125, 25),
                     labels = function (x) format(dormancyR::JDay_to_date(x, 2021), "%b %d")) +
  scale_y_continuous(breaks = seq(60, 150, 30),
                     labels = function (x) format(dormancyR::JDay_to_date(x, 2021), "%b %d")) +
  scale_shape_manual(values = 1) +
  scale_color_manual(values = "firebrick") +
  labs(x = "Observed bloom date",
       y = "Predicted bloom date",
       color = NULL,
       shape = NULL,
       fill = NULL) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.spacing = unit(-0.75, "cm"),
        legend.position = c(0.83, 0.92),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 8))
  
# Save the final plot to folder
ggsave("figures/model_performance_final_d_new_ana_supp.png", width = 12, height = 10, units = "cm", dpi = 600)


