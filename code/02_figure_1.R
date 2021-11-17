library(tidyverse)
library(chillR)
library(cowplot)

# This script reproduces Figure 1 of the main manuscript

# Import the weather data from the folder
data <- read.csv("data/final_weather_data_S1_S2_apple_hourly.csv")

# Generate a new column (Treat_season) to identify unique treatment ids
data[["Treat_season"]] <- ifelse(data$Treatment <= 33, paste0(data$Treatment, "_S1"),
                                 paste0(data$Treatment - 33, "_S2"))

# Add the date
data[["Date"]] <- as.Date(format(YEARMODA2Date(data$YEARMODA), format = "%Y-%m-%d")) 

# Substract 1 year to the second season for plotting compatibility
data[data$Treatment >= 34, "Date"] <- data[data$Treatment >= 34, "Date"] - 365


# Extract the last dates to remove the data according to the treatment
# Load the phenology data
# Import the phenology data from the repository
pheno <- read.csv("data/final_bio_data_S1_S2_apple.csv")

# Transform the JDay to date
pheno[["Date"]] <- dormancyR::JDay_to_date(pheno$pheno, year = 2019, na.rm = TRUE)

# Add the Treat_season column for further compatibility
pheno[["Treat_season"]] <- ifelse(pheno$Treatment <= 33, paste0(pheno$Treatment, "_S1"),
                                  paste0(pheno$Treatment - 33, "_S2"))


# Nest the weather dataframe to delete the rows after blooming
nested <- data %>% group_by(Treat_season) %>% nest()

# Merge the columns in the last_date and the nested data frame
nested[["last_date"]] <- pheno$Date

# Unnest the nested data frame to filter the rows
unnest <- unnest(nested, cols = c(data))

# Set 0 after the last date
unnest[["Temp_2"]] <- ifelse(unnest$Date > unnest$last_date + 5, as.numeric(NA), unnest$Temp)

# Define the order of the treatments according to temperature between the begining of the exp and the date of bloom
temp_treatments <- unnest %>% 
  filter(Date %in% c(as.Date("2018-10-01") : as.Date("2019-05-31"))) %>% 
  group_by(Treat_season) %>% summarise(mean_temp = mean(Temp_2, na.rm = TRUE),
                                       median_temp = median(Temp_2, na.rm = TRUE),
                                       sd_temp = sd(Temp_2, na.rm = TRUE)) %>% 
  filter(!(Treat_season %in% c("3_S1", "3_S2", "34_S2", "17_S1",
                               "18_S1", "23_S1", "24_S1", "28_S2")))

# Extract the ordered treatments
ordered_treatments_temp <- temp_treatments[order(temp_treatments$mean_temp, decreasing = TRUE),
                                           "Treat_season"][["Treat_season"]]

# Merge the date of bloom with the temperature dataframe for plotting reasons
temp_treat_labels <- left_join(temp_treatments, pheno, by = "Treat_season")

# Filter the data
unnest <- filter(unnest, !(Treat_season %in% c("3_S1", "3_S2", "34_S2", "17_S1",
                                               "18_S1", "23_S1", "24_S1", "28_S2")))

pheno <- filter(pheno, !(Treat_season %in% c("3_S1", "3_S2", "34_S2", "17_S1",
                                             "18_S1", "23_S1", "24_S1", "28_S2")))

# Create the plot
# Marginal seasons only
marginal_seasons <- ggplot() +
  geom_tile(data = filter(unnest, Treat_season %in% c("9_S2", "13_S1", "29_S2",
                                                      "25_S2", "13_S2")),
            aes(Date, Treat_season, fill = runn_mean(Temp_2, 9 * 24)), height = 0.8) +
  geom_point(data = filter(pheno, Treat_season %in% c("9_S2", "13_S1", "29_S2",
                                                      "25_S2", "13_S2")),
             aes(Date, Treat_season, shape = "Bloom date"), color = "navyblue", size = 2) +
  geom_text(aes(Date + 10, Treat_season,
                label = paste(sprintf("%0.2f", round(mean_temp, 2)), "°C")),
            data = filter(temp_treat_labels, Treat_season %in% c("9_S2", "13_S1", "29_S2",
                                                                 "25_S2", "13_S2")),
            hjust = 0, size = 2.5, color = "skyblue4") +
  scale_fill_gradient2(low = "blue4", mid = "white", high = "red4", midpoint = 8,
                       labels = function (x) paste(x, "°C"), na.value = "White") +
  scale_shape_manual(values = 1) +
  scale_x_date(limits = c(as.Date("2018-10-01"), as.Date("2019-06-05")),
               expand = expansion(mult = c(0, 0.12))) +
  scale_y_discrete(limits = rev(ordered_treatments_temp[1 : 5]),
                   labels = c(5 : 1), position = "right") +
  labs(x = NULL,
       y = NULL,
       fill = NULL,
       shape = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 0.2,
        plot.margin = unit(c(1, 4.5, 1, 1), "cm"),
        axis.text.x = element_text(hjust = 0))

# Remaining seasons
normal_seasons <- ggplot() +
  geom_tile(data = filter(unnest, !(Treat_season %in% c("9_S2", "13_S1", "29_S2",
                                                        "25_S2", "13_S2"))),
            aes(Date, Treat_season, fill = runn_mean(Temp_2, 9 * 24)), height = 0.8) +
  geom_point(data = filter(pheno, !(Treat_season %in% c("9_S2", "13_S1", "29_S2",
                                                        "25_S2", "13_S2"))),
             aes(Date, Treat_season, shape = "Bloom date"), color = "navyblue", size = 2) +
  geom_text(aes(Date + 10, Treat_season,
                label = paste(sprintf("%0.2f", round(mean_temp, 2)), "°C")),
            data = filter(temp_treat_labels, !(Treat_season %in% c("9_S2", "13_S1", "29_S2",
                                                                   "25_S2", "13_S2"))),
            hjust = 0, size = 2.5, color = "skyblue4") +
  scale_fill_gradient2(low = "blue4", mid = "white", high = "red4", midpoint = 8,
                       labels = function (x) paste(x, "°C"), na.value = "White") +
  scale_shape_manual(values = 1) +
  scale_x_date(limits = c(as.Date("2018-10-01"), as.Date("2019-06-05")),
               expand = expansion(mult = c(0, 0.12))) +
  scale_y_discrete(limits = rev(ordered_treatments_temp[6 : 59]),
                   labels = c(59 : 6), position = "right") +
  labs(x = NULL,
       y = NULL,
       fill = "Temperature",
       shape = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 2,
        legend.box.spacing = unit(0.6, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.title = element_text(vjust = 1, size = 10, color = "grey40"),
        legend.text = element_text(size = 8.5, color = "grey40"),
        legend.spacing.y = unit(0.1, "cm"),
        axis.text.x = element_text(hjust = 0))

# Histogram
histogram <- ggplot(temp_treatments,
                    aes(mean_temp, fill = Treat_season %in% c("9_S2", "13_S1", "29_S2",
                                                              "25_S2", "13_S2"))) +
  geom_histogram(binwidth = 1) +
  scale_fill_manual(values = c("cadetblue", "firebrick")) +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 4, 8, 12, 16),
                     trans = "reverse",
                     expand = expansion(mult = c(0.05, 0))) +
  scale_x_continuous(labels = function (x) paste(x, "°C"), 
                     expand = expansion(mult = c(0.005, 0.005))) +
  labs(x = "Mean temperature",
       y = "Number of\n seasons") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 7,
        axis.text = element_text(size = 8.5),
        # axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 10, color = "grey20"),
        axis.title.y = element_text(size = 10, color = "grey20"))

# Merge the plot with ggdraw
ggdraw() +
  draw_plot(marginal_seasons, 0.143, 0.40995, 1, 1, 0.8325) +
  draw_plot(normal_seasons, 0.138, -0.0749, 1, 1, 0.75) +
  draw_plot(histogram, -0.36, -0.0105, 1, 1, 0.95) + 
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  draw_line(x = c(0.1, 0.8), y = 0.828,
            size = 0.5, linetype = "dashed", color = "grey60") +
  draw_plot_label(c("Marginal seasons", "Non-marginal seasons"),
                  x = 0.8, y = c(0.842, 0.81), hjust = 1, vjust = 0,
                  size = 8, fontface = "italic", color = c("firebrick", "cadetblue")) +
  draw_plot_label("Experimental season", x = 0.812, y = 0.59, size = 9.5, hjust = 0, vjust = 0,
                  fontface = "plain", color = "grey20", angle = 270) +
  draw_plot_label(c("A)", "B)", "C)"), x = c(0.09, 0.28, 0.28), y = c(0.985, 0.985, 0.82),
                  size = 10, fontface = "plain")

# Save the plot to the folder
ggsave("figures/treatments_d.png",
       device = "png", dpi = 600, width = 16.6, height = 22.4, units = "cm")
