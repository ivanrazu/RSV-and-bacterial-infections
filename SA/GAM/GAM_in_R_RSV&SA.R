library(ggplot2)
library(gridExtra)
library(gratia)
library(cowplot)
library(mgcv)
theme_set(theme_bw())
library(dplyr)
library(tidymv)
library(readxl)

library(caret)


# There are no mothers in this group

# Load data
data_infants <- read_excel("GAM_RSV&SA_Infants.xlsx")



#########################################################
# Because without specifying k get over fitting issues
# I use AIC to choose from different k values

set.seed(123)  # Set a seed for reproducibility
folds <- createDataPartition(data_infants$Ct, p = 0.8, list = FALSE)

k_values <- seq(3, 5)  # Range of k values from 3 to 12. Had to increase and lower k
aic_values <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  cv_model <- gam(Ct ~ s(days, bs = "cr", k = k_values[i]), data = data_infants[folds, ])
  aic_values[i] <- AIC(cv_model)
}

best_k <- k_values[which.min(aic_values)]


#########################################################


# Define the first plot
gam1 <- gam(Ct ~ s(days, bs = "cr",k=4), data = data_infants )
plot_obj1 <- plot(ggeffects::ggpredict(gam1), facets = TRUE) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "white", alpha = 1) +
  coord_cartesian(xlim = c(0, 125), ylim = c(45, 20)) +  # Set the y-axis limits in reverse order
  labs(x = "Days after RSV+", y = "SA Ct") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120)) +  # Set the x-axis ticks
  theme(axis.text = element_text(size = 30)) +
  theme(text = element_text(size = 30)) +
  ggtitle("")  # Set a blank title


plot_obj1_with_data <- plot_obj1 + geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "red", alpha = 0.2)


# Increase the thickness of x and y axes
plot_obj1_with_data <- plot_obj1_with_data +
  theme(axis.text = element_text(size = 30),
        axis.line = element_line(linewidth = 2),  # Adjust the line width here
        text = element_text(size = 30))+
  guides(color = "none")  # Remove the color legend

plot_obj1_with_data <- plot_obj1_with_data + geom_line(color = "red", size = 1.5)

# Overlay data on smooths for the first plot
plot_obj1_with_data <- plot_obj1_with_data + geom_point(data = data_infants, aes(x = days, y = Ct), size = 8,alpha = 0.6)

plot1 <- grid.arrange(plot_obj1_with_data+ theme(strip.text.x = element_blank()), nrow = 1)



# # Save the plot as a PNG file
# ggsave(file = "GAM_infants_RSV&SA.png",
#        plot = plot1,
#        width = 8,  # Width in inches
#        height = 6,  # Height in inches
#        units = "in",  # Specify units as inches
#        dpi = 100)  # Adjust DPI as needed



summary(gam1)

