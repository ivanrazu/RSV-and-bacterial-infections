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


# Load data
data_infants <- read_excel("GAM_HI_only_infants.xlsx")
data_mothers <- read_excel("GAM_HI_only_mothers.xlsx")



# Define the first plot
gam1 <- gam(Ct ~ s(days, bs = "cr"), data = data_infants)
plot_obj1 <- plot(ggeffects::ggpredict(gam1), facets = TRUE) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "white", alpha = 1) +
  coord_cartesian(xlim = c(0, 125), ylim = c(45, 20)) +  # Set the y-axis limits in reverse order
  labs(x = "Infant age (days)", y = "HI Ct") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120)) +  # Set the x-axis ticks
  theme(axis.text = element_text(size = 30)) +
  theme(text = element_text(size = 30)) +
  ggtitle("")  # Set a blank title


plot_obj1_with_data <- plot_obj1 + geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "red", alpha = 0.2)

# Overlay data on smooths for the first plot
plot_obj1_with_data <- plot_obj1_with_data + geom_point(data = data_infants, aes(x = days, y = Ct), size = 8,alpha = 0.6)

# Increase the thickness of x and y axes
plot_obj1_with_data <- plot_obj1_with_data +
  theme(axis.text = element_text(size = 30),
        axis.line = element_line(linewidth = 2),  # Adjust the line width here
        text = element_text(size = 30))+
  guides(color = "none")  # Remove the color legend

plot_obj1_with_data <- plot_obj1_with_data + geom_line(color = "red", size = 1.5)


plot1 <- grid.arrange(plot_obj1_with_data+ theme(strip.text.x = element_blank()), nrow = 1)



# Save the plot as a PNG file
# ggsave(file = "GAM_infants_HI_only.png",
#        plot = plot1,
#        width = 8,  # Width in inches
#        height = 6,  # Height in inches
#        units = "in",  # Specify units as inches
#        dpi = 100)  # Adjust DPI as needed
# 

#########################################################
# Because without specifying k get over fitting issues
# I use AIC to choose from different k values

set.seed(123)  # Set a seed for reproducibility
folds2 <- createDataPartition(data_mothers$Ct, p = 0.8, list = FALSE)

k_values2 <- seq(3, 8)  # Range of k values from 3 to 12. Had to increase and lower k
aic_values2 <- numeric(length(k_values2))

for (i in seq_along(k_values2)) {
  cv_model2 <- gam(Ct ~ s(days, bs = "cr", k = k_values2[i]), data = data_mothers[folds2, ])
  aic_values2[i] <- AIC(cv_model2)
}

best_k2 <- k_values2[which.min(aic_values2)]


#########################################################

gam2 <- gam(Ct ~ s(days, bs = "cr",k=5), data = data_mothers )
plot_obj2 <- plot(ggeffects::ggpredict(gam2), facets = TRUE) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "white", alpha = 1) +
  coord_cartesian(xlim = c(0, 125), ylim = c(45, 20)) +  # Set the y-axis limits in reverse order
  labs(x = "Infant age (days)", y = "HI Ct") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120)) +  # Set the x-axis ticks
  theme(axis.text = element_text(size = 30)) +
  theme(text = element_text(size = 30)) +
  ggtitle("")  # Set a blank title


plot_obj2_with_data <- plot_obj2 + geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "red", alpha = 0.2)

# Overlay data on smooths for the first plot
plot_obj2_with_data <- plot_obj2_with_data + geom_point(data = data_mothers, aes(x = days, y = Ct), size = 8,alpha = 0.7,color = "#3399CC")


# Increase the thickness of x and y axes
plot_obj2_with_data <- plot_obj2_with_data +
  theme(axis.text = element_text(size = 30),
        axis.line = element_line(linewidth = 2),  # Adjust the line width here
        text = element_text(size = 30))+
  guides(color = "none")  # Remove the color legend

plot_obj2_with_data <- plot_obj2_with_data + geom_line(color = "red", size = 1.5)


plot2 <- grid.arrange(plot_obj2_with_data+ theme(strip.text.x = element_blank()), nrow = 1)


# Save the plot as a PNG file
# ggsave(file = "GAM_mothers_HI_only.png",
#        plot = plot2,
#        width = 8,  # Width in inches
#        height = 6,  # Height in inches
#        units = "in",  # Specify units as inches
#        dpi = 100)  # Adjust DPI as needed
# 
# 

comp <- compare_smooths(gam1, gam2)
draw(comp)


summary(gam1)
summary(gam2)


