library(ggplot2)
library(gridExtra)
library(gratia)
library(cowplot)
library(mgcv)
theme_set(theme_bw())
library(dplyr)
library(tidymv)
library(readxl)

# Load data
data_infants <- read_excel("GAM_MC_only_Infants.xlsx")
data_mothers <- read_excel("GAM_MC_only_Mothers.xlsx")



# Define the first plot
gam1 <- gam(Ct ~ s(days, bs = "cr"), data = data_infants )
plot_obj1 <- plot(ggeffects::ggpredict(gam1), facets = TRUE) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "white", alpha = 1) +
  coord_cartesian(xlim = c(0, 125), ylim = c(45, 20)) +  # Set the y-axis limits in reverse order
  labs(x = "Infant age (days)", y = "MC Ct") +
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
# ggsave(file = "GAM_infants_MC_only.png",
#        plot = plot1,
#        width = 8,  # Width in inches
#        height = 6,  # Height in inches
#        units = "in",  # Specify units as inches
#        dpi = 100)  # Adjust DPI as needed

# Define the second plot
gam2 <- gam(Ct ~ s(days, bs = "cr"),  data = data_mothers )
plot_obj2 <- plot(ggeffects::ggpredict(gam2), facets = TRUE) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "white", alpha = 1) +
  coord_cartesian(xlim = c(0, 125), ylim = c(45, 20)) +  # Set the y-axis limits in reverse order
  labs(x = "Infant age (days)", y = "MC Ct") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120)) +  # Set the x-axis ticks
  theme(axis.text = element_text(size = 30)) +
  theme(text = element_text(size = 30)) +
  ggtitle("")  # Set a blank title



plot_obj2_with_data <- plot_obj2 + geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "red", alpha = 0.2)


# Increase the thickness of x and y axes
plot_obj2_with_data <- plot_obj2_with_data +
  theme(axis.text = element_text(size = 30),
        axis.line = element_line(linewidth = 2),  # Adjust the line width here
        text = element_text(size = 30))+
  guides(color = "none")  # Remove the color legend

plot_obj2_with_data <- plot_obj2_with_data + geom_line(color = "red", size = 1.5)

# Overlay data on smooths for the first plot
plot_obj2_with_data <- plot_obj2_with_data + geom_point(data = data_mothers, aes(x = days, y = Ct), size = 8,alpha = 0.7,color = "#3399CC")

plot2 <- grid.arrange(plot_obj2_with_data+ theme(strip.text.x = element_blank()), nrow = 1)


# 
# # Save the plot as a PNG file
# ggsave(file = "GAM_mothers_MC_only.png",
#        plot = plot2,
#        width = 8,  # Width in inches
#        height = 6,  # Height in inches
#        units = "in",  # Specify units as inches
#        dpi = 100)  # Adjust DPI as needed


comp <- compare_smooths(gam1, gam2)
draw(comp)

summary(gam1)
summary(gam2)
