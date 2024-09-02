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
data_infants <- read_excel("C:/Users/ivanr/Documents/SAMIPS Project/RSV_SA_summarized/GAM_SA_RSV_infants.xlsx")

data_mothers <- read_excel("C:/Users/ivanr/Documents/SAMIPS Project/RSV_SA_summarized/GAM_SA_RSV_mothers.xlsx")


# Define the first plot
gam1 <- gam(Ct ~ s(days, bs = "cr",k=3), data = data_infants )
plot_obj1 <- plot(ggeffects::ggpredict(gam1), facets = TRUE) +
  coord_cartesian(xlim = c(0, 125), ylim = c(20, 42)) +
  labs(x = "Infant age (days)", y = "MC Ct") +
  theme(axis.text = element_text(size = 30)) +
  theme(text = element_text(size = 30))


gam2 <- gam(Ct ~ s(days, bs = "cr"), data = data_mothers )
plot_obj2 <- plot(ggeffects::ggpredict(gam2), facets = TRUE) +
  coord_cartesian(xlim = c(0, 125), ylim = c(20, 42)) +
  labs(x = "Infant age (days)", y = "MC Ct") +
  theme(axis.text = element_text(size = 30)) +
  theme(text = element_text(size = 30))


# Overlay data on smooths for the first plot
plot_obj1_with_data <- plot_obj1 + geom_point(data = data_infants, aes(x = days, y = Ct), size = 8,alpha = 0.6)

plot_obj2_with_data <- plot_obj2 + geom_point(data = data_mothers, aes(x = days, y = Ct), size = 8,alpha = 0.7,color = "#3399CC")



# Make the lines representing the smooths thicker with the original red color for the first plot
plot_obj1_with_data <- plot_obj1_with_data + geom_line(aes(color = "red"), linewidth = 1)
plot_obj2_with_data <- plot_obj2_with_data + geom_line(aes(color = "red"), size = 1)



# Increase the thickness of x and y axes
plot_obj1_with_data <- plot_obj1_with_data +
  theme(axis.text = element_text(size = 30),
        axis.line = element_line(linewidth = 2),  # Adjust the line width here
        text = element_text(size = 30))



# Increase the thickness of x and y axes
plot_obj2_with_data <- plot_obj2_with_data +
  theme(axis.text = element_text(size = 30),
        axis.line = element_line(size = 2),  # Adjust the line width here
        text = element_text(size = 30))


# Remove the title "days" from the top of the plot
plot1 <- grid.arrange(plot_obj1_with_data + theme(strip.text.x = element_blank()), nrow = 1)

plot2 <- grid.arrange(plot_obj2_with_data + theme(strip.text.x = element_blank()), nrow = 1)


# Display the combined plot
print(combined_plot)

summary(gam1)
summary(gam2)
