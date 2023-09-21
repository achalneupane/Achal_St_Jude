# Load required libraries
library(ggplot2)
library(ggthemes)  # For additional themes

# Create a data frame with your data
data <- data.frame(
  SN_types = c("Any SN", "SMN", "NMSC", "Breast cancer", "Thyroid cancer", "Meningioma", "Sarcoma"),
  SJLIFE = c(0.518, 0.455, 0.652, 0.672, 0.862, 0.438, 0.371),
  CCSS = c(0.419, 0.326, 0.555, 0.709, 0.659, 0.432, 0.368)
)

# Reshape the data for ggplot2
library(reshape2)
data_melted <- melt(data, id.vars = "SN_types")

# Create the grouped bar chart
ggplot(data_melted, aes(x = SN_types, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  
  # Customize the theme and appearance
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjust font size and angle
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),  # Customize Y-axis title
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 12),  # Adjust legend font size
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 1.5, face = "bold"),  # Title formatting
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
        
  ) +
  
  # Customize colors
  scale_fill_manual(values = c("SJLIFE" = "#0072B2", "CCSS" = "#D55E00")) +  # Use color codes
  
  # Add labels
  labs(title = "Overall-AF without lifestyle", y = "Attributable Fraction", x = NULL) +
  
  # Remove extra space at the bottom of the Y-axis
  scale_y_continuous(expand = c(0, 0)) +
  
  # Add data labels
  geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(width = 0.8), vjust = -0.25, size = 3.5) +
  
  # Adjust legend position
  theme(legend.position="top", legend.box = "horizontal") 
  
  ## Add caption
  # labs(caption = "Data source: Your Source")

# Save the plot as a high-resolution image
# ggsave("grouped_bar_chart.png", width = 10, height = 6, dpi = 300)
