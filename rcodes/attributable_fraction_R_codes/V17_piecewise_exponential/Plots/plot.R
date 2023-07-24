mydf  <- data.frame(
  SN_types_N_cases = c(rep("Any SN (605)", 6), rep("SMN (454)", 6), rep("NMSC (249)", 6), rep("Breast cancer (76)", 6)),
  Variables = c(rep(c("Radiation", "Chemo", "All treatments", "PRS", "Lifestyle", "Combined"), 4)),
  Overall = c(0.432, 0.04, 0.454, 0.118, NA, 0.518, 
              0.383, 0.003, 0.384, 0.115, NA, 0.455,
              0.392, NA, 0.392, 0.428, NA, 0.652,
              0.509, 0.075, 0.584, 0.276, NA, 0.569),
  Female = c(0.426, 0.041, 0.448, 0.115, NA, 0.513,
             0.382, 0.003, 0.384, 0.113, NA, 0.454,
             0.384, NA, 0.384, 0.417, NA, 0.644,
             0.532, 0.07, 0.623, 0.288, NA, 0.623),
  Male = c(0.441, 0.04, 0.461, 0.12, NA, 0.526,
           0.383, 0.003, 0.383, 0.117, NA, 0.455,
           0.4, NA, 0.4, 0.44, NA, 0.662,
           0.486, 0.103, 0.554, 0.239, NA, 0.536),
  Age_less_than_35 = c(0.322, 0.072, 0.368, 0.118, NA, 0.443,
                       0.268, 0.012, 0.275, 0.113, NA, 0.357,
                       0.31, NA, 0.31, 0.421, NA, 0.601,
                       0.446, NA ,0.446 ,NA ,NA ,NA),
  Age_greater_than_or_equal_to_35 = c(0.456 ,0.034 ,0.472 ,0.118 ,NA ,0.534,
                                      0.406 ,0.001 ,0.406 ,0.115 ,NA ,0.474,
                                      0.399 ,NA ,0.399 ,0.428 ,NA ,0.657,
                                      0.513 ,NA ,0.494 ,NA ,NA ,NA)
)


# library(ggplot2)
# 
# # Create a new dataframe with the data in long format
# mydf_long <- tidyr::pivot_longer(mydf, cols = c("Overall", "Female", "Male", "Age_less_than_35", "Age_greater_than_or_equal_to_35"), 
#                                  names_to = "Gender_and_Age_Group", values_to = "Value")
# 
# # Create a stacked bar chart
# ggplot(mydf_long, aes(x = SN_types_N_cases, y = Value, fill = Gender_and_Age_Group)) +
#   geom_bar(stat = "identity") +
#   labs(title = "Comparison of Different Variables by SN Types",
#        x = "SN Types (N Cases)", y = "Proportion",
#        fill = "Gender and Age Group") +
#   theme_minimal() +
#   theme(legend.position = "bottom")


mydf_long <- gather(mydf, key = "AgeGroup", value = "AttributableFraction", Overall:Age_greater_than_or_equal_to_35)

# Plot
ggplot(mydf_long, aes(x = SN_types_N_cases, y = AttributableFraction, fill = Variables)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(title = "Attributable Fraction for Different Subsequent Neoplasms",
       x = "Subsequent Neoplasms (N cases)",
       y = "Attributable Fraction",
       fill = "Variables") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")




