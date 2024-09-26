# Load required packages
library(shiny)
library(ggplot2)

## SEX 1 is Male; 2 is Female

# Define UI for the CMP prediction app
ui <- fluidPage(
  titlePanel("CMP Risk Prediction Model"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("age_diagnosis", "Age at Primary Diagnosis (years):", value = 0, min = 0, max = 200, step = 1),
      numericInput("follow_up_years", "Follow-up years:", value = 0, min = 0, max = 200, step = 1),
      selectInput("sex", "Sex:", choices = list("Female" = "female", "Male" = "male")),
      selectInput("anthracycline", "Cumulative Anthracycline Dose (mg/m2):",
                  choices = list("None" = "none", ">0-100" = "0_100", ">100-250" = "100_250", ">250" = "250")),
      selectInput("radiation", "Mean Heart Radiation Dose (Gray):",
                  choices = list("None" = "none", ">5" = "5", ">5-15" = "5_15", ">15-35" = "15_35", ">35" = "35")),
      selectInput("age_baseline", "Age at Baseline (years):",
                  choices = list("≤25" = "25", ">25-35" = "25_35", ">35-45" = "35_45", ">45" = "45")),
      selectInput("hypertension", "Hypertension:", choices = list("No" = "no", "Yes" = "yes")),
      selectInput("ancestry", "Genetic Ancestry:",
                  choices = list("European" = "european", "African" = "african", "Others" = "others")),
      numericInput("pr_score1", "Polygenic Risk Score HCM:", value = 0.02, min = 0.009631459, max = 0.06178067, step = 0.001),
      numericInput("pr_score2", "Polygenic Risk Score LVEVi:", value = 0.05, min = 0.04326087, max = 0.06923913, step = 0.001)
    ),
    
    mainPanel(
      h3("Predicted CMP Risk"),
      textOutput("risk_output"),
      uiOutput("risk_level"),
      plotOutput("risk_plot")
    )
  )
)

# Define server logic for CMP risk prediction
server <- function(input, output) {
  
  # Function to categorize continuous age into ranges
  categorize_age_diagnosis <- function(age) {
    if (age <= 5) {
      return("5")
    } else if (age > 5 && age <= 10) {
      return("5_10")
    } else if (age > 10 && age <= 15) {
      return("10_15")
    } else {
      return("15")
    }
  }
  
  # Function to calculate risk based on inputs
  calculate_risk <- function() {
    # Coefficients from the provided table (log(RR))
    coeffs <- list(
      "age_diagnosis" = c("5" = 0, "5_10" = -0.0382, "10_15" = -0.1181, "15" = -0.2631),
      "sex" = c("female" = 0, "male" = 0.4457),
      "anthracycline" = c("none" = 0, "0_100" = -0.0157, "100_250" = 0.9204, "250" = 1.7179),
      "radiation" = c("none" = 0, "5" = -0.1228, "5_15" = 0.3428, "15_35" = 0.9745, "35" = 2.818),
      "age_baseline" = c("25" = 0, "25_35" = 0.1775, "35_45" = 0.6221, "45" = 0.9118),
      "hypertension" = c("no" = 0, "yes" = 0.8247),
      "ancestry" = c("european" = 0, "african" = 0.6572, "others" = -0.5715)
    )
    
    # Categorize the age diagnosis
    categorized_age <- categorize_age_diagnosis(input$age_diagnosis)
    
    # Retrieve the relevant coefficients
    age_diagnosis_coeff <- coeffs$age_diagnosis[[categorized_age]]
    sex_coeff <- coeffs$sex[[input$sex]]
    anthracycline_coeff <- coeffs$anthracycline[[input$anthracycline]]
    radiation_coeff <- coeffs$radiation[[input$radiation]]
    age_baseline_coeff <- coeffs$age_baseline[[input$age_baseline]]
    hypertension_coeff <- coeffs$hypertension[[input$hypertension]]
    ancestry_coeff <- coeffs$ancestry[[input$ancestry]]
    
    # Incorporate the Polygenic Risk Scores
    pr_score_coeff1 <- input$pr_score1 * -0.1316  # Adjust based on model's impact
    pr_score_coeff2 <- input$pr_score2 * 0.1361   # Adjust based on model's impact
    
    # Calculate the overall log risk, starting from the intercept
    log_risk <- -4.7265 + age_diagnosis_coeff + sex_coeff + anthracycline_coeff +
      radiation_coeff + age_baseline_coeff + hypertension_coeff +
      ancestry_coeff + pr_score_coeff1 + pr_score_coeff2 + log(input$follow_up_years / 10) 
    
    # Convert log risk back to relative risk (RR)
    risk <- exp(log_risk)
    return(risk)
  }
  
  # Calculate the risk for the selected inputs
  output$risk_output <- renderText({
    risk <- calculate_risk()
    paste("The predicted relative risk is: ", round(risk, 4))
  })
  
  # Display risk level based on predicted risk
  output$risk_level <- renderUI({
    risk <- calculate_risk()
    
    if (risk < 1.5) {
      tagList(span("Predicted risk level: Low", style = "color: green;"))
    } else if (risk < 3.0) {
      tagList(span("Predicted risk level: Moderate", style = "color: orange;"))
    } else {
      tagList(span("Predicted risk level: High", style = "color: red;"))
    }
  })
  
  # Plot the risk across different age stages
  output$risk_plot <- renderPlot({
    # Define the age stages
    age_stages <- c("5", "5_10", "10_15", "15")
    
    # Calculate risks for each stage
    risks <- sapply(age_stages, function(stage) {
      local_age_diagnosis <- stage
      coeffs <- list(
        "age_diagnosis" = c("5" = 0, "5_10" = -0.0382, "10_15" = -0.1181, "15" = -0.2631),
        "sex" = c("female" = 0, "male" = 0.4457),
        "anthracycline" = c("none" = 0, "0_100" = -0.0157, "100_250" = 0.9204, "250" = 1.7179),
        "radiation" = c("none" = 0, "5" = -0.1228, "5_15" = 0.3428, "15_35" = 0.9745, "35" = 2.818),
        "age_baseline" = c("25" = 0, "25_35" = 0.1775, "35_45" = 0.6221, "45" = 0.9118),
        "hypertension" = c("no" = 0, "yes" = 0.8247),
        "ancestry" = c("european" = 0, "african" = 0.6572, "others" = -0.5715)
      )
      
      age_diagnosis_coeff <- coeffs$age_diagnosis[[local_age_diagnosis]]
      sex_coeff <- coeffs$sex[[input$sex]]
      anthracycline_coeff <- coeffs$anthracycline[[input$anthracycline]]
      radiation_coeff <- coeffs$radiation[[input$radiation]]
      age_baseline_coeff <- coeffs$age_baseline[[input$age_baseline]]
      hypertension_coeff <- coeffs$hypertension[[input$hypertension]]
      ancestry_coeff <- coeffs$ancestry[[input$ancestry]]
      
      # Adjusted Polygenic Risk Scores
      pr_score_coeff1 <- input$pr_score1 * -0.1316  # Adjust based on model's impact
      pr_score_coeff2 <- input$pr_score2 * 0.1361   # Adjust based on model's impact
      
      log_risk <- -4.7265 + age_diagnosis_coeff + sex_coeff + anthracycline_coeff +
        radiation_coeff + age_baseline_coeff + hypertension_coeff +
        ancestry_coeff + pr_score_coeff1 + pr_score_coeff2 + log(input$follow_up_years / 10) 
      
      risk <- exp(log_risk)
      return(risk)
    })
    
    # Create a data frame for plotting
    risk_data <- data.frame(
      Age_Stage = factor(c("≤5", ">5-10", ">10-15", ">15"), 
                         levels = c("≤5", ">5-10", ">10-15", ">15")),  # Specify order here
      Risk = risks
    )
    
    # Fixed limits for y-axis (can be adjusted as needed)
    y_min <- 0
    y_max <- 20  # Adjust this value based on the expected maximum relative risk
    
    # Fancy plot with custom colors
    ggplot(risk_data, aes(x = Age_Stage, y = Risk, fill = Age_Stage)) +
      geom_bar(stat = "identity", color = "black", size = 0.5) +
      scale_fill_manual(values = c("#3498DB", "#9B59B6", "#E74C3C", "#2ECC71")) +
      scale_y_log10(limits = c(1, 100), breaks = c(1, 10, 100)) +  # Set y-axis limits and breaks
      theme_minimal() +
      labs(title = "Predicted Risk Across Age Stages", x = "Age Stage", y = "Relative Risk") +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16)
      ) +
      ylim(y_min, y_max) +  # Fixed y-axis limit
      geom_text(aes(label = round(Risk, 2)), vjust = -0.5, size = 5)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)





# # ## Data check
# library(readxl)
# model13 <- read_xlsx("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/Kateryna_CMP_risk_prediction/figure_2_all_new_13.xlsx", sheet = "Model_13")
# dfgene <- read.table("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/Kateryna_CMP_risk_prediction/df_gene.csv", sep = ",", header = T)
# dfgene$pred_est <- model13$pred_est[match(dfgene$sjlid, model13$sjlid)]
# 
# # SJL5085717; HCM = -3.316509E-7; LVEV = -1.283923176 
# SJL5085717 Works
# SJL5088113 Works
# SJL5162009 Works
# SJL1751607 works (Afr)
# SJL4175007 works (Afr)