# Load required packages
library(shiny)
library(ggplot2)

## SEX 1 is Male; 2 is Female

# Define UI for the CMP prediction app
ui <- fluidPage(
  titlePanel("Cardiomyopathy 10-year Risk Prediction"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("age_diagnosis", "Age at Primary Diagnosis (years):", value = 0, min = 0, max = 200, step = 1),
      # numericInput("follow_up_years", "Follow-up years:", value = 0, min = 0, max = 200, step = 1),
      selectInput("sex", "Sex:", choices = list("Female" = "female", "Male" = "male")),
      selectInput("anthracycline", "Cumulative Anthracycline Dose (mg/m²):",
                  choices = list("None" = "none", ">0-100" = "0_100", ">100-250" = "100_250", ">250" = "250")),
      selectInput("radiation", "Mean Heart Radiation Dose (Gray):",
                  choices = list("None" = "none", ">5" = "5", ">5-15" = "5_15", ">15-35" = "15_35", ">35" = "35")),
      selectInput("age_baseline", "Age at Baseline (years):",
                  choices = list("≤25" = "25", ">25-35" = "25_35", ">35-45" = "35_45", ">45" = "45")),
      selectInput("hypertension", "Hypertension:", choices = list("No" = "no", "Yes" = "yes")),
      selectInput("ancestry", "Genetic Ancestry:",
                  choices = list("European" = "european", "African" = "african", "Other" = "others")),
      numericInput("pr_score1", "Polygenic Risk Score HCM:", value = 0, min = -3.167714, max = 4.547851, step = 0.001),
      numericInput("pr_score2", "Polygenic Risk Score LVEVi:", value = 0, min = -4.668462, max = 3.601426, step = 0.001),
      actionButton("predict", "Predict")
    ),
    
    # mainPanel(
    #   h3("Cardiomyopathy 10-year Risk Prediction"),
    #   textOutput("prob_output"),
    #   textOutput("rr_output"),
    #   uiOutput("risk_level"),
    #   plotOutput("prob_plot")
    # )
    mainPanel(
      h3("Cardiomyopathy 10-year risk prediction"),
      uiOutput("prob_output"),  # Use uiOutput instead of textOutput
      textOutput("rr_output"),  
      uiOutput("risk_level"),
      plotOutput("prob_plot")
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
  
  # Function to calculate predicted event probability and relative risk
  calculate_prob_and_rr <- function() {
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
    
    ## Remove offset for 10-year prediction; removing log(input$follow_up_years / 10) will give 10-year estimate
    # Calculate the log risk and relative risk; 
    log_risk <- -4.7265 + age_diagnosis_coeff + sex_coeff + anthracycline_coeff +
      radiation_coeff + age_baseline_coeff + hypertension_coeff +
      ancestry_coeff + pr_score_coeff1 + pr_score_coeff2
    
    relative_risk <- exp(log_risk)
    probability <- 1 - exp(-exp(log_risk))
    
    return(list(probability = probability, rr = relative_risk))
  }
  
  observeEvent(input$predict, {
    # Make sure all inputs are provided before calculating
    if (is.null(input$age_diagnosis) || 
        is.null(input$sex) || 
        is.null(input$anthracycline) ||
        is.null(input$radiation) || 
        is.null(input$age_baseline) || 
        is.null(input$hypertension) ||
        is.null(input$ancestry) ||
        is.null(input$pr_score1) || 
        is.null(input$pr_score2)) {
      output$prob_output <- renderText("Please provide all inputs to get the prediction.")
      return()
    }
    
    # Perform the prediction
    results <- calculate_prob_and_rr()
    
    # # Display predicted probability
    # output$prob_output <- renderText({
    #   paste("The predicted probability is: ", round(results$probability * 100, 2), "%")
    # })
    
    output$prob_output <- renderUI({
      results <- calculate_prob_and_rr()  # Assuming this function is defined
      prob_text <- paste("The predicted probability is: ", round(results$probability * 100, 2), "%")
      tags$span(prob_text, style = "font-size: 20px; font-weight: bold;")
    })
    
    # # Calculate and display the relative risk; comment out this to not display risk
    # output$rr_output <- renderText({
    #   results <- calculate_prob_and_rr()
    #   paste("The predicted risk is: ", round(results$rr, 4))
    # })
    
    # # Display risk level
    # output$risk_level <- renderUI({
    #   prob <- results$probability
    #   
    #   if (prob < 0.015) {
    #     tagList(span("Predicted risk level: Low", style = "color: green; font-size:18px;"))
    #   } else if (prob < 0.03) {
    #     tagList(span("Predicted risk level: Moderate", style = "color: orange; font-size:18px;"))
    #   } else {
    #     tagList(span("Predicted risk level: High", style = "color: red; font-size:18px;"))
    #   }
    # })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

## deploy
# library(rsconnect)
# rsconnect::setAccountInfo(name='ysapkota',
#                           token='EB411EA6A3322D871FC7B27306F85FD0',
#                           secret='cZ7GfWe0uG/7xoFgIpsSzlVNMqM+voTwU5fHLKNs')

# rsconnect::deployApp('Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/shiny_apps/cardiac_risk_prediction/Cardiomyopathy_Risk_Calculator/', appName = 'Cardiomyopathy_Risk_Calculator')

# convert_to_pred_est <- function(probability) {
#   # Ensure probability is between 0 and 1
#   if (probability >= 1 || probability <= 0) {
#     stop("Probability must be between 0 and 1")
#   }
#   # Calculate pred_est
#   pred_est <- log(-log(1 - probability))
#   return(pred_est)
# }
# 
# # Example usage:
# probability <- 0.8  # Replace with your probability value
# pred_est <- convert_to_pred_est(probability)
# pred_est

