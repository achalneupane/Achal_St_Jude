library(shiny)

# Define UI for the application
ui <- fluidPage(
  titlePanel(tags$h5("This risk assessment tool predicts the 10-year risk of cardiomyopathy among long-term survivors of childhood cancer. It uses information from the “Predicting 10-year risk of cardiomyopathy in long-term survivors of childhood cancer: a report from the St. Jude Lifetime Cohort and the Childhood Cancer Survivor Study” available at https://www.medrxiv.org/content/10.1101/2024.10.24.24316064v1, which created the prediction models with the highest-to-date performance, specifically designed for identifying risk of cardiomyopathy in long-term survivors of childhood cancer based on survivor characteristics, treatment exposures, and inherited genetic variation. Initially developed using data from the well-characterized and clinically assessed cohort of childhood cancer survivors in the St. Jude Lifetime (SJLIFE) cohort, these models were validated in the independent cohort of childhood cancer survivors from the Childhood Cancer Survivor Study (CCSS). Each model offers clinically relevant, easy-to-apply predictions, progressively incorporating more risk factors to refine the risk estimate according to the available clinical data.")),
  
  tabsetPanel(
    tabPanel("Model 1",
             titlePanel(tags$h5("This model estimates the 10-year risk of cardiomyopathy based on age at primary childhood cancer diagnosis, sex, and treatment exposures such as cumulative anthracycline dose (mg/m²) and average heart radiation dose (Gray).")),
             sidebarLayout(
               sidebarPanel(
                 numericInput("age_diagnosis5", "Age at childhood cancer diagnosis (>0 to ≤23.5 years):", value = 0, min = 0, max = 200, step = 1),
                 selectInput("sex5", "Sex:", choices = list("Female" = "female", "Male" = "male")),
                 selectInput("anthracycline5", "Cumulative anthracycline dose (mg/m²):",
                             choices = list("None" = "none", ">0-100" = "0_100", ">100-250" = "100_250", ">250" = "250")),
                 selectInput("radiation5", "Average heart radiation dose (Gray):",
                             choices = list("None" = "none", ">5" = "5", ">5-15" = "5_15", ">15-35" = "15_35", ">35" = "35")),
                 actionButton("calculate5", "Calculate"),
                 
                 
                 actionButton("reset_calc1", "Reset")
                 
                 
               ),
               mainPanel(uiOutput("prob_output5"))
             )),
    
    # Tab 2
    tabPanel("Model 2",
             titlePanel(tags$h5("This model estimates the 10-year risk of cardiomyopathy based on age at primary childhood cancer diagnosis, sex, treatment exposures such as cumulative anthracycline dose (mg/m²) and average heart radiation dose (Gray), and patient’s current age (age at assessment).")),
             sidebarLayout(
               sidebarPanel(
                 numericInput("age_diagnosis4", "Age at childhood cancer diagnosis (>0 to ≤23.5 years):", value = 0, min = 0, max = 200, step = 1),
                 selectInput("sex4", "Sex:", choices = list("Female" = "female", "Male" = "male")),
                 selectInput("anthracycline4", "Cumulative anthracycline dose (mg/m²):",
                             choices = list("None" = "none", ">0-100" = "0_100", ">100-250" = "100_250", ">250" = "250")),
                 selectInput("radiation4", "Average heart radiation dose (Gray):",
                             choices = list("None" = "none", ">5" = "5", ">5-15" = "5_15", ">15-35" = "15_35", ">35" = "35")),
                 selectInput("age_baseline4", "Current age (years):",
                             choices = list("≤25" = "25", ">25-35" = "25_35", ">35-45" = "35_45", ">45" = "45")),
                 actionButton("calculate4", "Calculate"),
                 
                 
                 actionButton("reset_calc2", "Reset")
               ),
               mainPanel(uiOutput("prob_output4"))
             )),
    
    # Tab 3
    tabPanel("Model 3",
             titlePanel(tags$h5("This model estimates the 10-year risk of cardiomyopathy based on age at primary childhood cancer diagnosis, sex, treatment exposures such as cumulative anthracycline dose (mg/m²) and average heart radiation dose (Gray), patient’s current age (age at assessment), and hypertension status (yes/no).")),
             sidebarLayout(
               sidebarPanel(
                 numericInput("age_diagnosis3", "Age at childhood cancer diagnosis (>0 to ≤23.5 years):", value = 0, min = 0, max = 200, step = 1),
                 selectInput("sex3", "Sex:", choices = list("Female" = "female", "Male" = "male")),
                 selectInput("anthracycline3", "Cumulative anthracycline dose (mg/m²):",
                             choices = list("None" = "none", ">0-100" = "0_100", ">100-250" = "100_250", ">250" = "250")),
                 selectInput("radiation3", "Average heart radiation dose (Gray):",
                             choices = list("None" = "none", ">5" = "5", ">5-15" = "5_15", ">15-35" = "15_35", ">35" = "35")),
                 selectInput("age_baseline3", "Current age (years):",
                             choices = list("≤25" = "25", ">25-35" = "25_35", ">35-45" = "35_45", ">45" = "45")),
                 selectInput("hypertension3", "Hypertension:", choices = list("No" = "no", "Yes" = "yes")),
                 actionButton("calculate3", "Calculate"),
                 
                 
                 actionButton("reset_calc3", "Reset")
               ),
               mainPanel(uiOutput("prob_output3"))
             )),
    
    # Tab 4
    tabPanel("Model 4",
             titlePanel(tags$h5("This model estimates the 10-year risk of cardiomyopathy based on age at primary childhood cancer diagnosis, sex, treatment exposures such as cumulative anthracycline dose (mg/m²) and average heart radiation dose (Gray), patient’s current age (age at assessment), hypertension status (yes/no), and genetic ancestry (European/African/Other).")),
             sidebarLayout(
               sidebarPanel(
                 numericInput("age_diagnosis2", "Age at childhood cancer diagnosis (>0 to ≤23.5 years):", value = 0, min = 0, max = 200, step = 1),
                 selectInput("sex2", "Sex:", choices = list("Female" = "female", "Male" = "male")),
                 selectInput("anthracycline2", "Cumulative anthracycline dose (mg/m²):",
                             choices = list("None" = "none", ">0-100" = "0_100", ">100-250" = "100_250", ">250" = "250")),
                 selectInput("radiation2", "Average heart radiation dose (Gray):",
                             choices = list("None" = "none", ">5" = "5", ">5-15" = "5_15", ">15-35" = "15_35", ">35" = "35")),
                 selectInput("age_baseline2", "Current age (years):",
                             choices = list("≤25" = "25", ">25-35" = "25_35", ">35-45" = "35_45", ">45" = "45")),
                 selectInput("hypertension2", "Hypertension:", choices = list("No" = "no", "Yes" = "yes")),
                 selectInput("ancestry2", "Genetic ancestry:",
                             choices = list("European" = "european", "African" = "african", "Other" = "others")),
                 actionButton("calculate2", "Calculate"),
                 
                 
                 actionButton("reset_calc4", "Reset")
               ),
               mainPanel(uiOutput("prob_output2"))
             )),
    
    # Tab 5
    tabPanel("Model 5",
             titlePanel(tags$h5("This model estimates the 10-year risk of cardiomyopathy based on age at primary childhood cancer diagnosis, sex, treatment exposures such as cumulative anthracycline dose (mg/m²) and average heart radiation dose (Gray), patient’s current age (age at assessment), hypertension status (yes/no), genetic ancestry (European/African/Other), and general population polygenic risk scores (PRS) for hypertrophic cardiomyopathy and left ventricular end-systolic volume index (z-score).")),
             sidebarLayout(
               sidebarPanel(
                 numericInput("age_diagnosis1", "Age at childhood cancer diagnosis (>0 to ≤23.5 years):", value = 0, min = 0, max = 200, step = 1),
                 selectInput("sex1", "Sex:", choices = list("Female" = "female", "Male" = "male")),
                 selectInput("anthracycline1", "Cumulative anthracycline dose (mg/m²):",
                             choices = list("None" = "none", ">0-100" = "0_100", ">100-250" = "100_250", ">250" = "250")),
                 selectInput("radiation1", "Average heart radiation dose (Gray):",
                             choices = list("None" = "none", ">5" = "5", ">5-15" = "5_15", ">15-35" = "15_35", ">35" = "35")),
                 selectInput("age_baseline1", "Current age (years):",
                             choices = list("≤25" = "25", ">25-35" = "25_35", ">35-45" = "35_45", ">45" = "45")),
                 selectInput("hypertension1", "Hypertension:", choices = list("No" = "no", "Yes" = "yes")),
                 selectInput("ancestry1", "Genetic ancestry:",
                             choices = list("European" = "european", "African" = "african", "Other" = "others")),
                 numericInput("pr_score11", "Polygenic risk score hypertrophic cardiomyopathy (z-score):", value = 0, min = -3, max = 4, step = 0.001),
                 numericInput("pr_score21", "Polygenic risk score left ventricular end-systolic volume index (z-score):", value = 0, min = -4, max = 3, step = 0.001),
                 actionButton("calculate1", "Calculate"),
                 
                 
                 actionButton("reset_calc5", "Reset")
               ),
               mainPanel(uiOutput("prob_output1"))
             ))
  )
)


# Define server logic
server <- function(input, output, session) {
  # Define intercepts
  var_1 <- -4.7265
  var_2 <- -4.7004
  var_3 <- -4.6164
  var_4 <- -4.5172
  var_5 <- -4.0984
  
  # Define coefficients for each calculator
  coeffs1 <- list(
    age_diagnosis = c("5" = 0, "5_10" = -0.0382, "10_15" = -0.1181, "15" = -0.2631),
    sex = c("female" = 0, "male" = 0.4457),
    anthracycline = c("none" = 0, "0_100" = -0.0157, "100_250" = 0.9204, "250" = 1.7179),
    radiation = c("none" = 0, "5" = -0.1228, "5_15" = 0.3428, "15_35" = 0.9745, "35" = 2.8180),
    age_baseline = c("25" = 0, "25_35" = 0.1775, "35_45" = 0.6221, "45" = 0.9118),
    hypertension = c("no" = 0, "yes" = 0.8247),
    ancestry = c("european" = 0, "african" = 0.6572, "others" = -0.5715)
  )
  
  coeffs2 <- list(
    # Replace with actual coefficients for Calculator 2
    age_diagnosis = c("5" = 0, "5_10" = -0.0227, "10_15" = -0.1097, "15" = -0.2367),
    sex = c("female" = 0, "male" = 0.4308),
    anthracycline = c("none" = 0, "0_100" = -0.0245, "100_250" = 0.9065, "250" = 1.6881),
    radiation = c("none" = 0, "5" = -0.1195, "5_15" = 0.3516, "15_35" = 0.9623, "35" = 2.7599),
    age_baseline = c("25" = 0, "25_35" = 0.2128, "35_45" = 0.6738, "45" = 0.9024),
    hypertension = c("no" = 0, "yes" = 0.7835),
    ancestry = c("european" = 0, "african" = 0.5298, "others" = -0.6232)
  )
  
  coeffs3 <- list(
    # Replace with actual coefficients for Calculator 3
    age_diagnosis = c("5" = 0, "5_10" = -0.0279, "10_15" = -0.0962, "15" = -0.2376),
    sex = c("female" = 0, "male" = 0.4346),
    anthracycline = c("none" = 0, "0_100" = -0.0666, "100_250" = 0.8594, "250" = 1.6396),
    radiation = c("none" = 0, "5" = -0.1147, "5_15" = 0.3533, "15_35" = 0.9757, "35" = 2.8130),
    age_baseline = c("25" = 0, "25_35" = 0.2176, "35_45" = 0.6718, "45" = 0.9038),
    hypertension = c("no" = 0, "yes" = 0.7870),
    
    ancestry = c("european" = 0, "african" = 0, "others" = 0)
  )
  
  coeffs4 <- list(
    # Replace with actual coefficients for Calculator 4
    age_diagnosis = c("5" = 0, "5_10" = -0.0255, "10_15" = -0.1271, "15" = -0.2121),
    sex = c("female" = 0, "male" = 0.4407),
    anthracycline = c("none" = 0, "0_100" = -0.0721, "100_250" = 0.8328, "250" = 1.5813),
    radiation = c("none" = 0, "5" = -0.0833, "5_15" = 0.3679, "15_35" = 0.9116, "35" = 2.6916),
    age_baseline = c("25" = 0, "25_35" = 0.2852, "35_45" = 0.8735, "45" = 1.2076),
    
    hypertension = c("no" = 0, "yes" = 0),
    ancestry = c("european" = 0, "african" = 0, "others" = 0)
  )
  
  coeffs5 <- list(
    # Replace with actual coefficients for Calculator 5
    age_diagnosis = c("5" = 0, "5_10" = 0.0596, "10_15" = 0.1454, "15" = 0.1873),
    sex = c("female" = 0, "male" = 0.4271),
    anthracycline = c("none" = 0, "0_100" = -0.2635, "100_250" = 0.5825, "250" = 1.4427),
    radiation = c("none" = 0, "5" = 0.0842, "5_15" = 0.4515, "15_35" = 0.9618, "35" = 2.8815),
    
    age_baseline = c("25" = 0, "25_35" = 0, "35_45" = 0, "45" = 0),
    hypertension = c("no" = 0, "yes" = 0),
    ancestry = c("european" = 0, "african" = 0, "others" = 0)
  )
  
  # Function to calculate risk
  calculate_risk <- function(age_diagnosis, sex, anthracycline, radiation, age_baseline = NULL, hypertension = NULL, ancestry = NULL, pr_score1 = NULL, pr_score2 = NULL, var, coeffs) {
    # Categorize age_diagnosis into radiationges
    categorized_age_diagnosis <- function(age) {
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
    
    categorized_age <- categorized_age_diagnosis(age_diagnosis)
    
    # Retrieve relevant coefficients
    age_diagnosis_coeff <- coeffs$age_diagnosis[[categorized_age]]
    sex_coeff <- coeffs$sex[[sex]]
    anthracycline_coeff <- coeffs$anthracycline[[anthracycline]]
    radiation_coeff <- coeffs$radiation[[radiation]]
    age_baseline_coeff <- if (!is.null(age_baseline)) coeffs$age_baseline[[age_baseline]] else 0
    hypertension_coeff <- if (!is.null(hypertension)) coeffs$hypertension[[hypertension]] else 0
    ancestry_coeff <- if (!is.null(ancestry)) coeffs$ancestry[[ancestry]] else 0
    pr_score1_coeff <- if (!is.null(pr_score1)) pr_score1 * -0.131 else 0
    pr_score2_coeff <- if (!is.null(pr_score2)) pr_score2 * 0.136 else 0
    
    # Calculate log risk, relative risk, and probability
    log_risk <- var + age_diagnosis_coeff + sex_coeff + anthracycline_coeff +
      radiation_coeff + age_baseline_coeff + hypertension_coeff + ancestry_coeff +
      pr_score1_coeff + pr_score2_coeff
    
    relative_risk <- exp(log_risk)
    probability <- 1 - exp(-exp(log_risk))
    
    return(list(probability = probability, rr = relative_risk))
  }
  
  # Observers for each tab
  observeEvent(input$calculate1, {
    results <- calculate_risk(
      age_diagnosis = input$age_diagnosis1,
      sex = input$sex1,
      anthracycline = input$anthracycline1,
      radiation = input$radiation1,
      age_baseline = input$age_baseline1,
      hypertension = input$hypertension1,
      ancestry = input$ancestry1,
      pr_score1 = input$pr_score11,
      pr_score2 = input$pr_score21,
      var = var_1,
      coeffs = coeffs1  # Use coeffs1 for Calculator 1
    )
    output$prob_output1 <- renderUI({
      tags$span(paste0("The predicted probability of cardiomayopathy in the next 10 years is: ", round(results$probability * 100, 1), "%"),
                style = "font-size: 20px; font-weight: bold;")
    })
  })
  
  observeEvent(input$calculate2, {
    results <- calculate_risk(
      age_diagnosis = input$age_diagnosis2,
      sex = input$sex2,
      anthracycline = input$anthracycline2,
      radiation = input$radiation2,
      age_baseline = input$age_baseline2,
      hypertension = input$hypertension2,
      ancestry = input$ancestry2,
      var = var_2,
      coeffs = coeffs2  # Use coeffs2 for Calculator 2
    )
    output$prob_output2 <- renderUI({
      tags$span(paste0("The predicted probability of cardiomayopathy in the next 10 years is: ", round(results$probability * 100, 1), "%"),
                style = "font-size: 20px; font-weight: bold;")
    })
  })
  
  observeEvent(input$calculate3, {
    results <- calculate_risk(
      age_diagnosis = input$age_diagnosis3,
      sex = input$sex3,
      anthracycline = input$anthracycline3,
      radiation = input$radiation3,
      age_baseline = input$age_baseline3,
      hypertension = input$hypertension3,
      var = var_3,
      coeffs = coeffs3  # Use coeffs3 for Calculator 3
    )
    output$prob_output3 <- renderUI({
      tags$span(paste0("The predicted probability of cardiomayopathy in the next 10 years is: ", round(results$probability * 100, 1), "%"),
                style = "font-size: 20px; font-weight: bold;")
    })
  })
  
  observeEvent(input$calculate4, {
    results <- calculate_risk(
      age_diagnosis = input$age_diagnosis4,
      sex = input$sex4,
      anthracycline = input$anthracycline4,
      radiation = input$radiation4,
      age_baseline = input$age_baseline4,
      var = var_4,
      coeffs = coeffs4  # Use coeffs4 for Calculator 4
    )
    output$prob_output4 <- renderUI({
      tags$span(paste0("The predicted probability of cardiomayopathy in the next 10 years is: ", round(results$probability * 100, 1), "%"),
                style = "font-size: 20px; font-weight: bold;")
    })
  })
  
  observeEvent(input$calculate5, {
    results <- calculate_risk(
      age_diagnosis = input$age_diagnosis5,
      sex = input$sex5,
      anthracycline = input$anthracycline5,
      radiation = input$radiation5,
      var = var_5,
      coeffs = coeffs5  # Use coeffs5 for Calculator 5
    )
    output$prob_output5 <- renderUI({
      tags$span(paste0("The predicted probability of cardiomayopathy in the next 10 years is: ", round(results$probability * 100, 1), "%"),
                style = "font-size: 20px; font-weight: bold;")
    })
  })
  
  
  
  # Reset functionality for Calculator 5
  observeEvent(input$reset_calc5, {
    updateNumericInput(session, "age_diagnosis1", value = 0)
    updateSelectInput(session, "sex1", selected = "female")
    updateSelectInput(session, "anthracycline1", selected = "none")
    updateSelectInput(session, "radiation1", selected = "none")
    updateSelectInput(session, "age_baseline1", selected = "25")
    updateSelectInput(session, "hypertension1", selected = "no")
    updateSelectInput(session, "ancestry1", selected = "european")
    updateNumericInput(session, "pr_score11", value = 0)
    updateNumericInput(session, "pr_score21", value = 0)
    
    # Clear the probability output
    output$prob_output1 <- renderUI({
      NULL  # Set output to NULL to make it disappear
    })
    
  })
  
  # Reset functionality for Calculator 4
  observeEvent(input$reset_calc4, {
    updateNumericInput(session, "age_diagnosis2", value = 0)
    updateSelectInput(session, "sex2", selected = "female")
    updateSelectInput(session, "anthracycline2", selected = "none")
    updateSelectInput(session, "radiation2", selected = "none")
    updateSelectInput(session, "age_baseline2", selected = "25")
    updateSelectInput(session, "hypertension2", selected = "no")
    updateSelectInput(session, "ancestry2", selected = "european")
    
    # Clear the probability output
    output$prob_output2 <- renderUI({
      NULL  # Set output to NULL to make it disappear
    })
  })
  
  # Reset functionality for Calculator 3
  observeEvent(input$reset_calc3, {
    updateNumericInput(session, "age_diagnosis3", value = 0)
    updateSelectInput(session, "sex3", selected = "female")
    updateSelectInput(session, "anthracycline3", selected = "none")
    updateSelectInput(session, "radiation3", selected = "none")
    updateSelectInput(session, "age_baseline3", selected = "25")
    updateSelectInput(session, "hypertension3", selected = "no")
    
    # Clear the probability output
    output$prob_output3 <- renderUI({
      NULL  # Set output to NULL to make it disappear
    })
  })
  
  # Reset functionality for Calculator 2
  observeEvent(input$reset_calc2, {
    updateNumericInput(session, "age_diagnosis4", value = 0)
    updateSelectInput(session, "sex4", selected = "female")
    updateSelectInput(session, "anthracycline4", selected = "none")
    updateSelectInput(session, "radiation4", selected = "none")
    updateSelectInput(session, "age_baseline4", selected = "25")
    
    # Clear the probability output
    output$prob_output4 <- renderUI({
      NULL  # Set output to NULL to make it disappear
    })
  })
  
  # Reset functionality for Calculator 1
  observeEvent(input$reset_calc1, {
    updateNumericInput(session, "age_diagnosis5", value = 0)
    updateSelectInput(session, "sex5", selected = "female")
    updateSelectInput(session, "anthracycline5", selected = "none")
    updateSelectInput(session, "radiation5", selected = "none")
    updateSelectInput(session, "age_baseline5", selected = "25")
    
    # Clear the probability output
    output$prob_output5 <- renderUI({
      NULL  # Set output to NULL to make it disappear
    })
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)

# # deploy
# library(rsconnect)
# rsconnect::setAccountInfo(name='ysapkota',
#                           token='EB411EA6A3322D871FC7B27306F85FD0',
#                           secret='cZ7GfWe0uG/7xoFgIpsSzlVNMqM+voTwU5fHLKNs')
# 
# rsconnect::deployApp('Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/shiny_apps/cardiac_risk_prediction/Cardiomyopathy_Risk_Calculator/', appName = 'Cardiomyopathy_Risk_Calculator')




