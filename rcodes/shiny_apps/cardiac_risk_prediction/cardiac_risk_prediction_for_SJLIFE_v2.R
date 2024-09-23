library(shiny)
library(ggplot2)

# Define UI for the CMP prediction app
ui <- fluidPage(
  titlePanel("CMP Risk Prediction"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Input Variables"),
      
      # Input: Age at primary diagnosis (using selectInput)
      selectInput("age_diagnosis", "Age at Primary Diagnosis (years):",
                  choices = list("≤5" = "5", 
                                 ">5-10" = "5_10", 
                                 ">10-15" = "10_15", 
                                 ">15" = "15")),
      
      # Input: Sex (select input)
      selectInput("sex", "Sex:", 
                  choices = list("Female" = "female", 
                                 "Male" = "male")),
      
      # Input: Cumulative anthracycline dose
      selectInput("anthracycline", "Cumulative Anthracycline Dose (mg/m2):",
                  choices = list("None" = "none",
                                 ">0-100" = "0_100",
                                 ">100-250" = "100_250",
                                 ">250" = "250")),
      
      # Input: Mean heart radiation exposure dose
      selectInput("radiation", "Mean Heart Radiation Dose (Gray):",
                  choices = list("None" = "none",
                                 ">5" = "5",
                                 ">5-15" = "5_15",
                                 ">15-35" = "15_35",
                                 ">35" = "35")),
      
      # Input: Age at baseline (select input)
      selectInput("age_baseline", "Age at Baseline (years):",
                  choices = list("≤25" = "25", 
                                 ">25-35" = "25_35", 
                                 ">35-45" = "35_45", 
                                 ">45" = "45"))
    ),
    
    mainPanel(
      h3("Predicted CMP Risk"),
      textOutput("age_stage"),
      verbatimTextOutput("risk_output"),
      uiOutput("risk_level"),  # This will show risk levels
      plotOutput("risk_plot")   # Output for the risk plot
    )
  )
)

# Define server logic for CMP risk prediction
server <- function(input, output) {
  
  # Dynamically show which age stage the user is in
  output$age_stage <- renderText({
    age_stage <- switch(input$age_diagnosis, 
                        "5" = "Childhood", 
                        "5_10" = "Early Childhood", 
                        "10_15" = "Pre-Teen", 
                        "15" = "Teenager")
    paste("Current stage: ", age_stage)
  })
  
  # Function to calculate risk based on age stage and other inputs
  calculate_risk <- function(age_stage_input) {
    # Coefficients from the provided table (log(RR))
    coeffs <- list(
      "age_diagnosis" = c("5" = 0, "5_10" = -0.0305, "10_15" = -0.1278, "15" = -0.2107),
      "sex" = c("female" = 0, "male" = 0.4383),
      "anthracycline" = c("none" = 0, "0_100" = -0.0726, "100_250" = 0.8329, "250" = 1.581),
      "radiation" = c("none" = 0, "5" = -0.0834, "5_15" = 0.3646, "15_35" = 0.9123, "35" = 2.6912),
      "age_baseline" = c("25" = 0, "25_35" = 0.2852, "35_45" = 0.8755, "45" = 1.209)
    )
    
    # Retrieve the relevant coefficients
    age_diagnosis_coeff <- coeffs$age_diagnosis[[age_stage_input]]
    sex_coeff <- coeffs$sex[[input$sex]]
    anthracycline_coeff <- coeffs$anthracycline[[input$anthracycline]]
    radiation_coeff <- coeffs$radiation[[input$radiation]]
    age_baseline_coeff <- coeffs$age_baseline[[input$age_baseline]]
    
    # Sum the coefficients (log(RR)) to get the overall log risk
    log_risk <- age_diagnosis_coeff + sex_coeff + anthracycline_coeff + radiation_coeff + age_baseline_coeff
    
    # Convert log risk back to relative risk (RR)
    risk <- exp(log_risk)
    return(risk)
  }
  
  # Calculate the risk for the selected age diagnosis stage
  output$risk_output <- renderText({
    risk <- calculate_risk(input$age_diagnosis)
    paste("The predicted relative risk is:", round(risk, 2))
  })
  
  # Display risk level using color coding based on predicted risk
  output$risk_level <- renderUI({
    risk <- calculate_risk(input$age_diagnosis)
    
    # Color-coded risk level
    if (risk < 1.5) {
      tagList(
        span("Predicted risk level: Low", style = "color: green; font-size: 18px;")
      )
    } else if (risk < 3) {
      tagList(
        span("Predicted risk level: Moderate", style = "color: orange; font-size: 18px;")
      )
    } else {
      tagList(
        span("Predicted risk level: High", style = "color: red; font-size: 18px;")
      )
    }
  })
  
  # Plot the risk across different age stages
  output$risk_plot <- renderPlot({
    # Define the age stages
    age_stages <- c("5", "5_10", "10_15", "15")
    
    # Calculate risks for each stage
    risks <- sapply(age_stages, calculate_risk)
    
    # Create a data frame for plotting
    risk_data <- data.frame(
      Age_Stage = c("≤5", ">5-10", ">10-15", ">15"),
      Risk = risks
    )
    
    # Plot the risks
    ggplot(risk_data, aes(x = Age_Stage, y = Risk)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(title = "Predicted Risk Across Age Stages", x = "Age Stage", y = "Relative Risk") +
      geom_text(aes(label = round(Risk, 2)), vjust = -0.5, size = 5)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
