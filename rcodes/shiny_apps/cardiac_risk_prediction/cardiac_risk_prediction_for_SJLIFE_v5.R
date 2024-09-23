# Load required packages
library(shiny)
library(shinythemes)  # For modern UI themes
library(shinyBS)      # For tooltips
library(ggplot2)

# Define UI for the stylish CMP prediction app
ui <- fluidPage(
  theme = shinytheme("cerulean"),  # Apply a modern theme
  
  tags$head(
    tags$style(HTML("
            .custom-title { font-size: 28px; font-weight: bold; color: #2C3E50; margin-bottom: 20px; }
            .custom-box { background-color: #ECF0F1; padding: 15px; border-radius: 10px; box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.1); }
            .custom-input { margin-bottom: 15px; }
            .risk-text { font-size: 24px; font-weight: bold; margin-top: 20px; }
            .low-risk { color: #27AE60; }
            .moderate-risk { color: #F39C12; }
            .high-risk { color: #E74C3E; }
        "))
  ),
  
  titlePanel("CMP Risk Prediction Model"),
  
  sidebarLayout(
    sidebarPanel(
      div(class = "custom-box",
          h3(class = "custom-title", "Input Variables"),
          
          # Input: Age at primary diagnosis (selectInput with tooltip)
          div(class = "custom-input",
              selectInput("age_diagnosis", "Age at Primary Diagnosis (years):",
                          choices = list("≤5" = "5", 
                                         ">5-10" = "5_10", 
                                         ">10-15" = "10_15", 
                                         ">15" = "15")),
              bsTooltip("age_diagnosis", "Select the age when the primary diagnosis occurred.", 
                        "right", options = list(container = "body"))
          ),
          
          # Input: Sex (select input)
          div(class = "custom-input",
              selectInput("sex", "Sex:", 
                          choices = list("Female" = "female", 
                                         "Male" = "male")),
              bsTooltip("sex", "Select the biological sex of the patient.", 
                        "right", options = list(container = "body"))
          ),
          
          # Input: Cumulative anthracycline dose
          div(class = "custom-input",
              selectInput("anthracycline", "Cumulative Anthracycline Dose (mg/m2):",
                          choices = list("None" = "none",
                                         ">0-100" = "0_100",
                                         ">100-250" = "100_250",
                                         ">250" = "250")),
              bsTooltip("anthracycline", "Select the total dose of anthracycline received.", 
                        "right", options = list(container = "body"))
          ),
          
          # Input: Mean heart radiation exposure dose
          div(class = "custom-input",
              selectInput("radiation", "Mean Heart Radiation Dose (Gray):",
                          choices = list("None" = "none",
                                         ">5" = "5",
                                         ">5-15" = "5_15",
                                         ">15-35" = "15_35",
                                         ">35" = "35")),
              bsTooltip("radiation", "Select the mean heart radiation dose.", 
                        "right", options = list(container = "body"))
          ),
          
          # Input: Age at baseline (select input)
          div(class = "custom-input",
              selectInput("age_baseline", "Age at Baseline (years):",
                          choices = list("≤25" = "25", 
                                         ">25-35" = "25_35", 
                                         ">35-45" = "35_45", 
                                         ">45" = "45")),
              bsTooltip("age_baseline", "Select the age at the time of baseline examination.", 
                        "right", options = list(container = "body"))
          )
      )
    ),
    
    mainPanel(
      div(class = "custom-box",
          h3(style = "font-size: 30px; font-weight: bold; color: black;", "Predicted CMP Risk"),
          textOutput("age_stage"),
          htmlOutput("risk_output"),
          uiOutput("risk_level"),
          plotOutput("risk_plot")
      )
    )
  )
)

# Define server logic for CMP risk prediction
server <- function(input, output) {
  
  # Dynamically show which age stage the user is in
  output$age_stage <- renderText({
    # Commented out as requested
    # age_stage <- switch(input$age_diagnosis, 
    #                     "5" = "Childhood", 
    #                     "5_10" = "Early Childhood", 
    #                     "10_15" = "Pre-Teen", 
    #                     "15" = "Teenager")
    # paste("Current stage: ", age_stage)
  })
  
  # Function to calculate risk based on age stage and other inputs
  calculate_risk <- function(age_stage_input) {
    # Coefficients from the provided table (log(RR))
    coeffs <- list(
      "age_diagnosis" = c("5" = 0, "5_10" = -0.0255, "10_15" = -0.1271, "15" = -0.2121),
      "sex" = c("female" = 0, "male" = 0.4407),
      "anthracycline" = c("none" = 0, "0_100" = -0.0726, "100_250" = 0.8329, "250" = 1.581),
      "radiation" = c("none" = 0, "5" = -0.0833, "5_15" = 0.3679, "15_35" = 0.9116, "35" = 2.6916),
      "age_baseline" = c("25" = 0, "25_35" = 0.2852, "35_45" = 0.8735, "45" = 1.2076)
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
    paste("<span style='font-size: 24px; font-weight: bold; color: blue;'>",
          "The predicted relative risk is: ", round(risk, 2), 
          "</span>", sep = "")
  })
  
  # Display risk level using color coding based on predicted risk
  output$risk_level <- renderUI({
    risk <- calculate_risk(input$age_diagnosis)
    
    # Color-coded risk level
    if (risk < 1.5) {
      tagList(
        span("Predicted risk level: Low", class = "risk-text low-risk")
      )
    } else if (risk < 3) {
      tagList(
        span("Predicted risk level: Moderate", class = "risk-text moderate-risk")
      )
    } else {
      tagList(
        span("Predicted risk level: High", class = "risk-text high-risk")
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
    
    # Fixed limits for y-axis
    y_min <- 0
    y_max <- 500  # Adjust this value based on the expected maximum relative risk
    
    # Fancy plot with custom colors
    ggplot(risk_data, aes(x = Age_Stage, y = Risk, fill = Age_Stage)) +
      geom_bar(stat = "identity", color = "black", size = 0.5) +
      scale_fill_manual(values = c("#3498DB", "#9B59B6", "#E74C3C", "#2ECC71")) +
      scale_y_log10(limits = c(1, 1000), breaks = c(1, 10, 100, 1000)) +  # Set y-axis limits and breaks
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
