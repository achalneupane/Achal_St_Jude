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
      h3("Predicted CMP Probability"),
      textOutput("prob_output"),
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
  
  # Function to calculate predicted event probability based on inputs
  calculate_prob <- function() {
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
    
    # Convert log risk into probability
    probability <- 1 - exp(-exp(log_risk))  # Probability using 1 - exp(-pred_est)
    
    return(probability)
  }
  
  # Calculate the probability for the selected inputs
  output$prob_output <- renderText({
    prob <- calculate_prob()
    paste("The predicted probability is: ", round(prob * 100, 2), "%")  # Showing up to 4 decimal places
  })
  
  # Display risk level based on predicted probability
  output$risk_level <- renderUI({
    prob <- calculate_prob()
    
    if (prob < 0.015) {
      tagList(span("Predicted risk level: Low", style = "color: green; font-size:18px;"))
    } else if (prob < 0.03) {
      tagList(span("Predicted risk level: Moderate", style = "color: orange; font-size:18px;"))
    } else {
      tagList(span("Predicted risk level: High", style = "color: red; font-size:18px;"))
    }
  })
  
  # Plot the probability across different age stages
  output$prob_plot <- renderPlot({
    # Define the age stages
    age_stages <- c("5", "5_10", "10_15", "15")
    
    # Calculate probabilities for each stage
    probs <- sapply(age_stages, function(stage) {
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
      
      probability <- 1 - exp(-exp(log_risk))
      return(probability)
    })
    
    # Create a data frame for plotting
    prob_data <- data.frame(
      Age_Stage = factor(c("≤5", ">5-10", ">10-15", ">15"), 
                         levels = c("≤5", ">5-10", ">10-15", ">15")),  # Specify order here
      Probability = probs * 100  # Convert to percentage
    )
    
    y_min <- 0
    y_max <- 100  # Adjust this value based on the 
    
    # Fancy plot with custom colors
    ggplot(prob_data, aes(x = Age_Stage, y = Probability, fill = Age_Stage)) +
      geom_bar(stat = "identity", color = "black", size = 0.5) +
      scale_fill_manual(values = c("#3498DB", "#9B59B6", "#E74C3C", "#2ECC71")) +
      scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Show percentages
      ylim(y_min, y_max) +  # Fixed y-axis limit
      geom_text(aes(label = round(Probability, 2)), vjust = -0.5, size = 5) +
      labs(x = "Age Strata", y = "Probability (%)", fill = "Age Strata") +  # Set axis and legend titles
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 18),  # Increase X-axis title size
        axis.title.y = element_text(size = 18),  # Increase Y-axis title size
        axis.text.x = element_text(size = 14),   # Increase X-axis text size
        axis.text.y = element_text(size = 14),   # Increase Y-axis text size
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Increase plot title size
        legend.text = element_text(size = 14),   # Increase legend text size
        legend.title = element_text(size = 20)   # Increase legend title size
      )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)





# # ## Data check
# library(readxl)
# model13 <- read_xlsx("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/Kateryna_CMP_risk_prediction/figure_2_all_new_13.xlsx", sheet = "Model_13")
# dfgene <- read.table("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Papers/Kateryna_CMP_risk_prediction/df_gene.csv", sep = ",", header = T)
# dfgene$pred_est <- model13$pred_est[match(dfgene$sjlid, model13$sjlid)]

## SJL5085717; HCM = -3.316509E-7; LVEV = -1.283923176