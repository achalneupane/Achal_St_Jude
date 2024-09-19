library(shiny)
library(glmnet)
library(car)  # For VIF calculation
library(shinythemes)

# Load the Framingham Heart Study dataset (you can replace this with the actual dataset you're using)
data <- read.csv("https://raw.githubusercontent.com/GauravPadawe/Framingham-Heart-Study/refs/heads/master/framingham.csv")

# Preprocess data (remove NA values for simplicity)
data <- na.omit(data)
data$sex <- data$male 

# Fit a logistic regression model using relevant features
model <- glm(TenYearCHD ~ age + sex + totChol + sysBP + currentSmoker, family = binomial, data = data)

# Check for multicollinearity
vif_values <- vif(model)
print(vif_values)  # If VIF values are too high (> 5 or 10), consider dropping or combining variables

# Define UI for application
ui <- fluidPage(
  theme = shinytheme("flatly"),  # Apply a nice theme to the app
  titlePanel(
    div(
      h1("Cardiac Risk Prediction", style = "color: #007bff; font-weight: bold;"),
      p("Estimate your 10-year risk of heart disease", style = "font-style: italic;")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      h3("Enter Your Health Information", style = "color: #17a2b8;"),
      
      # Input Fields Styled
      numericInput("age", "Age (years):", value = 50, min = 18, max = 100),
      selectInput("sex", "Gender:", choices = list("Male" = 1, "Female" = 0), selected = 1),
      numericInput("cholesterol", "Total Cholesterol (mg/dL):", value = 200, min = 100, max = 400),
      numericInput("sysBP", "Systolic Blood Pressure (mmHg):", value = 120, min = 80, max = 200),
      checkboxInput("currentSmoker", "Are you a current smoker?", FALSE),
      
      actionButton("predict", "Predict Risk", class = "btn btn-primary btn-lg btn-block"),
      br(),
      br(),
      div(icon("heartbeat", class = "fa-4x", style = "color: #dc3545; text-align: center;"))
    ),
    
    mainPanel(
      h3("Your Predicted Cardiac Risk", style = "color: #17a2b8;"),
      div(
        verbatimTextOutput("risk_output"),
        style = "background-color: #f8f9fa; padding: 20px; border-radius: 10px; font-size: 1.5em;"
      ),
      
      # Adding progress bar for visual effect
      uiOutput("progress_ui"),
      br(),
      p("This app is for informational purposes only. Please consult your healthcare provider for accurate diagnosis and medical advice.", 
        style = "font-size: 0.9em; font-style: italic; color: gray;")
    )
  )
)

# Define server logic
server <- function(input, output) {
  # Define reactive expression for prediction
  predict_risk <- eventReactive(input$predict, {
    # Capture user inputs
    new_data <- data.frame(
      age = input$age,
      sex = as.numeric(input$sex),  # Convert to numeric for modeling
      totChol = input$cholesterol,
      sysBP = input$sysBP,
      currentSmoker = as.numeric(input$currentSmoker)  # Convert logical to numeric
    )
    
    # Predict probability of cardiac event
    predicted_prob <- predict(model, newdata = new_data, type = "response")
    
    # Ensure prediction stays within 0-100% range
    predicted_prob <- pmin(pmax(predicted_prob, 0), 1)
    
    return(predicted_prob)
  })
  
  # Render the predicted risk when the button is clicked
  output$risk_output <- renderText({
    risk <- predict_risk()
    paste("Your estimated 10-year risk of developing heart disease is:", round(risk * 100, 2), "%")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
