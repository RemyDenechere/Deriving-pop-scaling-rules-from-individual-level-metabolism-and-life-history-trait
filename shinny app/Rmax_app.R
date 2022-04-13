library(latex2exp)
library(shiny)
source("functions.R")

# 
## Define UI
#
ui <- fluidPage(
  titlePanel("Max population growth rate scaling (Denechere et al., 2022 Am. Nat.)"),
  
  withMathJax(),
  # section below allows in-line LaTeX via $ in mathjax.

  sidebarLayout(
    sidebarPanel(
      helpText('$$r_{max} \\approx A (1-n) M^{n-1} ((1-a) ln(M/M_0) + ln(\\epsilon_R)) $$ '),
      helpText(' with $$ A = A_0 M^{b}$$'),
      
      # Input: Slider for A_0
      sliderInput(inputId = "A_0",
                  label = " Constant for somatic growth rate, A0",
                  min = 0.1,
                  max = 10,
                  value = 20),
      
      # Input: Slider for b
      sliderInput(inputId = "b",
                  label = "Exponent for somatic growth rate, b",
                  min = 0,
                  max = 0.3,
                  value = 10),
      
      # Input: Selector for choosing type of offspring size scaling 
      selectInput(inputId = "M_0scl",
                  label = "Choose the scaling between M_0 (offspring) and M (adult) mass",
                  choices = c("constant", "scaling"))
      
      
      ),
    
    # we describe what function runs in the 
    mainPanel(
      plotOutput("plot_rmax", height="500px", width = "600px")
    )
  )
)

#
# Define server logic
#
server <- function(input, output) {
  
  sim <- eventReactive({
    input$A_0
    input$b
    input$M_0scl
  },
  
  {
    # setup simulation
    res = rmax(A_0 = input$A_0, b = input$b, M_0scl = input$M_0scl)
    # Simulate
    return(res)
  }, 
  ignoreNULL = TRUE)
  
  # Make plots
  output$plot_rmax <- renderPlot(plot_rmax(sim()))
}

#
# Run the application 
#
shinyApp(ui = ui, server = server)
