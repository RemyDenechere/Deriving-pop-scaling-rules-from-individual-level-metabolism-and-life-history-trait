library(latex2exp)
library(shiny)
library(ggplot2)
# source("functions.R")

rmax = function(A_0, b, M_0scl = "scaling"){
    n = 3/4
    a = 0.42
    epsR = 0.03
    M_0 = 0.01
    
    M = exp(seq(log(M_0), log(1000000), length.out=1000))
    
    if (M_0scl == "constant"){
        ratio = exp(10.299)* M
    } else if(M_0scl == "scaling") {
        ratio =  1612106 # based on elasmobtanches
    }
    
    rmax = A_0 * M^b * (1 - n) * M^(n-1) * ((1-a)*log(ratio) + log(epsR) )
    
    res = cbind.data.frame(rmax, M)
    return(res)  
}

plot_rmax = function(res){
    
    ggplot(res, aes(x = M, y = rmax))+
        theme_bw(base_size=18)+
        geom_line(size=2.5)+ 
        xlab("Adult mass (g)")+
        ylab("r_max (1/yr)")+
        scale_x_log10(limits = c(0.001, 1000000))+
        scale_y_log10(limits = c(10^(-2), 10^2))+ 
        geom_vline(xintercept = 0.01, linetype = "dotted", 
                   color = "blue", size=2.5)+
        annotate("text", 0.1 , 50, label = "Constant \n offspring size",
                 color= "blue", size = 6.5)+
        geom_line(aes(x = M, y = 1*M^(-1/4)), linetype = "dashed", 
                  size = 2.5)+
        annotate("text", 100000 , 0.02, label = "n-1 scaling \n from MTE",
                 color= "black", size = 6.5)
}


# 
## Define UI
#
ui <- fluidPage(
    titlePanel("Max population growth rate scaling (Denechere et al., 2022 Am. Nat.)"),
    
    withMathJax(),
    # section below allows in-line LaTeX via $ in mathjax.
    
    sidebarLayout(
        sidebarPanel(
            helpText('Calculation of maximum population growth rate rmax (invasion rate) 
                    for a size-structure population. Rmax depends on 3 main traits: Adult size M,
                    the somatic growth rate A, and the adult:offspring size ratio M/M_0.'),
            helpText('$$r_{max} \\approx -AM^{n-1} (0.6ln(M/M_0) + ln(\\epsilon_R))/4 $$ '),
            helpText(' with $$ A = A_0 M^{b}$$'),
            helpText('More info: https://doi.org/10.1086/718642'),
            
            # Input: Slider for A_0
            sliderInput(inputId = "A_0",
                        label = " Constant for somatic growth, A0",
                        min = 0.1,
                        max = 10,
                        value = 20),
            
            # Input: Slider for b
            sliderInput(inputId = "b",
                        label = "Exponent for somatic growth, b",
                        min = 0,
                        max = 0.3,
                        value = 10),
            
            # Input: Selector for choosing type of offspring size scaling 
            selectInput(inputId = "M_0scl",
                        label = "Scaling of M_0 (offspring) with M (adult) mass",
                        choices = c("constant", "scaling"))
            
            
        ),
        
        # we describe what function runs in the 
        mainPanel(
            plotOutput("plot_rmax", height="570px", width = "700px")
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
