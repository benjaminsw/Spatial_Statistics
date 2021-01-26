#################################################
#### S1Y Rshiny app for visualising distributions
#################################################

library(shiny)
library(geoR)
library(spatstat)
ui <- fluidPage(title = "Point processes",
                tabsetPanel( 
                  #### HPP
                  tabPanel(title = "Homogeneous Poisson process (HPP)",
                           h3("This app illustrates the behaviour of a HPP."),
                           sliderInput("N.hpp", label= "Number of points", min = 1, max = 500, value = 100,step = 1),
                           actionButton(inputId = "clicks.hpp", 
                                        label = "Generate data"),
                           plotOutput("points.hpp"),
                           plotOutput("k.hpp")                           
                      ),
                  
                  
                  #### IPP
                  tabPanel(title = "Inhomogeneous Poisson process (IPP)",
                           h3("This app illustrates the behaviour of an IPP."),
                           sliderInput("N.ipp", label= "Mean intensity", min = 1, max = 8, value = 4,step = 0.01),
                           actionButton(inputId = "clicks.ipp", 
                                        label = "Generate data"),
                           plotOutput("points.ipp"),
                           plotOutput("k.ipp")                           
                  )
                )
)

server <- function(input, output) {
  
  #### HPP
  hpp.sim <- eventReactive(input$clicks.hpp, {
    rpoispp(lambda=input$N.hpp, win=owin(xrange=c(0,1), yrange=c(0,1)))
    
  })
  output$points.hpp <- renderPlot({
    par(mfrow=c(1,2))
        plot(hpp.sim(), pch=19, main="Realisation of a HPP")
        plot(density(hpp.sim(), bw.diggle), main="Estimated intensity for the HPP")
    })
  
  output$k.hpp <- renderPlot({
    Kc <- Kest(hpp.sim(), correction="Ripley")
    plot(envelope(hpp.sim(), Kest,nsim=39), main="Estimated K function for the HPP")
  })

  #### IPP
  ipp.sim <- eventReactive(input$clicks.ipp, {
  rLGCP(model="exp", mu=input$N.ipp, param=list(var=1, scale=0.1), win=owin(xrange=c(0,1), yrange=c(0,1)))
  })

    output$points.ipp <- renderPlot({
    par(mfrow=c(1,2))
    plot(ipp.sim(), pch=19, main="Realisation of a IPP")
    plot(density(ipp.sim(), bw.diggle), main="Estimated intensity for the IPP")
  })
  
  output$k.ipp <- renderPlot({
    Kc <- Kest(ipp.sim(), correction="Ripley")
    plot(envelope(ipp.sim(), Kest,nsim=39), main="Estimated K function for the IPP")
  })
  
  
}

shinyApp(server = server, ui = ui)