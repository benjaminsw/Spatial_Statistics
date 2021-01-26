#  Random fieds generator  App
library(shiny)
library(geoR)

ui <- fluidPage(title = "Variograms and spatial parameters",
sidebarLayout(position = "left",
      #### The sliders and input
      sidebarPanel(
        sliderInput("range", label= "Range", min = 0, max = 5, value = 1,step = 0.01),
        sliderInput("sill", label = "Partial Sill", min = 0, max = 10, value = 1,step = 0.1),
        sliderInput("nugget", label = "Nugget", min = 0, max = 10, value = 0,step = 0.1),
        selectInput("model",label = "Model",choices = c('exponential','gaussian','spherical'),selected = 'exponential'),
        actionButton(inputId = "clicks", label = "Generate data")),
                              
                              
      #### The plots and output
      mainPanel(
        h3("This app illustrates the behaviour of a spatial correlation function."),
        plotOutput("geoex"),
        plotOutput('vario')  
                )
)
)

server <- function(input, output) {
#### Set up the simulation
set.seed(1234)
nx <- 50
ny <- 50

#### Generate the random field
dat <- eventReactive(input$clicks, {
  grf(n=(nx*ny), grid="reg", nx=nx, ny=ny, 
      cov.model=input$model, cov.pars = c(input$sill,input$range), nugget=input$nugget)
  })
  
#### Plot the field  
output$geoex <- renderPlot({
  image(dat(), x.leg=c(1.1,1.8), y.leg=c(0.4,0.5), xlim=c(-0.2, 2), axes = FALSE)
})


#### Plot the semi-variogram  
output$vario <- renderPlot({
  plot(variog(dat()),ylim=c(0,1.5 * (input$sill + input$nugget)), pch=19)
  lines.variomodel(dat(), col="red")
})

}


shinyApp(server = server, ui = ui)


