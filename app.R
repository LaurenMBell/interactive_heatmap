# THE R PART OF THE HEATMAP VISUALIZATION
library(shiny) 
library(reticulate)
library(bslib)

use_condaenv("reticulate", required=TRUE)

source_python("shiny.py")

subnetworks <- c("C0_U", "C0_D", "C1_U", "C1_D", "C2_U", "C2_D", "C3_U", "C3_D", 
                 "C4_U", "C4_D", "C5_U", "C5_D", 
                 "P0_U", "P0_D", "P1_U", "P1_D", "P2_U", "P2_D", "P3_U", "P3_D", "P4_U", "P4_D", 
                 "S0_U", "S0_D","S1_U", "S1_D", "S2_U", "S2_D", "S3_U", "S3_D", "S4_U", "S4_D")

ui <- page_fluid(
  titlePanel("Subnetwork Heatmap Visualization"), 
  title = "heatmaps", 
  layout_sidebar(
    title = "Controls",
    sidebar = sidebar(
      sliderInput( 
          "sig_slider", "Significance Range (nLogP)", 
          min = 1.3, max = 5, 
          value = c(1.3, 5) 
                  ), 
      radioButtons( 
        "GO_level_select", 
        "GO Terms Shown", 
        choices = list("Summary GO Terms" = 1, "All GO Terms" = 2),
        selected = 1
      )
    )
  )
)
  
  card(
    card_header("Heatmap Visualization From Selected"),
    card_body(
      plotOutput("heatmap")
    )
  )



server <- function(input, output) {
  output$heatmap <- renderPlot({
    # CALL MAIN FROM THE PYTHON SCRIPT HERE 
    print(input$GO_level_select)
    
    startup(input$GO_level_select, input$sig_slider[1], input$sig_slider[2])
  })
}

shinyApp(ui, server)



#checkboxGroupInput('subnetwork_select',
                  #'Select Subnetworks\n(Will not work yet)',
                  #subnetworks
                  #), 


