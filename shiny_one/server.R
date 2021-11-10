#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(highcharter)
library(DT)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    re <- reactive({
        req(input$file1)
        data <- read.csv(input$file1$datapath, header = TRUE, sep = ";") %>% #read.csv("NKI-DE-results.csv", stringsAsFactors = FALSE)  
          mutate(
            minusLog10Pvalue = -log10(pval),
            pvalue_subset = minusLog10Pvalue > -log10(input$color_pvalue)
            
            
          )
          })

    
    ###### volcanoplot using ploty TEST   
    
    output$volcanoPlot_plotly <- renderPlotly({
      chgm_axe_x <- input$axe_x
      chgm_axe_y <- input$axe_y
      pvalue_color <- input$color_pvalue
      data_plot_plotly <- re()
      plot_plotly <- 
        data_plot_plotly %>%
        ggplot(aes(x = log2FC,
                   y = minusLog10Pvalue,
                   colour = pvalue_subset,
                   text = GeneName )) +
        geom_point(size = 0.25) +
        
        xlab("log fold change") +
        ylab("-log10(P-value)") +
        xlim(-chgm_axe_x, chgm_axe_x) +
        ylim(0, chgm_axe_y) +
        guides(color = guide_legend(override.aes = list(size=12))) 
        
      
      plot_plotly %>%
        ggplotly(GeneName = "GeneName") %>%
        layout(dragmode = "select")
    })
    
    output$table <- renderDataTable({ # affiche tableau des data
        D <- re()
        pvalue_color <- input$color_pvalue
        D_subset = subset(D, pvalue_subset == TRUE)
        DT::datatable(D_subset)
    }) # fin renderDataTable({
    
    } # end function(input, output) {
) # end shinyServer(