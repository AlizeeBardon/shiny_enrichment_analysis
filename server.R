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
            #probe.type = factor(ifelse(grepl(>7.12786E-14, padj), "sup", "inf")),
            minusLog10Pvalue = -log10(pval),
            #tooltip = ifelse(is.na(HUGO.gene.symbol), probe, paste(HUGO.gene.symbol, " (", probe, ")", sep = ""))
          )
          })

    # output$table <- renderDataTable({ # affiche tableau des data
    #     D <- re()
    #     DT::datatable(D)
    # }) # fin renderDataTable({

    
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
                   colour = minusLog10Pvalue < pvalue_color,
                   text = GeneName )) +
        geom_point(size = 0.25) +
        xlab("log fold change") +
        ylab("-log10(P-value)") +
        xlim(-chgm_axe_x, chgm_axe_x) +
        ylim(0, chgm_axe_y)
      
      plot_plotly %>%
        ggplotly(GeneName = "GeneName") %>%
        layout(dragmode = "select")
    })
    
    output$table <- renderDataTable({ # affiche tableau des data
        D <- re()
        pvalue_color <- input$color_pvalue
        D_subset = subset(D, minusLog10Pvalue >= pvalue_color)
        DT::datatable(D_subset)
    }) # fin renderDataTable({
    
    } # end function(input, output) {
) # end shinyServer(