#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    re <- reactive({
        req(input$file1)
        data <- read.csv(input$file1$datapath, header = TRUE, sep = ";")
        #print(head(data))
    })
    

    output$table <- renderDataTable({ # affiche tableau des data
        D <- re()
        DT::datatable(D)
    }) #renderDataTable({
    
###### test volcanoplot 10 octobre
    
    differentialExpressionResults <-
    
    output$volcanoPlot <- renderPlot({
      data_plot <- re()
      #probe.type = factor(ifelse(grepl("^Contig", probe), "EST", "mRNA"))
      df <- data.frame(
        log2FC = data_plot$log2FC,
        minusLog10Pvalue = -log10(data_plot$padj)
      )
      print(head(df))
      ggplot(df, aes(x = log2FC, y = minusLog10Pvalue)) +  #, colour = probe.type
        geom_point() +
        xlab("log fold change") +
        ylab("-log10(P-value)") +
        theme(legend.position = "bottom")
    })
    
    
    
    
    } # end function(input, output) {
) # end shinyServer(