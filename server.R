#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt", "pathview", "clusterProfiler")


organism = "org.Mm.eg.db" 

library(clusterProfiler) 
library(biomaRt)
library(pathview)
library(shiny)
library(shinythemes)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(highcharter)
library(DT)
library(ggplot2)



# Define server 
shinyServer(function(input, output) {

######################################################################### 
###### input data    
#########################################################################
  
re <- reactive({
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      req(file)
      print(head(file))
      validate(need(ext == "csv", "Invalid file. Please upload a .csv file"))
      data <- read.csv(file$datapath, header = TRUE, sep = ";") %>% 
        mutate(
          minusLog10Pvalue = -log10(padj)
          )
        })


    # BODY --------------------------------------------------------------------
    
        # BODY: tabPanel : whole Data Analysis --------------------------------


# Volcanoplot 


    output$volcanoPlot_plotly <- renderPlotly({

        data_plot_plotly <- re()
        data_plot_plotly$key <- row.names(data_plot_plotly)
        data_plot_plotly$col <- "black"
        
        click_data <- event_data("plotly_click")
        select_data <- event_data("plotly_selected")

        pvalue <- -log10(input$pvalue)
          
        data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC > 0, "col"] <- "green"
        data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC < 0, "col"] <- "red"
        
        p <- ggplot(data = data_plot_plotly, mapping = aes(x = log2FC, y = minusLog10Pvalue, col = I(col), key = key)) +
          geom_point(size = 0.45) + 
          theme_bw()  +
          xlab("log2 fold change") + ylab("-log10(p-value)")
        
        event_register(p, 'plotly_click')
        event_register(p, 'plotly_selected')
        
        plotly_object <- ggplotly(p,source = "source1") %>% 
          layout(dragmode = "lasso")  %>% 
          layout(showlegend = FALSE)
        
      })
    

# MA plot    

    output$MAPlot_plotly <- renderPlotly({
      
      data_plot_plotly <- re()
      data_plot_plotly$key <- row.names(data_plot_plotly)
      data_plot_plotly$col <- "black"
      
      click_data <- event_data("plotly_click")
      select_data <- event_data("plotly_selected")
      
      pvalue <- -log10(input$pvalue)
      
      data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC > 0, "col"] <- "green"
      data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC < 0, "col"] <- "red"
      
      p <- ggplot(data = data_plot_plotly, mapping = aes(
        x = log2(baseMean), y = log2FC, col = I(col), key = key)) +
        geom_point(size = 0.45) + 
        theme_bw()  +
        xlab("log2(baseMean)") + 
        ylab("log2FC")
      
      plotly_object <- ggplotly(p,source = "source1") %>% 
        layout(dragmode = "lasso")  %>% 
        layout(showlegend = FALSE)
      
    })
    

# Interactive Table

    output$Table_subset_data_selected <- renderDataTable({ 
      D <- re() %>%  
        mutate(
          indice = row_number()
        )
      data_selected_from_graph <- event_data("plotly_selected",source = "source1")

      D_subset = subset(D, indice %in%  data_selected_from_graph$key)
      DT::datatable(D_subset) 
    }) # fin renderDataTable({

    
    proxy <- DT::dataTableProxy("Table_subset_data_selected")
    
    observe({
      data <- re()
      subset_data <- data.frame(
        ID_test = data$ID,
        row_id = 1:length(data$ID),
        stringsAsFactors = FALSE
      )
      s <- event_data("plotly_click",source = "source1")
      req(!is.null(s))
      # map point number to subset_data
      row_clicked <- subset_data[s$pointNumber + 1,"row_id"]
      proxy %>%
        selectRows(NULL) %>%
        selectRows(row_clicked)
      
    }) # fin observe
    
# Annotation Table

    annot <- eventReactive(input$Run_Annotation, {
        data <- re() 
        gene_list = data$GeneName                 
        organism = input$espece
        generef = bitr(gene_list, fromType = "SYMBOL", toType= "GO", OrgDb=organism)
    })

    output$annotation <- renderDataTable({
      D <- annot()
      DT::datatable(D)
    })#fin renderDataTable
    
    
    } # end function(input, output) {
) # end shinyServer(