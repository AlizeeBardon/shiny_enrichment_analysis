#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
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

# Define UI for application that draws a histogram
shinyUI(dashboardPage(
    
    dashboardHeader(title = "Enrichment Analysis"),
    
    dashboardSidebar(
        sidebarMenu(id="tabs",
                    fileInput("file1", "Choose CSV File",
                              accept = c(
                                  "text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")
                    ),
                    textInput("text", ""),
                    actionButton("printtext","RUN"),
                    menuItem("Data", tabName = "data", icon = icon("table"),startExpanded = TRUE,
                             menuSubItem("CSV file", tabName = "data1"),
                             menuSubItem("TSV file", tabName = "data")
                             
                    )
        )
    ), # fin dashboardSidebar
    
    
    dashboardBody(
        tabsetPanel(
            tabPanel("Whole Data Analysis",
                     
                     # fluidRow(
                     #     box(dataTableOutput("table")),
                     # ),
                     #     hr(),
                     #     hr(),
                         
                         ###### volcanoplot using ploty TEST   
                         
                          #titlePanel("Volcano Plotly"),
                         fluidRow(
                           column( #column1
                              width = 2,
                              sliderInput(inputId = "axe_x",
                                       label = "limite log2 fold change (axe x)",
                                       min = 0,
                                       max = 7,
                                       value = 3),
                              sliderInput(inputId = "axe_y",
                                          label = "limite -log10 p value (axe y)",
                                          min = 0,
                                          max = 60,
                                          value = 6),
                              
                              numericInput("color_pvalue", label = h3("pvalue limite selection"), value = 0.05, step = 0.01),
                           ), #fin column1
                              
                           column( #column2
                             width = 10,
                             plotlyOutput("volcanoPlot_plotly", height = "500px",width = "100%")
                           )#fin column2
                         ),#fluidRow(
                     
                     fluidRow(
                         box(dataTableOutput("table")),
                     )
                     
            ), # tabPanel("Whole Data Analysis"
            
            
            tabPanel("GO Term Enrichment", "text"
            ), # tabPanel("GO Term Enrichment"
            
            
            tabPanel("Pathway Enrichment", "text"
            ), # tabPanel("Pathway Enrichment"
        
            
            tabPanel("Protein Domain Enrichment", "text"
            ), # tabPanel("Protein Domain Enrichment"
            
            
            tabPanel("About", "text"
            ) # tabPanel("About"
            
            
            
            
            
            
        ) # fin tabsetPanel
        
    ) # fin dashboardBody(
    
    
) # fin dashboardPage(
) # fin shinyUI(