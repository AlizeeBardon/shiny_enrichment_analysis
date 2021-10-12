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
library(plotly)

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
                     fluidRow(
                         box(dataTableOutput("table")),
                     
                         headerPanel('Example'),
                         mainPanel(
                             plotlyOutput('plot')
                         ),
                      
                         hr(),
                         hr(),
                         #test volcanoplot 10 octobre  
                         titlePanel("Volcano Plot using ggplot "),
                         fluidRow(
                           column(
                             width = 6,
                             plotOutput("volcanoPlot", click = "volcanoPlotSelection", height = "500px")
                           )
                         ),
                         
                         hr(),
                         hr(),
                         hr()
                         
                     )#fluidRow(
                     
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