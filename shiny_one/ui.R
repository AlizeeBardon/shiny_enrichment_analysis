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
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(highcharter)
library(DT)

# Define UI for application that draws a histogram
shinyUI(dashboardPage(
  skin = "purple",
  title = "Enrichment Analysis",
  
  # HEADER ------------------------------------------------------------------
  
  dashboardHeader(
    title = span(img(src = "DNA.png", height = 25), "Enrichment Analysis"),
    titleWidth = 250,
    dropdownMenu(
      type = "notifications", 
      headerText = strong("HELP"), 
      icon = icon("question"), 
      badgeStatus = NULL,
      
      
      notificationItem(
        text = ("plus d informations ?"),
        icon = icon("spinner")
      ),
      notificationItem(
        text = "nous joindre",
        icon = icon("address-card")
      ),
      notificationItem(
        text = "nous joindre",
        icon = icon("calendar")
      ),
      notificationItem(
        text = "nous joindre",
        icon = icon("user-md")
      ),
      notificationItem(
        text = "nous joindre",
        icon = icon("ambulance")
      ),
      notificationItem(
        text = "nous joindre",
        icon = icon("flask")
      ),
      notificationItem(
        text = strong("attention!"),
        icon = icon("exclamation")
      )
    ), # dropdownMenu
    
    tags$li(
      a(
        strong("ABOUT US"),
        height = 40,
        href = "https://github.com/ceefluz/radar/blob/master/README.md",
        title = "",
        target = "_blank"
      ),
      class = "dropdown"
    )
  ),  # fin dashboardHeader
  
  
  # SIDEBAR -----------------------------------------------------------------
  
  dashboardSidebar(
    width = 250,
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
  
  
  # BODY --------------------------------------------------------------------
  
  
  dashboardBody(

    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "css_shinyapp_first_test.css")
    ),
    
    tags$style(
    
      HTML("
    .tabbable > .nav > li > a                  {background-color: transparent;  color:black; font-size: 12px}
    .tabbable > .nav > li[class=active]    > a {background-color: purple; color:white; font-size: 12px}
      ")#html
      
    ), # tags$style
    
    tabsetPanel(
      tabPanel("Whole Data Analysis",
            
               br(),
               br(),
               
               box(title = "Parameters", status = "warning", solidHeader = TRUE, width = 12,
                   sliderInput(inputId = "pvalue",
                               label = "pvalue",
                               min = 0,
                               max = 0.25,
                               value = 0.05),
               ), # fin box
               
               
               box(title = "Volcano Plot", status = "warning", solidHeader = TRUE, width = 6, height = "550px",
                    plotlyOutput("volcanoPlot_plotly", height = "450px")
               ), #fin box
               
               box(title = "MA Plot", status = "warning", solidHeader = TRUE, width = 6, height = "550px",
                   plotlyOutput("MAPlot_plotly",  height = "450px")
               ), #fin box
               
  

               
               br(),
               br(),               br(),
               br(),               br(),
               br(),               br(),
               br(),
           
                 
 
               box (dataTableOutput("Table_subset_data_selected"), 
                    width = 12),
               
               
            
             
      ), # tabPanel("Whole Data Analysis"
      
      
      tabPanel("GO Term Enrichment", 
               
               
               
               
               "text",
               fluidRow(
                        tabBox(
                          title = "First tabBox",
                          # The id lets us use input$tabset1 on the server to find the current tab
                          id = "tabset1", height = "250px",
                          tabPanel("Tab1", "First tab content"),
                          tabPanel("Tab2", "Tab content 2")
                        )
               ) #fluidrow
 
      ), # tabPanel("GO Term Enrichment"
      
      
      tabPanel("Pathway Enrichment", 
               
               "text"
               
      ), # tabPanel("Pathway Enrichment"
      
      
      tabPanel("Protein Domain Enrichment", "text"
      ) # tabPanel("Protein Domain Enrichment"
      
      
    ) # fin tabsetPanel
    
  ) # fin dashboardBody(
  
  
) # fin dashboardPage(
) # fin shinyUI(