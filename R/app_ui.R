#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinythemes
#' @import shinydashboard
#' @import tidyverse
#' @import plotly
#' @import highcharter
#' @import DT
#' @import shinyWidgets
#' @noRd


app_ui <- function(request) {
  
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
                  
                  
                  selectInput(
                    "espece",
                    label = h4("choose organism:"),
                    choices = list(
                      "Anophele" = "Anophele",
                      "Arabidopsis" = "Arabidopsis",
                      "Bovine" = "Bovine",
                      "Worm" = "Worm",
                      "Canine" = "Canine", 
                      "Fly" = "Fly",
                      "ZebraFish" = "ZebraFish",
                      "Ecoli K12" = "Ecoli K12",
                      "Ecoli Sakai" = "Ecoli Sakai",
                      "Chicken" = "Chicken",
                      "Human" = "Human",
                      "Mouse" = "Mouse",
                      "Rhesus" = "Rhesus",
                      "Myxococcus xanthus DK 1622" = "Myxococcus xanthus DK 1622",
                      "Malaria" = "Malaria" ,
                      "Chimpanzee" = "Chimpanzee",
                      "Rat" =  "Rat",
                      "Yeast" = "Yeast",
                      "Pig" = "Pig",
                      "Xenopus" = "Xenopus"
                    ),
                    selected = "Human"
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
        tags$link(rel = "stylesheet", type = "text/css", href = "custom_shinyapp.css")
      ),
      
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
                 
                 br(), br(),
                 
                 radioGroupButtons("method", label = h3("Analysis method"),
                                   choices = list("Over epresentation analysis (ORA)" = 1, "Gene Set Enrichment Analysis (GSEA)" = 2), direction = "horizontal"), 
                 radioGroupButtons("db", label = h3("DataBase"),
                                   choices = list("KEGG" = 1, "REACTOME (you can try, but it doesn't work...)" = 2), 
                                   selected = 1, checkIcon = list(
                                     yes = icon("ok", justufued = TRUE, lib = "glyphicon")), direction = "horizontal"),
                 radioGroupButtons("type", label = h3("DEG type:"), 
                                   choices = list("Over expressed DEG only" = 1, "Under expressed DEG only" = 2, "Both" = 3), 
                                   selected = 1, direction = "horizontal")
                 
                 
        ), # tabPanel("Pathway Enrichment"
        
        
        tabPanel("Protein Domain Enrichment", "text"
        ) # tabPanel("Protein Domain Enrichment"
        
        
      ) # fin tabsetPanel
      
    ) # fin dashboardBody(
    
    
  ) # fin dashboardPage(

    #fluidPage(
    #h1("testgolem"),
      

    #)
  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'testgolem'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

