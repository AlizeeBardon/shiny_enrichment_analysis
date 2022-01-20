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
        href = "https://github.com/AlizeeBardon/shiny_enrichment_analysis",
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

                menuItem(
                  "need any help to import your data ?",
                  h5("please choose a .csv document"),
                  h5("the document must be composed ", br(), " of 6 columns(GeneName, ID, baseMean, ", br(), " log2FC, pval, padj)")
                ),
                
                selectInput( "espece",
                  label = h4("choose organism:"),
                  choices = list(
                                  "Human (org.Hs.eg.db)"="org.Hs.eg.db",
                                  "Mouse (org.Mm.eg.db)"="org.Mm.eg.db",
                                  "Rat (org.Rn.eg.db)"="org.Rn.eg.db",
                                  "Yeast (org.Sc.sgd.db)"="org.Sc.sgd.db",
                                  "Fly (org.Dm.eg.db)"="org.Dm.eg.db",
                                  "Arabidopsis (org.At.tair.db)"="org.At.tair.db",
                                  "Zebrafish (org.Dr.eg.db)"="org.Dr.eg.db",
                                  "Bovine (org.Bt.eg.db)"="org.Bt.eg.db",
                                  "Worm (org.Ce.eg.db)"="org.Ce.eg.db",
                                  "Chicken (org.Gg.eg.db)"="org.Gg.eg.db",
                                  "Canine (org.Cf.eg.db)"="org.Cf.eg.db",
                                  "Pig (org.Ss.eg.db)"="org.Ss.eg.db",
                                  "Rhesus (org.Mmu.eg.db)"="org.Mmu.eg.db",
                                  "E coli strain K12 (org.EcK12.eg.db)"="org.EcK12.eg.db",
                                  "Xenopus (org.Xl.eg.db)"="org.Xl.eg.db",
                                  "Chimp (org.Pt.eg.db)"="org.Pt.eg.db",
                                  "Anopheles (org.Ag.eg.db)"="org.Ag.eg.db",
                                  "Malaria (org.Pf.plasmo.db)"="org.Pf.plasmo.db",
                                  "E coli strain Sakai (org.EcSakai.eg.db)"="org.EcSakai.eg.db"
                                )
                 
                ) , 
                
                actionButton("Run_Annotation","Run Annotation"),

                br(),
                
                menuItem(
                  "Annotation Information ?",
                  h5("Make sure you select the ", br(), "species that matches  your dataset"),
                  h5("If it is not in the list, we advise", br(), " to choose the most phylogenetically ", br(), "related species")


                ) #menuItem

                
    )
  ), # fin dashboardSidebar
  
  
  # BODY --------------------------------------------------------------------
  
  
  dashboardBody(

    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom_shinyapp.css")
    ),
    
    
    
    tabsetPanel(
      
      # BODY: tabPanel : whole Data Analysis --------------------------------
      
      tabPanel("Whole Data Analysis",
            
               br(),

               h1("whole Data Analysis"), 
               
               br(),
               
               box(title = "Parameters", status = "warning", solidHeader = TRUE, width = 12,
                   sliderInput(inputId = "pvalue",
                               label = "pvalue",
                               min = 0,
                               max = 0.25,
                               value = 0.05)
               ), # fin box
               
               
               box(title = "Volcano Plot", status = "warning", solidHeader = TRUE, width = 6, height = "550px",
                    plotlyOutput("volcanoPlot_plotly", height = "450px")
               ), #fin box
               
               box(title = "MA Plot", status = "warning", solidHeader = TRUE, width = 6, height = "550px",
                   plotlyOutput("MAPlot_plotly",  height = "450px")
               ), #fin box
               
 
               box (dataTableOutput("Table_subset_data_selected"), 
                    width = 12),
               
               box (dataTableOutput("annotation"), 
                    width = 12)
            
             
      ), # tabPanel("Whole Data Analysis"
    
      
      # BODY: tabPanel : GO Term Enrichment --------------------------------
      
      tabPanel("GO Term Enrichment", 
               
               br(),
               
               h1("GO Term Enrichment"), 
               
               br(),
               
               
               selectInput( "filtre_annotation",
                            label = h4("choose organism:"),
                            choices = list(
                              "BP"="BP",
                              "CC"="CC",
                              "MF"="MF",
                              "all"=""
                            )
               ),
              
                            
               br(),
                              
               box (dataTableOutput("Table_go_enrichment"), 
                    width = 12),
               
               br(),
               

	       fluidRow(
                        tabBox(
                          title = "First tabBox",
                          # The id lets us use input$tabset1 on the server to find the current tab
                          id = "tabset1", height = "250px",
                          tabPanel("Tab1", "First tab content"),
                          tabPanel("Tab2", "Tab content 2")
                        )
               ) #fluidrow
 
      ), # tabPanel("GO Term Enrichm:qent"
      

      # BODY: tabPanel : Pathway Enrichment -------------------------------- 
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
                                 selected = 1, direction = "horizontal"),
               box (dataTableOutput("Table_kegg"),
                    width = 12),
               box(title = "Dot Plot", status = "warning", solidHeader = TRUE, width = 6, height = "550px",
                   plotlyOutput("dotplot_kegg", height = "450px")
               ),
               box(title = "Bar Plot", status = "warning", solidHeader = TRUE, width = 6, height = "550px",
                   plotlyOutput("barplot_kegg", height = "450px")
               )
               #fluidRow(DT::dataTableOutput('table_kegg'))
               
               # div(
               # box(title = "Parametres", status = "warning", solidHeader = TRUE, width = 4,
               #     sliderInput(inputId = "axe_x",
               #                 label = "log2 fold change (axe x)",
               #                 min = 0,
               #                 max = 7,
               #                 value = 3),
               #     br(),
               #     sliderInput(inputId = "axe_y",
               #                 label = "-log10(p value) (axe y)",
               #                 min = 0,
               #                 max = 60,
               #                 value = 6),
               #     br(),
               #     numericInput("color_pvalue", label = h3("pvalue limite selection"), value = 0.05, step = 0.01),
               #     
               # ) # fin box
               
               
      ), # tabPanel("Pathway Enrichment"
      
      
      tabPanel("Protein Domain Enrichment", "text"
      ) # tabPanel("Protein Domain Enrichment"
      
      
    ) # fin tabsetPanel
    
  ) # fin dashboardBody(
  
  
) # fin dashboardPage(
) # fin shinyUI(
