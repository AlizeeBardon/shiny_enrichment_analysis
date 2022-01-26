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
                                ),
                  selected = NULL
                 
                ) , 
                
                

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
               
               box(title = "Parameters", status = "info", solidHeader = TRUE, width = 12,
                   
                   box (                   
                     sliderInput(inputId = "pvalue",
                                 label = "pvalue",
                                 min = 0,
                                 max = 0.25,
                                 value = 0.05),
                     width = 6),
                   box (
                     sliderInput(inputId = "tresholdLog2FoldChange",
                               label = "treshold for Log2FoldChange",
                               min = 0,
                               max = 4,
                               step = 0.1 ,
                               value = 0.4),
                     width = 6 )
                     
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
               
               box(title = "Parameters", status = "info", solidHeader = TRUE, width = 12,  
               
                   box(
                     radioButtons("radio_filtre_ontology", label = h4("Choose ontology"),
                                  choices = list(
                                    "all (BP + CC + MF)" = "all",
                                    "Biological Process (BP)"= "BP",
                                    "Cellular Component (CC)"= "CC",
                                    "Molecular Function (MF)"= "MF"
                                  ),
                                  selected = 1),
                     width = 4
                   
                     ),
                   
                   box(
                     h4("pvalue"),
                     verbatimTextOutput("pvalue_go_enrich"), 
                     h4("log2FoldChange threshold"),
                     verbatimTextOutput("log2foldchange_go_enrich"), 
                     h6("if you want to change the pvalue log2FC threshold, modified the corresponding slider in the Whole Data Analysis part"),
                     
                     width = 4
                   ), 
                   
                   box(
                   
                     radioButtons("subset_up_down_regulated", label = h4("Expression selection (upregulated or downregulated)"),
                                choices = list(
                                  "both up and down regulated" = "both",
                                  "up regulated (> log2foldchange)"= "up",
                                  "down regulated (< -log2foldchange)"= "down"
                                ),
                                
                                selected = 1),
                   width = 4
                   )
                   
               
                     
               ),
              
               
               br(),
                              
               box (dataTableOutput("Table_go_enrichment"), 
                    width = 12)
               
 
      ), # tabPanel("GO Term Enrichm:qent"
      
      
      tabPanel("Pathway Enrichment", 
               
               br(),
               
               h1("Pathway Enrichment"), 
               
               br(),
               

               
               
      ), # tabPanel("Pathway Enrichment"
      
      
      tabPanel("Protein Domain Enrichment", "text"
      ) # tabPanel("Protein Domain Enrichment"
      
      
    ) # fin tabsetPanel
    
  ) # fin dashboardBody(
  
  
) # fin dashboardPage(
) # fin shinyUI(
