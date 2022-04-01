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
library(shinyBS)
library(ggridges)
#library(upsetplot)

# Define UI for application that draws a histogram
shinyUI(dashboardPage(
  skin = "purple",
  title = "Enrichment Analysis",
  
  # HEADER ------------------------------------------------------------------
  
  dashboardHeader(
    title = span(img(src = "cute-daisies.png", height = 25), "Enrichment Analysis"),
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
                fileInput("file1", h4("Choose CSV File"),
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
                  label = h4("Choose organism:"),
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
      
      # BODY: tabPanel : Whole Data Analysis --------------------------------
      
      tabPanel("Whole Data Inspection",
            
               br(), br(),
               
               box(
                 title = "Parameters - Whole Data Inspection", 
                 closable = TRUE, 
                 status = "primary",
                 width = 12,
                 solidHeader = FALSE, 
                 collapsible = TRUE,
                 box (                   
                   sliderInput(inputId = "pvalue",
                               label = "pvalue",
                               min = 0,
                               max = 1,
                               value = 0.05),
                   width = 6),
                 box (
                   sliderInput(inputId = "tresholdLog2FoldChange",
                               label = "treshold for Log2FoldChange",
                               min = 0,
                               max = 5,
                               step = 0.1 ,
                               value = 1),
                   width = 6 )
                 
                 
               ),
               
               
               box(
                 title = "Volcano Plot",
                 status = "primary",
                 plotlyOutput("volcanoPlot_plotly"),
                 width = 6
               ), #fin box
               
               box(
                 title = "MA Plot",
                 status = "primary",
                 plotlyOutput("MAPlot_plotly"),
                 width = 6
               ), #fin box
               
               
               
               box (
                 title = "Selected data",
                 status = "primary",
                 dataTableOutput("Table_subset_data_selected"), 
                 downloadButton('download',"Download"),
                 width = 12)
            
             
      ), # tabPanel("Whole Data Analysis"
    
      
      # BODY: tabPanel : GO Term Enrichment --------------------------------
      
      tabPanel("GO Term Enrichment", 
               
               br(),
               
               h1("GO Term Enrichment"), 
               
               br(),
               
               
               box(
                 selectInput( "Ontology",
                              label = h4("choose ontology:"),
                              choices = c(
                                "BP"="BP",
                                "CC"="CC",
                                "MF"="MF",
                                "all"=""
                              ),
                              selected = "BP"),
                 
                 width = 4
               ),
               box(
                 selectInput( "Ajustement",
                              label = h4("choose ajustement method:"),
                              choices = c("holm", "hochberg", 
                                          "hommel", "bonferroni", 
                                          "BH", "BY", "fdr", 
                                          "none"), 
                              selected = "none"),
                 width = 4
                 
                 
                 
                 
               ),#fin box
               
               box(
                 numericInput("showCategory_enrichmap", "number of categories to show", value = 5),
                 width = 4
                 
               ),#fin box
               
               box(title = "Dot Plot GSEA", status = "warning", solidHeader = TRUE, width = 12, height = "600px",
                   
                   plotlyOutput("dotplot",  height = "500px")
               ),#fin box
               
               box(title = "ridge Plot GSEA", status = "warning", solidHeader = TRUE, width = 12, height = "600px",
                   
                   plotlyOutput("ridgeplot",  height = "500px")
               ),#fin box
               
               box(title = "gsea Plot", status = "warning", solidHeader = TRUE, width = 12, height = "1000px",
                   
                   plotlyOutput("gsea_plot",  height = "900px")
               ),#fin box
               
               box (dataTableOutput("goGse_annot_table"), width = 12, style = "overflow-x: scroll;"),
               
               box(title = "Bar Plot SEA", status = "warning", solidHeader = TRUE, width = 12, height = "600px",
                   
                   plotlyOutput("barplot",  height = "500px")
               ),#fin box
               box(title = "Dot Plot SEA", status = "warning", solidHeader = TRUE, width = 12, height = "600px",
                   
                   plotlyOutput("dotplot_sea",  height = "500px")
               ),#fin box
               box(title = "upsetplot SEA", status = "warning", solidHeader = TRUE, width = 12, height = "600px",
                   
                   plotlyOutput("usetplot",  height = "500px")
               ),#fin box
               
               box(title = "goplot SEA", status = "warning", solidHeader = TRUE, width = 12, height = "600px",
                   
                   plotlyOutput("goplot",  height = "500px")
               ),#fin box
               box (dataTableOutput("goGse_enrich_table"), width = 12, style = "overflow-x: scroll;"),
               
               
 
      ), # tabPanel("GO Term Enrichm:qent"
      

      # BODY: tabPanel : Pathway Enrichment -------------------------------- 
      
      tabPanel("Pathway Enrichment", 
               
               
               br(),
               
               h1("Pathway Enrichment"), 
               
               br(),
               
               
               box(title = "Parameters", status = "info", solidHeader = TRUE, width = 12,  
                   
                   box(
                     radioGroupButtons("method", label = h3("Analysis method"),
                                       choices = list(
                                         "Over epresentation analysis (ORA)" = 1, 
                                         "Gene Set Enrichment Analysis (GSEA)" = 2), 
                                       direction = "vertical"),
                     width = 4
                   ),
                   
                   box(
                     
                     radioGroupButtons("db", label = h3("DataBase"),
                                       choices = list(
                                         "KEGG" = 1, 
                                         "REACTOME (you can try, but it doesn't work...)" = 2), 
                                       direction = "vertical"),
                     width = 5
                   ),
                   
                   box(
                     radioGroupButtons("type", label = h3("DEG type:"), 
                                       choices = list(
                                         "Over expressed DEG only" = 1, 
                                         "Under expressed DEG only" = 2, 
                                         "Both" = 3), 
                                       direction = "vertical"),
                     width = 3
                   ),
                   
                   box(
  
                   sliderInput(inputId = "pvalue_gsea",
                               label = "pvalue cutoff for enrichment analysis",
                               min = 0,
                               max = 0.25,
                               value = 0.05),
                   width = 9
                   
                   ),
                   
                   box(
                     actionButton("Run_Annotation_ENSEMBL_to_GO","Run Annotation"),
                     width = 3
                   )
                   
                   
               ),
               
               box (dataTableOutput("enrichKEGG_table"), 
                    width = 12, style = "overflow-x: scroll;"),
               
               box(tstatus = "warning", solidHeader = TRUE, width = 12, height = "550px",
                   plotlyOutput("dotplot_kegg", height = "450px")
               ),
               
               
               mainPanel(title = "GSE Plot", status = "warning", solidHeader = TRUE, width = 12, height = "550px",
                   plotlyOutput("method_kegg", height = "450px")
               ),
              
               sidebarLayout(
                 sidebarPanel(selectInput("paths", label = h4("Choose a pathway"),choices = "",selected = NULL),actionButton("go", "Generate pathway with pathview")),
                 mainPanel(
                   bsModal("modalExample", "KEGG PATHWAY", "go", imageOutput("pathview_kegg"))
                 )
               )
               ), # tabPanel("Pathway Enrichment"


      
      

      tabPanel("Protein Domain Enrichment",
                
               br(),
               
               h1("Protein Domain Enrichment"), 
               
               br(),
               
               
               box(title = "Parameters", status = "info", solidHeader = TRUE, width = 12,
                   box(
                     radioGroupButtons("method", label = h4("Analysis method"),
                                       choices = list(
                                         "Over epresentation analysis (ORA)" = 1, 
                                         "Gene Set Enrichment Analysis (GSEA)" = 2), 
                                       direction = "vertical"),
                     width = 4
                   ),
                   
                   box (                   
                     sliderInput(inputId = "pvalue_domains",
                                 label = h4("Select a adjusted p-value cutoff"),
                                 min = 0,
                                 max = 0.25,
                                 value = 0.05),
                     width = 4),
                   
                   box(
                     radioGroupButtons("type", label = h4("DEG type:"), 
                                       choices = list(
                                         "Over expressed DEG only" = 1, 
                                         "Under expressed DEG only" = 2, 
                                         "Both" = 3), 
                                       direction = "vertical"),
                     width = 4
                   ),
                   
                   box(

                     textInput("biomart_listMarts", label = h4("biomart_listMarts"), value = "ensembl"),
                     width = 4
                   ),
                   
                   box(
                     textInput("biomart_dataset", label = h4("biomart_dataset"), value = "mmusculus_gene_ensembl" ),
                     width = 4
                   ),
                   
                   box(
                     actionButton("Run_protein_domains",h4("Run Protein Domains")),
                     width = 4
                   )
               ),
               
               
               br(),
               
               box (dataTableOutput("Table_domains_enrichment"), 
                    width = 12),
               
               box(
                 title = "Bar Plot", status = "warning", solidHeader = TRUE, width = 6, height = "550px",
                 plotlyOutput("barplot_domains_enrichment", height = "450px")
               ),
                   
      ) # tabPanel("Protein Domain Enrichment"
      
      
    ) # fin tabsetPanel
    
  ) # fin dashboardBody(
  
  
) # fin dashboardPage(
) # fin shinyUI(