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
                               value = 0.1),
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
      
               br(), br(),
               
               box(
                 title = "Parameters - GO Term Enrichment", 
                 closable = TRUE, 
                 status = "primary",
                 width = 12,
                 solidHeader = FALSE, 
                 collapsible = TRUE
               )
               
      ),

      # BODY: tabPanel : Pathway Enrichment -------------------------------- 
      
      tabPanel("Pathway Enrichment", 
               
               br(), br(),
               
               box(
                   title = "Parameters - Pathway Enrichment", 
                   closable = TRUE, 
                   status = "primary",
                   width = 12,
                   solidHeader = FALSE, 
                   collapsible = TRUE,
                   
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


      
      # BODY: tabPanel : Protein Domain Enrichment -------------------------------- 

      tabPanel("Protein Domain Enrichment",
               
               br(),
               
               br(),
               
               box(
                 title = tags$b("Parameters - Protein Domain Enrichment"), 
                 closable = TRUE, 
                 status = "primary",
                 width = 12,
                 solidHeader = FALSE, 
                 collapsible = TRUE,
                 

                 br(),
                   
                 column(4, 
                        
                        radioButtons("method_prt_domain", 
                                     label = h4(
                                       "Analysis method",
                                       tipify(actionBttn(
                                         inputId = "prt5",
                                         icon = icon("question"),
                                         style = "jelly",
                                         size = "xs",
                                         block = FALSE,
                                         no_outline = TRUE
                                       ), "ORA : hypergeom test with BH adjustment", placement="bottom", trigger = "hover") 
                                     ),
                                     choices = list(
                                       "Over epresentation analysis (ORA)" = 1, 
                                       "Gene Set Enrichment Analysis (GSEA)" = 2), 
                                     selected = 1),
                        br(),
                        radioButtons("type_prt_domain", 
                                     label = h4("DEG type",
                                                tipify(actionBttn(
                                                  inputId = "prt5",
                                                  icon = icon("question"),
                                                  style = "jelly",
                                                  size = "xs",
                                                  block = FALSE,
                                                  no_outline = TRUE
                                                ), "Diffrential Expression Group", placement="bottom", trigger = "hover") 
                                                ),
                                     choices = list(
                                       "Both" = "Both",
                                       "Over expressed DEG only" = "Over", 
                                       "Under expressed DEG only" = "Under" 
                                       ),
                                     selected = NULL)
                 ),
                 column(4,
                        textInput("biomart_listMarts", 
                                  label = h4(
                                    "BioMart database ",
                                    tipify(actionBttn(
                                      inputId = "prt1",
                                      icon = icon("question"),
                                      style = "jelly",
                                      size = "xs",
                                      block = FALSE,
                                      no_outline = TRUE
                                    ), "This must be BioMart databases to which biomaRt can connect to (cf listMarts)", placement="bottom", trigger = "hover") 
                                    ), value = "ensembl"),
                        br(),br(),
                      
                               
                        textInput("biomart_dataset", 
                                  label = h4(
                                    "BioMart Dataset",
                                    tipify(actionBttn(
                                      inputId = "prt2",
                                      icon = icon("question"),
                                      style = "jelly",
                                      size = "xs",
                                      block = FALSE,
                                      no_outline = TRUE
                                    ), "Enter a valid BioMart Dataset for the BioMart database (for example, within the Ensembl genes mart every species is a different dataset)", placement="bottom", trigger = "hover") 
                                    ), value = "mmusculus_gene_ensembl" )
                 ),
                 
                 column(4,
                        sliderInput(inputId = "pvalue_prt_domain",
                                    label = h4(
                                      "Adjusted p-value cutoff",
                                      tipify(actionBttn(
                                        inputId = "prt3",
                                        icon = icon("question"),
                                        style = "jelly",
                                        size = "xs",
                                        block = FALSE,
                                        no_outline = TRUE
                                      ), "Output parameters for graph creation (thresholding of domains on their adjusted p-value", placement="bottom", trigger = "hover") 
                                      ),
                                    min = 0,
                                    max = 1,
                                    value = 0.05),
                        
                        selectInput( inputId = "pvalue_adjustment_prt_domain",
                                     label = h4("Adjustment method:"),
                                     choices = list(
                                       'Holm (1979) ("holm")' = "holm", 
                                       'Hochberg (1988) ("hochberg")' = "hochberg", 
                                       'Hommel (1988) ("hommel")' = "hommel", 
                                       'Bonferroni correction ("bonferroni")' = "bonferroni", 
                                       'Benjamini & Hochberg (1995) ("BH" or its alias "fdr")' = "BH", 
                                       'Benjamini & Yekutieli (2001) ("BY")' = "BY",
                                       'none' = "none"
                                     ),
                                     selected = NULL
                        )
                        
                        
                 ),

                 
                 box (
                   status = "primary",
                   actionButton("Run_protein_domains",h4("Run Protein Domains")),
                   width = 12)
               ),
               
               
               br(), br(), 
               
               
               conditionalPanel(
                 condition = "input.method_prt_domain == 1",
                
                 h1(strong("Over-representation (or enrichment) analysis"),  br(),
                    tags$a(href="https://github.com/AlizeeBardon/shiny_enrichment_analysis","coded method"),
                    align = "center"),

                 br(), 
                 h2("Protein domain enrichment - Result table"),
                 dataTableOutput("Table_domains_enrichment"), 
                 br(), br(),
                 
                 
                 h2("Protein domain enrichment - BarPlot"),
                 sliderInput(inputId = "nb_barplot_ora_coder",
                             label = h4("number of protein to see in the barplot :"  ),
                             min = 0,
                             max = 200,
                             value = 10),
                 
                 plotlyOutput("barplot_domains_enrichment"),
                 br(), br(),
                 
                 h1( strong("Over-representation (or enrichment) analysis"), 
                     br(),
                     tags$a(href="https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/enricher","enricher() function from clusterProfiler"),
                     align = "center"),
                 br(),
                 h2("Protein domain enrichment - Result table"),
                 dataTableOutput("Table_domains_enrichment_enricher"), 
                 br(),
                 h2("Protein domain enrichment - BarPlot"),
                 

                 plotlyOutput("barplot_domains_enrichment_enricher"),
                 br(),
                 h2("Protein domain enrichment - DotPlot"),
                 plotlyOutput("dotplot_domains_enrichment_enricher")
                 
               ),
               
               
               conditionalPanel(
                 condition = "input.method_prt_domain == 2",
                 h1(strong("Gene Set Enrichment Analysis (GSEA) - using the method enricher() from clusterProfiler method "), 
                    br(), 
                    tags$a(href="https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/GSEA","GSEA() function from clusterProfiler"),
                    align = "center"),

                 
                 br(), br(),
                 dataTableOutput("Table_domain_enrichment_GSEA"), 
                 br(), br(),
                 
                 "Barplot :",
                 plotlyOutput("barplot_domain_enrichment_GSEA"),
                 
                 br(), br(),

                     
                 selectInput(
                       "protein_id_list", 
                       label = h4("choose the protein id to display"),
                       choices = "",
                       selected = NULL),

                 
                 br(), br(),
                 plotlyOutput("gseaplot_domain_enrichment_GSEA")
               )  
               

               
                   
      ) # tabPanel("Protein Domain Enrichment"
      
      
    ) # fin tabsetPanel
    
  ) # fin dashboardBody(
  
  
) # fin dashboardPage(
) # fin shinyUI(