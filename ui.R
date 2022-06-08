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
library(shinycustomloader)
library(viridis)

# Define UI for application that draws a histogram
shinyUI(dashboardPage(
  skin = "green",
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
                fileInput("file1", 
                          label = h4( 
                                     tipify(actionLink(
                                       inputId = "choose_csv_file",
                                       label = "Choose CSV File" ,
                                       style = "simple",
                                       size = "md",
                                       block = FALSE,
                                       no_outline = TRUE
                                     ), "please choose a .csv document the document must be composed  of 5 columns (ID, baseMean, log2FC, pval, padj)", 
                                     placement="bottom", 
                                     trigger = "hover"), 
                                     ),
                          
                          accept = c(
                            "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
                ),
               
                menuItem(
                 "Input exemple",
                 icon = icon("table"),
                 br(),
                 img(src = "cadre_exemple_input.png", width = 400),
                 br(),br(),
                 downloadButton('download_exemple',  "Exemple", status = "primary",)
                 
                ),
                
                br(),br(),
                img(src = "banner_daisy.png", width = "240 e"),
                
                selectInput( "espece_id",
                             label = h4("Choose organism:"),
                             choices = list(
                               "Human (org.Hs.eg.db)"=1,
                               "Mouse (org.Mm.eg.db)"=2,
                               "Rat (org.Rn.eg.db)"=3,
                               "Yeast (org.Sc.sgd.db)"=4,
                               "Fly (org.Dm.eg.db)"=5,
                               "Arabidopsis (org.At.tair.db)"=6,
                               "Zebrafish (org.Dr.eg.db)"=7,
                               "Bovine (org.Bt.eg.db)"=8,
                               "Worm (org.Ce.eg.db)"=9,
                               "Chicken (org.Gg.eg.db)"=10,
                               "Canine (org.Cf.eg.db)"=11,
                               "Pig (org.Ss.eg.db)"=12,
                               "Rhesus (org.Mmu.eg.db)"=13,
                               "E coli strain K12 (org.EcK12.eg.db)"=14,
                               "Xenopus (org.Xl.eg.db)"=15,
                               "Chimp (org.Pt.eg.db)"=16,
                               "Anopheles (org.Ag.eg.db)"=17,
                               "Malaria (org.Pf.plasmo.db)"=18,
                               "E coli strain Sakai (org.EcSakai.eg.db)"=19
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
                               label = "Adjusted p-value (padj) cutoff",
                               min = 0,
                               max = 1,
                               value = 0.05),
                   width = 6),
                 box (
                   sliderInput(inputId = "tresholdLog2FoldChange",
                               label = "Log2FoldChange (log2FC) cutoff",
                               min = 0,
                               max = 5,
                               step = 0.1 ,
                               value = 2),
                   width = 6 ),
                 
                 box (
                   actionButton(
                     "Run_whole_data_inspection",
                     icon = icon("seedling"),
                     strong("Run Whole Data Inspection")),
                   width = 12)
                 
                 
               ),
               
               box (
                 title = "Input Data",
                 closable = FALSE, 
                 collapsible = TRUE,
                 collapsed = TRUE, 
                 status = "success",
                 dataTableOutput("Table_whole_data"), 
                 width = 12), 
               
               box(
                 title = "Volcano Plot",
                 status = "success",
                 shinycustomloader::withLoader(plotlyOutput("volcanoPlot_plotly", height = "450px"), type = "image", loader = "wait.gif"),
                 width = 6
               ), #fin box
               
               
               box(
                 title = "MA Plot",
                 status = "success",
                 shinycustomloader::withLoader(plotlyOutput("MAPlot_plotly", height = "450px"), type = "image", loader = "wait.gif"),
                 width = 6
               ), #fin box

               

               box (
                 title = "DEG (using padj and log2FC cutoff selected) used for subsequent analysis",
                 closable = TRUE, 
                 collapsible = TRUE,
                 collapsed = TRUE,
                 status = "success",
                 dataTableOutput("Table_DEG"), 
                 downloadButton('download_Table_DEG', icon = icon("leaf"),"Download"),
                 width = 12),
               
               box (
                 title = "Selected data (points selected on the plot (Volcanoplot or MaPlot)",
                 closable = TRUE, 
                 collapsible = TRUE,
                 collapsed = TRUE,
                 status = "success",
                 dataTableOutput("Table_data_selected"), 
                 downloadButton('download_Table_data_selected', icon = icon("leaf"), "Download"),
                 width = 12),
            
               dataTableOutput("id_toto")
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
                  width = 7
                  ),
                   
                   
                   box(
                   selectInput( inputId = "kegg_adj_method",
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
                                selected = 'none'),
                   width = 3
                   ),
                   
                   box(
                     actionButton("Run_Annotation_ENSEMBL_to_GO","Run"),
                     width = 2
                   )
                   
                   
               ),
               
               box (title = "Gene annotation with KEGG",
                 dataTableOutput("enrichKEGG_table"), 
                    width = 12, style = "overflow-x: scroll;"),
               
               box(tstatus = "warning", solidHeader = TRUE, width = 12, height = "550px",
                   shinycustomloader::withLoader(plotlyOutput("dotplot_kegg", height = "450px"), type = "image", loader = "wait.gif")
               ),
               
               
               mainPanel(title = "GSE Plot", status = "warning", solidHeader = TRUE, width = 12, height = "550px",
                         shinycustomloader::withLoader(plotlyOutput("method_kegg", height = "450px"), type = "image", loader = "wait.gif")
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
                 title = strong("Parameters - Protein Domain Enrichment"), 
                 closable = TRUE, 
                 status = "success",
                 width = 12,
                 solidHeader = TRUE, 
                 collapsible = TRUE,
                 

                 br(),
                 
                 box (
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
                   
                   conditionalPanel(
                     
                     
                     
                     condition = "input.method_prt_domain == 1", 
                     
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
                                  selected = NULL),
                     ), 
                   width = 3
                 ),
                 
                 
                 box (
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
                   ),
                   width = 3
                 ),
                 
                 
                 box (
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
                                 ), "Thresholding for subsequent analysis (ORA or GSEA)", placement="bottom", trigger = "hover") 
                               ),
                               min = 0,
                               max = 1,
                               value = 0.05),
                   width = 6
                 ),
                 
                 
                 
                 box (
                   status = "success",
                   actionButton("Run_protein_domains",h4("Run Protein Domains")),
                   width = 12)
               ),
               
               
               

               
               conditionalPanel(
                 
                 
                 
                 condition = "input.method_prt_domain == 1",
                 
                 box(
                   title = h3(strong("Over-representation (or enrichment) analysis"),  br(),
                              tags$a(href="https://github.com/AlizeeBardon/shiny_enrichment_analysis","coded method" , style="color:white"),
                              align = "center"),
                   background = "black", solidHeader = TRUE,
                   width = 12
                 ),

                 
                 
                 box (
                   title = strong("Result table"),
                   dataTableOutput("Table_domains_enrichment"),
                   br(), br(), 
                   
                   sliderInput(inputId = "nb_barplot_ora_coder",
                               label = strong("Number of protein domains to see in the barplot and dotplot" ),
                               min = 0,
                               max = 200,
                               value = 10),
                   width = 12),
                 
                 box (
                   title = strong("BarPlot"), 
                   shinycustomloader::withLoader(plotlyOutput("barplot_domains_enrichment", height = "450px"), type = "image", loader = "wait.gif"),
                   width = 6),
                
                 box (
                   title = strong("DotPlot"),
                   shinycustomloader::withLoader(plotlyOutput("dotplot_domains_enrichment", height = "450px"), type = "image", loader = "wait.gif"),
                   width = 6),
                 
                 box (
                   title = strong("PieChart"),
                   shinycustomloader::withLoader(plotlyOutput("piechart_domains_enrichment", height = "450px"), type = "image", loader = "wait.gif"),
                   width = 6),

                 box (
                   title = strong("WordCloud"),
                   plotOutput("cloud"),
                   width = 6),
                 
                 
                 br(), br(),
                 

                 box(
                   title = h2( strong("Over-representation (or enrichment) analysis"), 
                               br(),
                               tags$a(href="https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/enricher","enricher() function from clusterProfiler", style="color:white"),
                               align = "center"),
                     
                   background = "black", solidHeader = TRUE,
                   width = 12
                 ),
                 
                 box (
                   title = strong("Result Table"), 
                   dataTableOutput("Table_domains_enrichment_enricher"), 
                   width = 12),
                 
                 box (
                   title = strong("BarPlot"), 
                   shinycustomloader::withLoader(plotlyOutput("barplot_domains_enrichment_enricher", height = "450px"), type = "image", loader = "wait.gif"),
                   width = 6),
                 
                 box (
                   title = strong("DotPlot"), 
                   shinycustomloader::withLoader(plotlyOutput("dotplot_domains_enrichment_enricher", height = "450px"), type = "image", loader = "wait.gif"),
                   width = 6)
                 
                 
               ),
               
               
               conditionalPanel(
                 condition = "input.method_prt_domain == 2",
                 
                 box(
                   title = h2(strong("Gene Set Enrichment Analysis (GSEA) - using the method enricher() from clusterProfiler method"),  br(),
                              tags$a(href="https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/GSEA","GSEA() function from clusterProfiler","coded method", style="color:white"),
                              align = "center"),
                   background = "black", solidHeader = TRUE,
                   width = 12
                 ),
                 
                 br(), br(),
                 
                 box (
                   title = strong("Result table"),
                   dataTableOutput("Table_domain_enrichment_GSEA"), 
                   width = 12),
                 
                 box (
                   title = strong("BarPlot"), 
                   shinycustomloader::withLoader(plotlyOutput("barplot_domain_enrichment_GSEA", height = "450px"), type = "image", loader = "wait.gif"),
                   width = 6),
                 
                 box (
                   title = strong("BarPlot"), 
                   selectInput(
                     "protein_id_list", 
                     label = h4("choose the protein id to display"),
                     choices = "",
                     selected = NULL),
                   shinycustomloader::withLoader(plotlyOutput("gseaplot_domain_enrichment_GSEA", height = "450px"), type = "image", loader = "wait.gif"),
                   width = 12),
                     
                 

                 
                 br(), br(),
                 
               )  
               

               
                   
      ) # tabPanel("Protein Domain Enrichment"
      
      
    ) # fin tabsetPanel
    
  ) # fin dashboardBody(
  
  
) # fin dashboardPage(
) # fin shinyUI(
