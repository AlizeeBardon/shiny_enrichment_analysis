#' pathways UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_pathways_ui <- function(id){
  ns <- NS(id)
  
  tagList(
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
 
  )
}
    
#' pathways Server Functions
#'
#' @noRd 
mod_pathways_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_pathways_ui("pathways_ui_1")
    
## To be copied in the server
# mod_pathways_server("pathways_ui_1")
