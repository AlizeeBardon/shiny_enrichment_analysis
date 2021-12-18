#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
  ######################################################################### 
  ###### input data    
  #########################################################################
  
  re <- reactive({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    req(file)
    print(head(file))
    validate(need(ext == "csv", "Invalid file. Please upload a .csv file"))
    data <- read.csv(file$datapath, header = TRUE, sep = ";") %>% 
      mutate(
        minusLog10Pvalue = -log10(padj)
      )
  })
  
  ######################################################################### 
  ###### volcanoplot using ploty    
  #########################################################################
  
  output$volcanoPlot_plotly <- renderPlotly({
    
    data_plot_plotly <- re()
    data_plot_plotly$key <- row.names(data_plot_plotly)
    data_plot_plotly$col <- "black"
    
    click_data <- event_data("plotly_click")
    select_data <- event_data("plotly_selected")
    
    pvalue <- -log10(input$pvalue)
    
    data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC > 0, "col"] <- "green"
    data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC < 0, "col"] <- "red"
    
    p <- ggplot(data = data_plot_plotly, mapping = aes(x = log2FC, y = minusLog10Pvalue, col = I(col), key = key)) +
      geom_point(size = 0.45) + 
      theme_bw()  +
      xlab("log2 fold change") + ylab("-log10(p-value)")
    
    event_register(p, 'plotly_click')
    event_register(p, 'plotly_selected')
    
    plotly_object <- ggplotly(p,source = "source1") %>% 
      layout(dragmode = "lasso")  %>% 
      layout(showlegend = FALSE)
    
  })
  
  #########################################################################
  ###### MA plot using ploty TEST   
  #########################################################################
  
  output$MAPlot_plotly <- renderPlotly({
    
    data_plot_plotly <- re()
    data_plot_plotly$key <- row.names(data_plot_plotly)
    data_plot_plotly$col <- "black"
    
    click_data <- event_data("plotly_click")
    select_data <- event_data("plotly_selected")
    
    pvalue <- -log10(input$pvalue)
    
    data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC > 0, "col"] <- "green"
    data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC < 0, "col"] <- "red"
    
    p <- ggplot(data = data_plot_plotly, mapping = aes(
      x = log2(baseMean), y = log2FC, col = I(col), key = key)) +
      geom_point(size = 0.45) + 
      theme_bw()  +
      xlab("log2(baseMean)") + 
      ylab("log2FC")
    
    plotly_object <- ggplotly(p,source = "source1") %>% 
      layout(dragmode = "lasso")  %>% 
      layout(showlegend = FALSE)
    
  })
  
  #########################################################################
  ### TABLEAU INTERACTIF
  ######################################################################### 
  
  output$Table_subset_data_selected <- renderDataTable({ 
    D <- re() %>%  
      mutate(
        indice = row_number()
      )
    data_selected_from_graph <- event_data("plotly_selected",source = "source1")
    
    D_subset = subset(D, indice %in%  data_selected_from_graph$key)
    DT::datatable(D_subset) 
  }) # fin renderDataTable({
  
  
  proxy <- DT::dataTableProxy("Table_subset_data_selected")
  
  observe({
    data <- re()
    subset_data <- data.frame(
      ID_test = data$ID,
      row_id = 1:length(data$ID),
      stringsAsFactors = FALSE
    )
    s <- event_data("plotly_click",source = "source1")
    req(!is.null(s))
    # map point number to subset_data
    row_clicked <- subset_data[s$pointNumber + 1,"row_id"]
    proxy %>%
      selectRows(NULL) %>%
      selectRows(row_clicked)
    
  }) # fin observe
  #mod_pathways_server("pathways_ui_1")
  # Your application server logic 
  
}
