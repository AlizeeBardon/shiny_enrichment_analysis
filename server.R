#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("pasilla") 
# BiocManager::install("biomaRt")


organism = "org.Mm.eg.db" 
#BiocManager::install(organism, character.only = TRUE) 
library(organism, character.only = TRUE) 
library(org.Mm.eg.db)
library(DOSE)
library(clusterProfiler) 
library(biomaRt)
library(pathview)
library(shiny)
library(shinythemes)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(highcharter)
library(DT)
library(ggplot2)
library(enrichplot)


# Define server 
shinyServer(function(input, output, session) {

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


    # BODY --------------------------------------------------------------------
    
        # BODY: tabPanel : whole Data Analysis --------------------------------


# Volcanoplot 


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
    

# MA plot    

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
    

# Interactive Table

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
    
# Annotation Table
    
    annot <- eventReactive(input$Run_Annotation, {
        data <- re() 
        gene_list = data$GeneName                 
        organism = input$espece
        generef = bitr(gene_list, fromType = "SYMBOL", toType= "GO", OrgDb=organism)
    })

    output$annotation <- renderDataTable({
      D <- annot()
      DT::datatable(D)
    })#fin renderDataTable

    # -------------------------------------------------------------------
    # BODY: tabPanel :GO Term Enrichment --------------------------------
    # -------------------------------------------------------------------
    ####################################################
    ### GSEA pour GO
    #####################################################
    
    goGse_annot<-reactive({
      # reading in data
      go <- re()
      
      # we want the log2 fold change
      original_gene_list <- go$log2FC
      
      # name the vector
      names(original_gene_list) <- go$ID
      
      # omit any NA values 
      gene_list<-na.omit(original_gene_list)
      
      # sort the list in decreasing order (required for clusterProfiler)
      gene_list <- sort(gene_list, decreasing=TRUE)
      
      
      gse <- gseGO(geneList=gene_list, 
                   ont ="BP", 
                   keyType = 'ENSEMBL', 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = "org.Mm.eg.db", 
                   pAdjustMethod = "none")
    })
    
    output$goGse_annot_table <- renderDataTable({
      data <- goGse_annot()
      updateSelectInput(session, "paths", choices = data$Description)
      data.df <- as.data.frame(data)
      DT::datatable(data.df)
    })#fin renderDataTable 
  
    output$dotplot <-renderPlotly({
      gse<-goGse_annot()
      require(DOSE)
      dotplot(gse, showCategory = 7, title = "gsea dotplot" , split=".sign") + facet_grid(.~.sign)
    })

    output$ridgeplot <-renderPlotly({
      gse<-goGse_annot()
      ridgeplot(gse, showCategory = 7)
    })
    
    output$gsea_plot <-renderPlotly({
      gse<-goGse_annot()
      gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
    })

    ####################################################
    ### ORA pour GO
    #####################################################
    goGse_enrich<-reactive({
      
      # reading in data
      go <- re()
      
      # we want the log2 fold change 
      original_gene_list <- df$log2FC
      
      # name the vector
      names(original_gene_list) <- df$ID
      
      # omit any NA values 
      gene_list<-na.omit(original_gene_list)
      
      # sort the list in decreasing order (required for clusterProfiler)
      gene_list = sort(gene_list, decreasing = TRUE)
      
      head(gene_list)
      # Exctract significant results (padj < 0.05)
      sig_genes_df = subset(df, padj < 0.05)
      
      # From significant results, we want to filter on log2fold change
      genes <- sig_genes_df$log2FC
      
      # Name the vector
      names(genes) <- sig_genes_df$ID
      
      # omit NA values
      genes <- na.omit(genes)
      
      genes <- names(genes)
      
      go_enrich <- enrichGO(gene = genes,
                            universe = names(gene_list),
                            OrgDb = "org.Mm.eg.db", 
                            keyType = 'ENSEMBL',
                            readable = T,
                            ont = "BP",
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.10)
    })
    
    output$goGse_enrich_table <- renderDataTable({
      data <- goGse_enrich()
      updateSelectInput(session,"paths", choices = data$Description)
      data.df <- as.data.frame(data)
      DT::datatable(data.df)
    })#fin renderDataTable 
    
    output$barplot <-renderPlotly({
      gse<-goGse_enrich()
      barplot(gse, showCategory = 7)+ ggtitle("barplot for SEA")
    })
    
    output$dotplot_sea <-renderPlotly({
      gse<-goGse_enrich()
      dotplot(gse, showCategory = 7)+ ggtitle("barplot for SEA")
    })
    
    output$usetplot <-renderPlotly({
      gse<-goGse_enrich()
      upsetplot(gse, n=2)+ ggtitle("usetplot for SEA")
    })
    
    output$goplot <-renderPlotly({
      gse<-goGse_enrich()
      goplot(go_enrich, 
             #drop = TRUE, 
             showCategory = 6, #nombre de pathway à afficher
             #title = "GO Biological Pathways",
             font.size = 8,
             split=".sign")
    })

    
    } # end function(input, output) {
) # end shinyServer(
