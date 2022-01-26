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
#BiocManager::install("biomaRt")
  

organism = "org.Hs.eg.db" 
#BiocManager::install(organism, character.only = TRUE) 
library(organism, character.only = TRUE) 
library(org.Mm.eg.db)
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



# Define server 
shinyServer(function(input, output) {

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

pvalue <- reactive({
  pvalue <- input$pvalue
  req(pvalue)
})

tresholdLog2FoldChange <- reactive({
  tresholdLog2FoldChange <- input$tresholdLog2FoldChange
  req(tresholdLog2FoldChange)
})

radio_filtre_ontology <- reactive({
  radio_filtre_ontology <- input$radio_filtre_ontology
  req(radio_filtre_ontology)
})

subset_up_down_regulated <- reactive({
  subset_up_down_regulated <- input$subset_up_down_regulated
  req(subset_up_down_regulated)
})


espece <- reactive({
  espece <- input$espece
  req(espece)
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
        
        pvalue <- pvalue()
        pvalue_log10 <- -log10(pvalue)
        tresholdlog2foldchange <- input$tresholdLog2FoldChange
          
        data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue_log10 & data_plot_plotly$log2FC > tresholdlog2foldchange, "col"] <- "green"
        data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue_log10 & data_plot_plotly$log2FC < -tresholdlog2foldchange, "col"] <- "red"
        
        p <- ggplot(data = data_plot_plotly, mapping = aes(x = log2FC, y = minusLog10Pvalue, col = I(col), key = key)) +
          geom_point(size = 0.45) + 
          theme_bw()  + 
          geom_hline(yintercept=pvalue_log10, linetype="dashed", color = "black", size=0.1) +
          geom_vline(xintercept=tresholdlog2foldchange, linetype="dashed", color = "black", size=0.1) +
          geom_vline(xintercept=-tresholdlog2foldchange, linetype="dashed", color = "black", size=0.1) +
          xlab("log2 fold change") + 
          ylab("-log10(p-value)") 
        
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
    
    output$pvalue_go_enrich <- renderPrint({ 
      pvalue() })
    
    output$log2foldchange_go_enrich <- renderPrint({ input$tresholdLog2FoldChange })
    
    output$subset_annotation <- renderDataTable({ 
      resOrdered <- re()
      pvalue <- pvalue()
      tresholdLog2FoldChange <- tresholdLog2FoldChange()
      subset_up_down_regulated <- subset_up_down_regulated()
      
      GeneList = resOrdered[which(resOrdered$padj<=pvalue),]$ID
      # genes considered DE with treshold log2foldchange (up or down regulated define by user (default= 0.4)
      if (input$subset_up_down_regulated == "both"){
        GeneList = resOrdered[which(resOrdered$log2FC< -tresholdLog2FoldChange | resOrdered$log2FC > tresholdLog2FoldChange),]$ID
      } else if(input$subset_up_down_regulated == "up"){
        GeneList = resOrdered[which(resOrdered$log2FC > tresholdLog2FoldChange),]$ID
      } else if(input$subset_up_down_regulated == "down"){
        GeneList = resOrdered[which(resOrdered$log2FC < -tresholdLog2FoldChange),]$ID
      }
      genes = resOrdered$ID
      GeneRef =  bitr(genes, fromType="ENSEMBL", toType="GO", OrgDb="org.Mm.eg.db")
      filtre_ontology = input$radio_filtre_ontology
      if ( filtre_ontology != "all") {
        GeneRef <- subset(GeneRef, ONTOLOGY == filtre_ontology)
      }
      
      DT::datatable(GeneRef)
      
      })
    
    
    
    output$Table_go_enrichment <- renderDataTable({
      resOrdered <- re()
      pvalue <- pvalue()
      tresholdLog2FoldChange <- tresholdLog2FoldChange()
      subset_up_down_regulated <- subset_up_down_regulated()
      
      # get the interest list
      # genes considered DE with treshold alpha define by user (default= 0.05)
      GeneList = resOrdered[which(resOrdered$padj<=pvalue),]$ID
      # genes considered DE with treshold log2foldchange (up or down regulated define by user (default= 0.4)
      if (input$subset_up_down_regulated == "both"){
        GeneList = resOrdered[which(resOrdered$log2FC< -tresholdLog2FoldChange | resOrdered$log2FC > tresholdLog2FoldChange),]$ID
      } else if(input$subset_up_down_regulated == "up"){
        GeneList = resOrdered[which(resOrdered$log2FC > tresholdLog2FoldChange),]$ID
      } else if(input$subset_up_down_regulated == "down"){
        GeneList = resOrdered[which(resOrdered$log2FC < -tresholdLog2FoldChange),]$ID
      }
        
        
      GeneList = data.frame(Gene = GeneList)

      # get gene annotation (for all genes)
      genes = resOrdered$ID
      GeneRef =  bitr(genes, fromType="ENSEMBL", toType="GO", OrgDb="org.Mm.eg.db")
      filtre_ontology <- radio_filtre_ontology()
      if ( filtre_ontology != "all") {
        GeneRef <- subset(GeneRef, ONTOLOGY == filtre_ontology)
      }
      
      
      #################################################
      #  prepare data for enrichment                  #
      #################################################
      get_Gene_and_Bg_ratio = function(GeneList, GeneRef) {
        # reference list
        # m : nb of annotated genes in the reference list (for each term)
        m = table(GeneRef$GO)
        # n : nb of non annotated genes in the reference list
        n = length(unique(GeneRef$ENSEMBL)) - m
        
        # experience (interest list)
        # x : nb of annotated genes in the interest list
        experience = merge(GeneList, GeneRef, by.x = "Gene", by.y = "ENSEMBL")
        x = table(factor(experience$GO, rownames(m)))
        # k : total nb of genes in the interest list
        k = length(unique(GeneList$Gene))
        
        Term = unique(GeneRef$GO)
        x = as.numeric(x)
        m = as.numeric(m)
        k = as.numeric(k)
        n = as.numeric(n)
        
        return(list(Term = Term, 
                    x = x, 
                    k = k, 
                    m = m, 
                    n = n))
      }
      Gene.Bg.ratio = get_Gene_and_Bg_ratio(GeneList = GeneList, GeneRef = GeneRef)
      Bg.ratio = signif(100 * Gene.Bg.ratio$m/(Gene.Bg.ratio$m + Gene.Bg.ratio$n), 3)
      Gene.ratio = signif(100 * Gene.Bg.ratio$x / Gene.Bg.ratio$k, 3)
      
      
      
      #################################################
      #  hypergeometric test                          #
      #################################################
      hypergeom_test = function(x, k, m, n){
        # calculate p-value and adjusted p-value
        pvalue = phyper(x-1,m,n,k,lower.tail=FALSE)
        padj = p.adjust(pvalue, n=length((pvalue)))
        
        return (list(pvalue = pvalue, padj = padj))
      }
      res_hypergeom_test = hypergeom_test(x = Gene.Bg.ratio$x, 
                                          k = Gene.Bg.ratio$k,
                                          m = Gene.Bg.ratio$m, 
                                          n = Gene.Bg.ratio$n)
      
      #################################################
      #  create results table                         #
      #################################################
      create_table_enrichment = function(GeneList, GeneRef){
        # call function get_Gene_and_Bg_ratio() to get BgRatio and GeneRatio 
        Gene.Bg.ratio = get_Gene_and_Bg_ratio(GeneList = GeneList, GeneRef = GeneRef)
        Bg.ratio = signif(100 * Gene.Bg.ratio$m/(Gene.Bg.ratio$m + Gene.Bg.ratio$n), 3)
        Gene.ratio = signif(100 * Gene.Bg.ratio$x / Gene.Bg.ratio$k, 3)
        
        # call function hypergeom_test to get p-value and adjusted p-value
        test = hypergeom_test(x = Gene.Bg.ratio$x, 
                              k = Gene.Bg.ratio$k,
                              m = Gene.Bg.ratio$m, 
                              n = Gene.Bg.ratio$n)
        
        # create dataframe
        table.enrich = data.frame(Term = Gene.Bg.ratio$Term, 
                                  GeneRatio = Gene.ratio, 
                                  BgRatio = Bg.ratio,
                                  pval = test$pvalue, 
                                  padj = test$padj, 
                                  count = Gene.Bg.ratio$x)
        
        return (table.enrich[order(table.enrich$pval), ])
      }
      res.enrich.hypergeom.GO = create_table_enrichment(GeneList = GeneList, GeneRef = GeneRef)
      res.enrich.hypergeom.GO[which(res.enrich.hypergeom.GO$padj<0.05),]
      
      DT::datatable(res.enrich.hypergeom.GO[which(res.enrich.hypergeom.GO$padj<0.05),])
      
    })#fin renderDataTable
    
    
    
    
    # -------------------------------------------------------------------
    # BODY: tabPanel :KEGG --------------------------------
    # -------------------------------------------------------------------
    
    
    
    
    
    } # end function(input, output) {
) # end shinyServer(
