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
library(shinyalert)



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
      #list_df <- c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj")
      data <- read.csv(file$datapath, header = TRUE, sep = ";") %>% 
        mutate(
          minusLog10Pvalue = -log10(padj)
          )

      required_columns <- c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj")
      column_names <- colnames(data)
      min_columns <- 6
      
      shiny::validate(
        need(ncol(data) >= min_columns, "Your data has not enought columns. Your data must contain : GeneName, ID, baseMean, log2FC, pval, padj"),
        need(all(required_columns %in% column_names), "You don't have the right data.  Your data must contain : GeneName, ID, baseMean, log2FC, pval, padj")
      )
      
      data
      
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

biomart_listMarts <- reactive({
  biomart_listMarts <- input$biomart_listMarts
  req(biomart_listMarts)
})

biomart_dataset <- reactive({
  biomart_dataset <- input$biomart_dataset
  req(biomart_dataset)
})


    # BODY --------------------------------------------------------------------
    
        # BODY: tabPanel : whole Data Analysis --------------------------------

### Volcanoplot 

    output$volcanoPlot_plotly <- renderPlotly({
        #input data
        data_plot_plotly <- re()
        data_plot_plotly$col <- "black"
        pvalue <- pvalue()
        pvalue_log10 <- -log10(pvalue)
        tresholdlog2foldchange <- tresholdLog2FoldChange()
        data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue_log10 & data_plot_plotly$log2FC > tresholdlog2foldchange, "col"] <- "green"
        data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue_log10 & data_plot_plotly$log2FC < -tresholdlog2foldchange, "col"] <- "red"
        data_plot_plotly$key <- row.names(data_plot_plotly)
        
        # volcanoPlot using ggplot and plotly
        p <- ggplot(
          data = data_plot_plotly, 
          mapping = aes(x = log2FC, y = minusLog10Pvalue, gene= GeneName, ID= ID, col = I(col), key = key)) +
          geom_point(size = 0.45) + 
          theme_bw()  + 
          geom_hline(yintercept=pvalue_log10, linetype="dashed", color = "black", size=0.1) +
          geom_vline(xintercept=tresholdlog2foldchange, linetype="dashed", color = "black", size=0.1) +
          geom_vline(xintercept=-tresholdlog2foldchange, linetype="dashed", color = "black", size=0.1) +
          xlab("log2 fold change") + 
          ylab("-log10(p-value)") 
        
        # preparation of the interactive table (save the selected point into variables)
        event_register(p, 'plotly_click')
        event_register(p, 'plotly_selected')
        plotly_object <- ggplotly(p,source = "source1") %>% 
          layout(dragmode = "lasso")  %>% 
          layout(showlegend = FALSE)
        
        # updating the name of the saved plot 
        config(plotly_object,
               toImageButtonOptions= list(filename = paste0("VolcanoPlot_pvalue_", pvalue, "_log2FC_", tresholdlog2foldchange)))

      })

### MA plot    

    output$MAPlot_plotly <- renderPlotly({
      #input data
      data_plot_plotly <- re()
      data_plot_plotly$key <- row.names(data_plot_plotly)
      data_plot_plotly$col <- "black"
      pvalue <- pvalue()
      pvalue_log10 <- -log10(pvalue)
      tresholdlog2foldchange <- tresholdLog2FoldChange()
      data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC > tresholdlog2foldchange, "col"] <- "green"
      data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC < -tresholdlog2foldchange, "col"] <- "red"
      
      # MAplot using ggplot and plotly
      p <- ggplot(
        data = data_plot_plotly, 
        mapping = aes(
          x = log2(baseMean), y = log2FC, gene= GeneName, ID= ID, col = I(col), key = key)) +
        geom_point(size = 0.45)+ 
        theme_bw()  +
        xlab("log2(baseMean)") + 
        ylab("log2FC")
      
      # preparation of the interactive table (save the selected point into variables)
      event_register(p, 'plotly_click')
      event_register(p, 'plotly_selected')      
      plotly_object <- ggplotly(p,source = "source1") %>% 
        layout(dragmode = "lasso")  %>% 
        layout(showlegend = FALSE)

      # updating the name of the saved plot 
      config(plotly_object,
             toImageButtonOptions= list(filename = paste0("MAPlot_pvalue_", pvalue, "_log2FC_", tresholdlog2foldchange)))
      
    })
    

### Interactive Table

    output$Table_subset_data_selected <- renderDataTable({ 
      D <- re() %>%  
        mutate(
          indice = row_number()
        )
      data_selected_from_graph <- event_data("plotly_selected",source = "source1")
      D_subset = subset(D, indice %in%  data_selected_from_graph$key)
      DT::datatable(D_subset) 
    }) # fin renderDataTable({
    
    #surligner la ligne correspondant au point clique sur le graphique
    proxy <- DT::dataTableProxy("Table_subset_data_selected")
    observe({
      data <- re()
      subset_data <- data.frame(
        ID_test = data$ID,
        row_id = 1:length(data$ID),
        stringsAsFactors = FALSE )
      s <- event_data("plotly_click",source = "source1")
      req(!is.null(s))
      # map point number to subset_data
      row_clicked <- subset_data[s$pointNumber + 1,"row_id"]
      proxy %>%
        selectRows(NULL) %>%
        selectRows(row_clicked)
    }) # fin observe
    




    # -------------------------------------------------------------------
    # BODY: tabPanel :GO Term Enrichment --------------------------------
    # -------------------------------------------------------------------
         
    
    
    
    # -------------------------------------------------------------------
    # BODY: tabPanel :KEGG --------------------------------
    # -------------------------------------------------------------------

    # Annotation Table
    
    kegg_data <- eventReactive(input$Run_Annotation_ENSEMBL_to_GO, {
      data <- re() 
      orga = input$espece
      generef = bitr(data$ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=orga)
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
      
      df2 = df[df$ID %in% dedup_ids$ENSEMBL,]
      
      # Create a new column in df2 with the corresponding ENTREZ IDs
      df2$Y = dedup_ids$ENTREZID
      
      # Create a vector of the gene unuiverse
      kegg_gene_list <- df2$log2FC

      # Name vector with ENTREZ ids
      names(kegg_gene_list) <- df2$Y
      # omit any NA values 
      kegg_gene_list<-na.omit(kegg_gene_list)
      # sort the list in decreasing order (required for clusterProfiler)
      kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
      # if (input$type == 1){
      #   #kegg_gene_list <- names(kegg_gene_list)[abs(kegg_gene_list) > 0]
      #   deg_list <- kegg_gene_list[kegg_gene_list > 0]
      # }
      # else if(input$type == 2) {
      #   deg_list <- kegg_gene_list[kegg_gene_list < 0]
      # }
      # else {
      #   deg_list <- kegg_gene_list
      # }
      list(kegg_gene_list=kegg_gene_list, df2=df2)
    })
    
   
    gse_kegg <- reactive({
      gse_kegg_data <- kegg_data()
      gse_kegg_annal <- gseKEGG(geneList     = gse_kegg_data$kegg_gene_list,
                          organism     = "mmu",
                          nPerm        = 10000,
                          minGSSize    = 3,
                          maxGSSize    = 800,
                          pvalueCutoff = input$pvalue_gsea,
                          pAdjustMethod = "none",
                          keyType       = "ncbi-geneid")

    })
     
    
    ora_kegg <- reactive({
      ora_kegg_data <- kegg_data()
      # Exctract significant results from df2
      kegg_sig_genes_df = subset(ora_kegg_data$df2, padj < 0.05)

      # From significant results, we want to filter on log2fold change

      kegg_genes <- kegg_sig_genes_df$log2FC

      # Name the vector with the CONVERTED ID!
      names(kegg_genes) <- kegg_sig_genes_df$Y

      # omit NA values
      kegg_genes <- na.omit(kegg_genes)
      browser()
      # filter on log2fold change (under or over expressed DEG)
      if (input$type == 1){
        kegg_genes <- names(kegg_genes)[kegg_genes > 0]
      }
      else if(input$type == 2) {
        kegg_genes <- names(kegg_genes)[kegg_genes < 0]
      }
       kk <- enrichKEGG(gene=kegg_genes, universe=names(ora_kegg_data$kegg_gene_list),organism='mmu', pvalueCutoff = input$pvalue_gsea, keyType = "ncbi-geneid")
    })
    
    output$enrichKEGG_table <- renderDataTable({
      if (input$method == 1){
        data <- ora_kegg()
      }
      else {
        data <- gse_kegg()
      }
      updateSelectInput(session, "paths", choices = data$Description)
      data.df <- as.data.frame(data)
      DT::datatable(data.df)
    })#fin renderDataTable
    
     output$dotplot_kegg <- renderPlotly({
       if(input$method == 1){
         result_kegg <- ora_kegg()
       }
       else {
         result_kegg <- gse_kegg()
       }
       dotplot(result_kegg, showCategory = 10, title = "Enriched Pathways")
     })


     output$method_kegg <- renderPlotly({
       if (input$method == 1){
         result_kegg <- ora_kegg()
         barplot(result_kegg, showCategory = 10)
       }
       else {
         result_kegg <- gse_kegg()
         gseaplot(result_kegg, by = "all", title =input$paths, geneSetID = result_kegg[result_kegg$Description==input$paths,]$ID)
       }
     })


     output$pathview_kegg <- renderImage({
       if (input$method == 1){
         result_kegg <- ora_kegg()
       }
       else {
         result_kegg <- gse_kegg()
       }
       path_id<-result_kegg[result_kegg$Description==input$paths,]$ID
       # Make and save picture of pathway with pathview  
       pathview(cpd.data=kegg_gene_list, pathway.id=path_id, species = "mmu")
       path_img<-paste("./",path_id,".png", sep="")
       # Return picture path to load it on popup window
       list(src = path_img)
     }, deleteFile = TRUE)


     # -------------------------------------------------------------------
     # BODY: tabPanel : Protein DOmains   --------------------------------
     # -------------------------------------------------------------------
     

     domain_enrichment <-  eventReactive(input$Run_protein_domains, {
       #input data
       resOrdered <- re()
       pvalue <- pvalue()
       espece <- espece()
       biomart_dataset <- biomart_dataset()
       biomart_listMarts <- biomart_listMarts()
       
       # recuperation des domain ID pour les ensembl ID de notre jeu de donnees 
       ensembl = useMart(biomart_listMarts,dataset=biomart_dataset)
       interpro_id <- getBM(
         attributes=c('interpro', 'ensembl_gene_id'), # namespace_1003 = for go domain
         filters = 'ensembl_gene_id', #ensembl_gene_id
         values = resOrdered$ID, 
         mart = ensembl)
       #background
       GeneRef = interpro_id
       
       #liste d'interet
       GeneList = resOrdered[which(resOrdered$padj<=pvalue),]$ID
       GeneList = data.frame(GeneList)

       ###  I - prepare data for enrichment      

       get_Gene_and_Bg_ratio = function(GeneList, GeneRef) {
         # reference list
         # m : nb of annotated genes in the reference list (for each term)
         m = table(GeneRef$interpro)
         # n : nb of non annotated genes in the reference list
         n = length(unique(GeneRef$ensembl_gene_id)) - m
         # experience (interest list)
         # x : nb of annotated genes in the interest list
         experience = merge(GeneList, GeneRef, by.x = "GeneList", by.y = "ensembl_gene_id")
         x = table(factor(experience$interpro, rownames(m)))
         # k : total nb of genes in the interest list
         k = length(unique(GeneList$GeneList))
         
         Term = unique(GeneRef$interpro)
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
       
       ###  II -  hypergeometric test                        

       hypergeom_test = function(x, k, m, n){
         # calculate p-value and adjusted p-value
         pvalue = phyper(x-1,m,n,k,lower.tail=FALSE)
         padj = p.adjust(pvalue, n=length((pvalue)))
         return (list(pvalue = pvalue, padj = padj))
       }
       
       ### III - create results table                        
       
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
  # /!\ ajouter : interpro description + pourcentage 
       } )  
     
     
     output$Table_domains_enrichment <- renderDataTable({ 
       D <- domain_enrichment()
       DT::datatable(D) 
     }) # fin renderDataTable({
     
     
     output$barplot_domains_enrichment <- renderPlotly({
       D <- domain_enrichment()
       ggplot(data = D, aes(x= count , y = reorder(Term, count)  ) ) +
         geom_bar(stat = "identity", aes(fill = padj))  +
         theme(axis.text.x = element_text(
           angle = 90,
           hjust = 1,
           vjust = 0.5
         ))
       
     })
     
    
    } # end function(input, output) {
) # end shinyServer(
