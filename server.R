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
      #list_df <- c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj")
      data <- read.csv(file$datapath, header = TRUE, sep = ";") %>% 
        mutate(
          minusLog10Pvalue = -log10(padj)
          )
      validate(need(ncol(data) == 7, "Invalid file. Please upload a file with 6 col"))
      
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


    # BODY --------------------------------------------------------------------
    
        # BODY: tabPanel : whole Data Analysis --------------------------------


# Volcanoplot 


    output$volcanoPlot_plotly <- renderPlotly({
        #input data
        data_plot_plotly <- re()
        data_plot_plotly$col <- "black"
        pvalue <- pvalue()
        pvalue_log10 <- -log10(pvalue)
        tresholdlog2foldchange <- tresholdLog2FoldChange()
        data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue_log10 & data_plot_plotly$log2FC > tresholdlog2foldchange, "col"] <- "green"
        data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue_log10 & data_plot_plotly$log2FC < -tresholdlog2foldchange, "col"] <- "red"
        
        # preparation of the interactive table (save the selected point into variables)
        data_plot_plotly$key <- row.names(data_plot_plotly)
        click_data <- event_data("plotly_click")
        select_data <- event_data("plotly_selected")
          

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
        
        
        config(plotly_object,
               toImageButtonOptions= list(filename = paste0("VolcanoPlot_pvalue_", pvalue, "_log2FC_", tresholdlog2foldchange)))

      })
    

# MA plot    

    output$MAPlot_plotly <- renderPlotly({
      
      data_plot_plotly <- re()
      data_plot_plotly$key <- row.names(data_plot_plotly)
      data_plot_plotly$col <- "black"
      
      click_data <- event_data("plotly_click")
      select_data <- event_data("plotly_selected")
      
      pvalue <- pvalue()
      pvalue_log10 <- -log10(pvalue)
      tresholdlog2foldchange <- tresholdLog2FoldChange()
      
      data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC > tresholdlog2foldchange, "col"] <- "green"
      data_plot_plotly[data_plot_plotly$minusLog10Pvalue > pvalue & data_plot_plotly$log2FC < -tresholdlog2foldchange, "col"] <- "red"
      
      p <- ggplot(data = data_plot_plotly, mapping = aes(
        x = log2(baseMean), y = log2FC, col = I(col), key = key)) +
        geom_point(size = 0.45) + 
        theme_bw()  +
        xlab("log2(baseMean)") + 
        ylab("log2FC")
      
      plotly_object <- ggplotly(p,source = "source1") %>% 
        layout(dragmode = "lasso")  %>% 
        layout(showlegend = FALSE)
      
      config(plotly_object,
             toImageButtonOptions= list(filename = paste0("MAPlot_pvalue_", pvalue, "_log2FC_", tresholdlog2foldchange)))
      
      
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

    # Annotation Table
    
    annot <- eventReactive(input$Run_Annotation_ENSEMBL_to_GO, {
      data <- re() 
      orga = input$espece
      generef = bitr(data$ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=orga)
    })
    
    kegg_data <- reactive({
      df <- re()

      ids <- annot()
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
      list(kegg_gene_list=kegg_gene_list, df2=df2)
    })
   
    gse_kegg <- reactive({
      gse_kegg_data <- kegg_data()
      gse_kegg_annal <- gseKEGG(geneList     = gse_kegg_data$kegg_gene_list,
                          organism     = "mmu",
                          nPerm        = 10000,
                          minGSSize    = 3,
                          maxGSSize    = 800,
                          pvalueCutoff = 0.05,
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
      
      # filter on log2fold change (PARAMETER)
      kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 2]
      
      kk <- enrichKEGG(gene=kegg_genes, universe=names(kegg_gene_list),organism='mmu', pvalueCutoff = 0.05, keyType = "ncbi-geneid")
      
    })
    
    
    output$enrichKEGG_table <- renderDataTable({
      data <- ora_kegg()
      data.df <- as.data.frame(data)
      DT::datatable(data.df)
    })#fin renderDataTable
   
#     annotation_gsea <- eventReactive(input$GSEA_Annotation, {
#       pval_threshold <- input$pvalue_gsea
#       if (renderPrint({ input$method }) == 1){
#         
#         output$barplot_kegg <- renderPlotly({
#           result_gse_kegg <- ora_kegg()
#           barplot(result_gse_kegg, showCategory = 10)
#         })
#         
#         output$pathview_kegg <- renderImage({
#           result_gse_kegg <- ora_kegg()
#           
#           path_id<-result_gse_kegg$ID[1]
#           pathview(cpd.data=kegg_gene_list, pathway.id=path_id, species = "mmu")
#           path_img<-paste("./",path_id,".pathview.png", sep="")
#           list(src = path_img,
#                alt = "This is alternate text")
#         }, deleteFile = TRUE)
#       }
#       else {
#         
#         
#         output$gseaplot_kegg <- renderPlotly({
#           result_gse_kegg <- gse_kegg()
#           
#           gseaplot(result_gse_kegg, by = "all", title = kk2$Description[2], geneSetID = 1)
#         })  
#         
#         output$dotplot_kegg <- renderPlotly({
#           toto <- gse_kegg()
#           dotplot(toto, showCategory = 10, title = "Enriched Pathways" , split=".sign")
#         })
#       }
#     })
#   
#     
# ###########################A FAIRE FONCTIONNER################################  
#     output$gsea_annot <- renderDataTable({
#       D <- ora_kegg()
#       DT::datatable(D)
#     })#fin renderDataTable
# ###########################A FAIRE FONCTIONNER################################
    

     output$dotplot_kegg <- renderPlotly({
     toto <- gse_kegg()
     dotplot(toto, showCategory = 10, title = "Enriched Pathways" , split=".sign")
     })




     output$barplot_kegg <- renderPlotly({
       result_gse_kegg <- ora_kegg()
       barplot(result_gse_kegg, showCategory = 10)
     })


     output$gseaplot_kegg <- renderPlotly({
       result_gse_kegg <- gse_kegg()
       gseaplot(result_gse_kegg, by = "all", title =kk$Description[1], geneSetID = 1)
     })


     output$pathview_kegg <- renderImage({
       result_gse_kegg <- ora_kegg()

    
       path_id<-result_gse_kegg$ID[1]
       pathview(cpd.data=kegg_gene_list, pathway.id=path_id, species = "mmu")
       path_img<-paste("./",path_id,".pathview.png", sep="")
       list(src = path_img,
            alt = "This is alternate text")
     }, deleteFile = TRUE)


     # -------------------------------------------------------------------
     # BODY: tabPanel : Protein DOmains   --------------------------------
     # -------------------------------------------------------------------
     

     domain_enrichment <- reactive({
       
       pvalue <- pvalue()
       espece <- espece()
       
       # or
       ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
       

       # recuperation des domain ID pour les ensembl ID de notre jeu de donnees 
       interpro_id <- getBM(
         attributes=c('interpro', 'ensembl_gene_id'), # namespace_1003 = for go domain
         filters = 'ensembl_gene_id', #ensembl_gene_id
         values = data$ID, 
         mart = ensembl)
       
       #resultats de l'analyse differentielle:
       resOrdered <- re()
       
       #liste d'interet
       GeneList = resOrdered[which(resOrdered$padj<=0.05),]$ID
       GeneList = data.frame(GeneList)
       
       genes = resOrdered$ID
       interpro_id <- getBM(
         attributes=c('ensembl_gene_id', 'interpro'), # namespace_1003 = for go domain
         filters = 'ensembl_gene_id', #ensembl_gene_id
         values = genes, 
         mart = ensembl)
       
       GeneRef <- subset(interpro_id, interpro_id$interpro != "")
       
       
       ###  I - prepare data for enrichment      

       get_Gene_and_Bg_ratio = function(GeneList, GeneRef) {
         # reference list
         # m : nb of annotated genes in the reference list (for each term)
         m = table(GeneRef$interpro)
         m
         # n : nb of non annotated genes in the reference list
         n = length(unique(GeneRef$ensembl_gene_id)) - m
         n
         # experience (interest list)
         # x : nb of annotated genes in the interest list
         experience = merge(GeneList, GeneRef, by.x = "GeneList", by.y = "ensembl_gene_id")
         x = table(factor(experience$interpro, rownames(m)))
         x
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
       
       Gene.Bg.ratio = get_Gene_and_Bg_ratio(GeneList = GeneList, GeneRef = GeneRef)
       head(Gene.Bg.ratio)
       Bg.ratio = signif(100 * Gene.Bg.ratio$m/(Gene.Bg.ratio$m + Gene.Bg.ratio$n), 3)
       Gene.ratio = signif(100 * Gene.Bg.ratio$x / Gene.Bg.ratio$k, 3)
       
       
       ###  II -  hypergeometric test                        

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
       
       ggplot(data = D, aes(x = count, y = Term)) +
         geom_bar(stat = "identity", aes(fill = padj)) +
         theme(axis.text.x = element_text(
           angle = 90,
           hjust = 1,
           vjust = 0.5
         ))
       
     })
     
    
    } # end function(input, output) {
) # end shinyServer(
