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
#BiocManager::install("AnnotationForge")

# #organism = "org.Hs.eg.db" 
# #BiocManager::install(organism, character.only = TRUE) 
# library(organism, character.only = TRUE) 
# library(org.Mm.eg.db)

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
library(ggrepel)
#library(wordcloud2)
library(wordcloud)

# Define server 
shinyServer(function(input, output, session) {

######################################################################### 
###### input data    
#########################################################################

ex_data <- reactive({
  data <- read.csv("www/exemple.csv")
  req(data)
  })
  
output$download_exemple <- downloadHandler(
  filename = function(){"exemple.csv"},
  content = function(fname){
    write.csv(ex_data(), fname)
  }
)
  
re <- reactive({
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      req(file)
      validate(need(ext == "csv", "Invalid file. Please upload a .csv file"))
      data <- read.csv(file$datapath, header = TRUE, sep = ";") 
      data <- na.omit(data)
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

espece_id <- reactive({
  espece_id <- input$espece_id
  req(espece_id)
})

biomart_dataset <- reactive({
  biomart_dataset <- input$biomart_dataset
  req(biomart_dataset)
})

orga_translate_table <- reactive({
  orga_translate_table  <- read.csv("www/orga_translate_table.csv", sep=";")
})

data_pour_plot <- eventReactive(input$Run_whole_data_inspection, {
  df_brut <- re()
  df <- na.omit(df_brut)
  pvalue <- pvalue()
  pvalue_log10 <- -log10(pvalue)
  tresholdlog2foldchange <- tresholdLog2FoldChange()
  
  df$key <- row.names(df)
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
  # add a column of NAs
  df$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  df$diffexpressed[df$log2FC > tresholdlog2foldchange & df$padj < pvalue] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  df$diffexpressed[df$log2FC < -tresholdlog2foldchange & df$padj < pvalue] <- "DOWN"
  
  req(df)
})

table_DEG_data <- reactive({
  D <- data_pour_plot()  
  pvalue <- pvalue()
  tresholdlog2foldchange <- tresholdLog2FoldChange()
  DEG = D[which(D$padj<=pvalue & abs(D$log2FC)> tresholdlog2foldchange ),]
  req(DEG)
})



    # BODY --------------------------------------------------------------------
    
        # BODY: tabPanel : whole Data Analysis --------------------------------

### Volcanoplot 

    output$volcanoPlot_plotly <- renderPlotly({
        #input data
        data_plot_plotly <- data_pour_plot()
        pvalue <- pvalue()
        tresholdlog2foldchange <- tresholdLog2FoldChange()

        # volcanoPlot using ggplot and plotly
        p <- ggplot(
          data = data_plot_plotly, 
          #mapping = aes(x = log2FC, y=-log10(padj), gene= GeneName, ID= ID, col=diffexpressed, key = key)) +
          mapping = aes(x = log2FC, y=-log10(padj), col=diffexpressed, key = key)) +
          geom_point(size = 0.45) + 
          theme_bw()  + 
          scale_color_manual(values=c("red", "black", "green")) +
          geom_hline(yintercept=-log10(pvalue), linetype="dashed", color = "black", size=0.1) +
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
      df <- data_pour_plot()
      pvalue <- pvalue()
      pvalue_log10 <- -log10(pvalue)
      tresholdlog2foldchange <- tresholdLog2FoldChange()
      
      # MAplot using ggplot and plotly
      p <- ggplot(
        data = df, 
        mapping = aes(
          #x = log2(baseMean), y = log2FC, gene= GeneName, ID= ID, col=diffexpressed, key = key)) +
          x = log2(baseMean), y = log2FC, col=diffexpressed, key = key)) +
        geom_point(size = 0.45)+ 
        theme_bw()  +
        scale_color_manual(values=c("red", "black", "green")) +
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
    
    selected_data <- reactive({
      D <- data_pour_plot() %>%  
        mutate(
          indice = row_number()
        )
      data_selected_from_graph <- event_data("plotly_selected",source = "source1")
      D_subset = subset(D, indice %in%  data_selected_from_graph$key)
    })

    output$Table_whole_data <- renderDataTable({ 
      D <- data_pour_plot()
      DT::datatable(D) 
    }) # fin renderDataTable({
    
    output$Table_DEG <- renderDataTable({ 
      DEG <- table_DEG_data()
      DT::datatable(DEG) 
    }) # fin renderDataTable({
    
    output$Table_data_selected <- renderDataTable({ 
      D <- selected_data()
      DT::datatable(D) 
    }) # fin renderDataTable({
    
    output$download_Table_DEG <- downloadHandler(
      filename = function(){"DEG.csv"}, 
      content = function(fname){
        write.csv(table_DEG_data(), fname)
      }
    )
    
    output$download_Table_data_selected <- downloadHandler(
      filename = function(){"selected_data.csv"}, 
      content = function(fname){
        write.csv(selected_data(), fname)
      }
    )
    
    #surligner la ligne correspondant au point clique sur le graphique
    proxy <- DT::dataTableProxy("Table_subset_data_selected")
    observe({
      data <- data_pour_plot()
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
         
    
    ####################################################
    ### GSEA pour GO
    #####################################################
    

    
    ####################################################
    ### ORA pour GO
    #####################################################

    
    # -------------------------------------------------------------------
    # BODY: tabPanel :KEGG --------------------------------
    # -------------------------------------------------------------------

    # Annotation Table
    
    kegg_data <- eventReactive(input$Run_Annotation_ENSEMBL_to_GO, {
      data <- re() 
      organism = espece()
      adj_method <- input$kegg_adj_method
      ids = bitr(data$ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
      
      df2 = data[data$ID %in% dedup_ids$ENSEMBL,]
      
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
      list(kegg_gene_list=kegg_gene_list, df2=df2, adj_method=adj_method)
    })
    
   
    gse_kegg <- reactive({
      gse_kegg_data <- kegg_data()
      gse_kegg_annal <- gseKEGG(geneList     = gse_kegg_data$kegg_gene_list,
                          organism     = "mmu",
                          nPerm        = 10000,
                          minGSSize    = 3,
                          maxGSSize    = 800,
                          pvalueCutoff = input$pvalue_gsea,
                          pAdjustMethod = gse_kegg_data$adj_method,
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
      # filter on log2fold change (under or over expressed DEG)
      if (input$type == 1){
        kegg_genes <- names(kegg_genes)[kegg_genes > 0]
      }
      else if(input$type == 2) {
        kegg_genes <- names(kegg_genes)[kegg_genes < 0]
      }
      else {
        kegg_genes <- names(kegg_genes)
      }
       kk <- enrichKEGG(gene=kegg_genes, universe=names(ora_kegg_data$kegg_gene_list),organism='mmu', pvalueCutoff = input$pvalue_gsea, keyType = "ncbi-geneid", pAdjustMethod = ora_kegg_data$adj_method)
    })
    
    output$enrichKEGG_table <- renderDataTable({
      if (input$method == 1){
        data <- ora_kegg()
      }
      else {
        data <- gse_kegg()
      }
      #updateSelectInput(session, "paths", choices = data$Description)
      updateSelectInput(session, "paths", choices = data$interpro_ID)
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
       kegg_data <- kegg_data()
       path_id<-result_kegg[result_kegg$Description==input$paths,]$ID

       # Make and save picture of pathway with pathview
       pathview(cpd.data=kegg_data$kegg_gene_list, pathway.id=path_id, species = "mmu")
       path_img<-paste("./",path_id,".png", sep="")
       # Return picture path to load it on popup window
       list(src = path_img)
     }, deleteFile = TRUE)


     # -------------------------------------------------------------------
     # BODY: tabPanel : Protein Domains   --------------------------------
     # -------------------------------------------------------------------
     
     ## recuperation des annotations
     Protein_Domains_data <- eventReactive(input$Run_protein_domains, {
       #input data
       resOrdered <- re()
       orga_translate_table <- orga_translate_table()
       espece_id <- espece_id()
       biomart_dataset <- orga_translate_table[espece_id,3]
       biomart_listMarts <- 'ensembl'
       
       # recuperation des domain ID pour les ensembl ID de notre jeu de donnees 
       ensembl =useMart(biomart=biomart_listMarts, 
                        dataset = biomart_dataset 
                        #host = "useast.ensembl.org"
                        )
       
       interpro_id_raw <- getBM(
         attributes=c('interpro', 'interpro_description', 'ensembl_transcript_id', 'ensembl_gene_id'),
         filters = 'ensembl_gene_id', 
         values = resOrdered$ID, 
         uniqueRows = FALSE,
         mart = ensembl)
       #enlever les valeurs manquantes
       interpro_id  <- subset(interpro_id_raw, interpro_id_raw['interpro'] != "")
     })
     
     comptage_Domains <- eventReactive(input$Run_protein_domains, {
       interpro_id <- Protein_Domains_data()
       # combien de fois ce domaine est present dans le set experimental, y compris lorsqu il est present plusieurs fois par proteine
       nb_domain = as.data.frame(table(interpro_id$interpro))
       colnames(nb_domain) <- c('interpro','nb_domain') 
       # combien de proteines du set experimental possedent au moins une fois un domaine X
       prt_par_domain = as.data.frame(table(unique(interpro_id[c('ensembl_transcript_id','interpro')])$interpro))
       colnames(prt_par_domain) <- c('interpro','prt_par_domain') 
       comptage_domains = merge(prt_par_domain, nb_domain, by.x = "interpro", by.y = "interpro")
     })
     
     DEG <-  eventReactive(input$Run_protein_domains, {
       resOrdered <- re()
       tresholdLog2FoldChange <- tresholdLog2FoldChange()
       pvalue <- pvalue()
       
       if(input$type_prt_domain == "Both"){
         GeneList = resOrdered[which(resOrdered$padj<=pvalue & abs(resOrdered$log2FC)>tresholdLog2FoldChange),]$ID
       }
       else if(input$type_prt_domain == "Under"){
         GeneList = resOrdered[which(resOrdered$padj<=pvalue & resOrdered$log2FC < tresholdLog2FoldChange),]$ID
       }
       else if(input$type_prt_domain == "Over"){
         GeneList = resOrdered[which(resOrdered$padj<=pvalue & resOrdered$log2FC > -tresholdLog2FoldChange),]$ID
       }
       GeneList
       
     })
     

     ## code fonction ORA 
     domain_enrichment_ORA <-  eventReactive(input$Run_protein_domains, if(input$method_prt_domain == 1){
       #input data
       resOrdered <- re()
       pvalue <- input$pvalue_prt_domain
       adjustment_method <- input$pvalue_adjustment_prt_domain
      
       # recuperation des domain ID pour les ensembl ID de notre jeu de donnees 
       interpro_id <- Protein_Domains_data()
       comptage_Domains <- comptage_Domains()
       #background
       GeneRef = unique(na.omit(interpro_id[c('interpro','ensembl_gene_id')]))
       summary(GeneRef)
       
       #liste d'interet
       GeneList <- DEG()
       GeneList <- data.frame(GeneList)
       
       ###  I - prepare data for enrichment      
       
       get_Gene_and_Bg_ratio = function(GeneList, GeneRef) {
         # experience (interest list)
         experience = merge(GeneList, GeneRef, by.x = "GeneList", by.y = "ensembl_gene_id")
         # x : nb of annotated genes in the interest list
         x = table(experience$interpro)
         
         # reference list
         #interpro en commun entre background et liste d'intéret
         #il faut créer un liste avec que les id de la liste d'interret, pas la peine de compter les autres
         #interpro_liste_int = merge(experience[2:3], GeneRef[1:2], by.x = "interpro", by.y = "interpro")
         interpro_liste_int = subset(GeneRef, GeneRef$interpro %in% experience$interpro)
         # m : nb of annotated genes in the reference list (for each term)
         m = table(interpro_liste_int$interpro)
         
         # n : nb of non annotated genes in the reference list
         n = length(unique(GeneRef$ensembl_gene_id))
         
         # k : total nb of genes in the interest list
         k = length(GeneList$GeneList)
         
         Term = unique(row.names(m))
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
       
       
       #  II -  hypergeometric test                        
       
       hypergeom_test = function(x, k, m, n){
         # calculate p-value and adjusted p-value
         pvalue_fisher.test <- c()
         # Test for over-representation (enrichment)
         # phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)
         # fisher.test(matrix(c(Overlap, group2-Overlap, group1-Overlap, Total-group2-group1 +Overlap), 2, 2), alternative='greater')$p.value
         for (i in 1:length(x)) {
             plop = ((fisher.test(matrix(c(x[i], m[i]-x[i], k-x[i], n-(k-x[i])),2,2), alternative='greater'))$p.value)
             pvalue_fisher.test <- c(pvalue_fisher.test, plop)
             }
         
         padj = p.adjust(pvalue_fisher.test, method = adjustment_method, n=length((pvalue_fisher.test)))
         return (list(pvalue_fisher.test = pvalue_fisher.test, padj = padj))
       }
       
       # III - create results table                        
       
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
         table.enrich = data.frame(interpro_ID = Gene.Bg.ratio$Term, 
                                   pvalue = signif(test$pvalue_fisher.test, 3),
                                   padj = signif(test$padj, 3),
                                   BgRatio = Bg.ratio,
                                   BgRatio_count = paste0(" (", Gene.Bg.ratio$m, "/", Gene.Bg.ratio$n, ")"),
                                   GeneRatio = Gene.ratio,
                                   GeneRatio_count = paste0(" (", Gene.Bg.ratio$x, "/", Gene.Bg.ratio$k, ")"),
                                   count = Gene.Bg.ratio$x, 
                                   count_and_padj = paste0("count: ", Gene.Bg.ratio$x, "\n padj: ", signif(test$padj, 3) )
                                   )
         
         return (table.enrich[order(table.enrich$pval), ])
       }
       res.enrich.hypergeom.prt_dommain = create_table_enrichment(GeneList = GeneList, GeneRef = GeneRef)
       result = res.enrich.hypergeom.prt_dommain[which(res.enrich.hypergeom.prt_dommain$padj<pvalue),]
       updateSliderInput(session, "nb_barplot_ora_coder", max = nrow(result))
       
       res_with_description = merge(result, unique(interpro_id[c('interpro_description','interpro')]), by.x = "interpro_ID", by.y = "interpro")
       res_with_domains_comptage = merge(res_with_description, comptage_Domains, by.x = "interpro_ID", by.y = "interpro")
       } )  
     
     ## avec la fonction enricher (a titre de comparaison)
     domain_enrichment_ORA_enricher <-  eventReactive(input$Run_protein_domains, if(input$method_prt_domain == 1){
       #input data
       resOrdered <- re()
       pvalue <- input$pvalue_prt_domain
       adjustment_method <- input$pvalue_adjustment_prt_domain
       # recuperation des domain ID pour les ensembl ID de notre jeu de donnees 
       interpro_id <- Protein_Domains_data()
       # recuperation des domain ID pour les ensembl ID de notre jeu de donnees 

       table_TERM2GENE = unique(na.omit(interpro_id[c('interpro','ensembl_gene_id')]))
       #gene_list_enricher = resOrdered[which(resOrdered$padj<=pvalue),]$ID
       DEG_under = DEG()
       result_enricher <- enricher(
         gene = unique(DEG_under), 
         pvalueCutoff = pvalue, 
         pAdjustMethod = adjustment_method, 
         universe = unique(resOrdered$ID), 
         minGSSize = 10, 
         maxGSSize = 500, 
         qvalueCutoff = 0.2, 
         TERM2GENE =table_TERM2GENE)
       })
     

     domain_enrichment_GSEA <-  eventReactive(input$Run_protein_domains, if(input$method_prt_domain == 2){
       #input data
       resOrdered <- re()
       pvalue <- input$pvalue_prt_domain
       adjustment_method <- input$pvalue_adjustment_prt_domain
       
       # recuperation des domain ID pour les ensembl ID de notre jeu de donnees 
       interpro_id <- Protein_Domains_data()
       table_TERM2GENE = unique(na.omit(interpro_id[c('interpro','ensembl_gene_id')]))
       
       # we want the log2 fold change 
       original_gene_list <- resOrdered$log2FC
       # name the vector
       names(original_gene_list) <- resOrdered$ID
       # omit any NA values 
       gene_list<-na.omit(original_gene_list)
       # sort the list in decreasing order (required for clusterProfiler)
       gene_list = sort(gene_list, decreasing = TRUE)
       
       result_GSEA = GSEA(
         geneList =gene_list,
         exponent = 1,
         minGSSize = 10,
         maxGSSize = 500,
         pvalueCutoff = pvalue,
         pAdjustMethod = "BH",
         TERM2GENE=table_TERM2GENE)
     })
     
     output$Table_domains_enrichment <- DT::renderDataTable(DT::datatable({
       data <- domain_enrichment_ORA() %>%
         mutate(interpro_link = paste0("<a href='https://www.ebi.ac.uk/interpro/entry/InterPro/", interpro_ID,"' target='_blank'>", interpro_ID,"</a>"))
       col_a_afficher = c('interpro_link','interpro_description', 'pvalue', 'padj', 'BgRatio', 'BgRatio_count', 'GeneRatio', 'GeneRatio_count', 'nb_domain', 'prt_par_domain')
       data[col_a_afficher]
     },
     escape = FALSE))
     
     output$barplot_domains_enrichment <- renderPlotly( if(input$method_prt_domain == 1){
         max <- input$nb_barplot_ora_coder
         D <- domain_enrichment_ORA()
         ggplot(data = D[1:max,], aes(x= count , y = reorder(interpro_ID, count)  ) ) +
           geom_bar(stat = "identity", aes(fill = padj))  +
           scale_fill_viridis_c(option = 'magma') +
           theme(axis.text.x = element_text(
             angle = 90,
             hjust = 1,
             vjust = 0.5
           )) +
           labs(y = "Protein Domain (Interpro ID)", x = "Count")

     })

     output$dotplot_domains_enrichment <- renderPlotly( if(input$method_prt_domain == 1){
         max <- input$nb_barplot_ora_coder
         D <- domain_enrichment_ORA()

         ggplot(D[1:max,], aes(x=interpro_ID, y=GeneRatio)) +
           geom_point(stat='identity', aes(col=padj, size=count), alpha=0.75) +   # Draw points
           scale_colour_viridis_c(option = 'magma') +
           geom_segment(aes(x=interpro_ID,
                            xend=interpro_ID,
                            y=min(GeneRatio),
                            yend=max(GeneRatio)
                            ),
                        linetype="none",
                        size=0.1) +   # Draw dashed lines
           labs(x="Protein Domain",
                y="Gene Ratio",
                title = "Dotplot",
                fill="padj") +
           coord_flip()

     })
     
     
     output$piechart_domains_enrichment <- renderPlotly( if(input$method_prt_domain == 1){
       D <- domain_enrichment_ORA()
       data <- D[,c('interpro_ID', 'count')]
       
       fig <- plot_ly(D, 
                      labels = ~interpro_ID, 
                      values = ~count, 
                      type = 'pie',
                      text = ~count_and_padj,
                      hovertemplate = paste("%{text}")
                      )
       fig <- fig %>% layout(title = 'Pie chart ',
                             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

       
       fig
       
     })
     
     output$cloud <- renderPlot({
       D <- domain_enrichment_ORA()
       data <- D[,c('interpro_ID', 'count')]
       
       wordcloud(data$interpro_ID, data$count, scale=c(2,0.1),
                     colors=brewer.pal(8, "Dark2"))
     })
     
     output$Table_domains_enrichment_enricher <- renderDataTable( if(input$method_prt_domain == 1){ 
       result_enricher <- domain_enrichment_ORA_enricher()
       D <- as.data.frame(result_enricher)[c(1, 3:7,9)]
       DT::datatable(D) 
     }) # fin renderDataTable({
     
     
     output$dotplot_domains_enrichment_enricher <- renderPlotly( if(input$method_prt_domain == 1){ 
       result_ORA_enricher <- domain_enrichment_ORA_enricher()
       dotplot(result_ORA_enricher)
     }) # fin renderDataTable({
     
     output$barplot_domains_enrichment_enricher <- renderPlotly( if(input$method_prt_domain == 1){ 
       result_ORA_enricher <- domain_enrichment_ORA_enricher()
       barplot(result_ORA_enricher)
     }) # fin renderDataTable({
     
     
     output$Table_domain_enrichment_GSEA <- renderDataTable( if(input$method_prt_domain == 2){ 
       result_GSEA <- domain_enrichment_GSEA()
       D <- as.data.frame(result_GSEA)[c(1, 3:7,9)]
       
       updateSelectInput(session, "protein_id_list", choices = D$ID)
       DT::datatable(D)
     }) # fin renderDataTable({
     
     output$barplot_domain_enrichment_GSEA <- renderPlotly( if(input$method_prt_domain == 2){ 
       result_GSEA <- domain_enrichment_GSEA()
       #dotplot(result_GSEA)
       dotplot(result_GSEA, showCategory=10, split=".sign") + facet_grid(.~.sign)
       }) # fin renderDataTable({
    
     output$gseaplot_domain_enrichment_GSEA <- renderPlotly( if(input$method_prt_domain == 2){ 
       result_GSEA <- domain_enrichment_GSEA()
       gseaplot(result_GSEA, input$protein_id_list)
     }) # fin renderDataTable({
     
     
     
    } # end function(input, output) {
) # end shinyServer(
