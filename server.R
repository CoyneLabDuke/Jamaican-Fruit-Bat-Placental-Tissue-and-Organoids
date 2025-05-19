library(shiny)
library(hdf5r)
library(ggplot2)
library(reshape2)
library(viridisLite)
library(ggrepel)
library(heatmaply)
library(pheatmap)
library(DT)

sList <- c(Small = 18, Medium = 24, Large = 30)
pList <- c(Small = 5, Medium = 6, Large = 7)

datasets <- lapply(1:4, function(i) {
  prefix <- paste0("sc", i)
  list(
    meta = readRDS(paste0(prefix, "meta.rds")),
    conf = readRDS(paste0(prefix, "conf.rds")),
    gene = readRDS(paste0(prefix, "gene.rds")),
    def  = readRDS(paste0(prefix, "def.rds")),
    h5   = H5File$new(paste0(prefix, "gexpr.h5"), mode = "r")
  )
})
names(datasets) <- paste0("sc", 1:4)

bulk_rpkm <- read.csv("Bat_All_Samples_RPKMcorrect.csv", row.names = 1)
bulk_meta <- data.frame(Sample = colnames(bulk_rpkm))
bulk_meta$Group <- gsub(".*_(TO|DO|Fib).*", "\\1", bulk_meta$Sample)

shinyServer(function(input, output, session) {
  observe({
    updateSelectInput(session, "bulk_gene", choices = rownames(bulk_rpkm), selected = "GATA3")
  })
  
  for (scx in names(datasets)) {
    local({
      sc <- scx
      meta <- datasets[[sc]]$meta
      conf <- datasets[[sc]]$conf
      gene <- datasets[[sc]]$gene
      def  <- datasets[[sc]]$def
      h5   <- datasets[[sc]]$h5
      
      updateSelectInput(session, paste0(sc, "_meta"), choices = conf$UI, selected = def$meta1)
      updateSelectizeInput(session, paste0(sc, "_gene_input"), choices = names(gene), selected = def$gene1, server = FALSE)
      
      output[[paste0(sc, "_umap")]] <- renderPlot({
        dimX <- meta[[conf[conf$UI == def$dimred[1], "ID"][[1]]]]
        dimY <- meta[[conf[conf$UI == def$dimred[2], "ID"][[1]]]]
        val  <- meta[[conf[conf$UI == input[[paste0(sc, "_meta")]], "ID"][[1]]]]
        df <- data.frame(x = dimX, y = dimY, val = val)
        gg <- ggplot(df, aes(x, y, color = val)) +
          geom_point(size = 0.6) +
          theme_minimal(base_size = sList[[input[[paste0(sc, "_font_size")]]]]) +
          labs(x = def$dimred[1], y = def$dimred[2])
        if (is.numeric(df$val)) {
          gg + scale_color_viridis_c()
        } else {
          clrs <- strsplit(conf[conf$UI == input[[paste0(sc, "_meta")]], "fCL"][[1]], "\\|")[[1]]
          names(clrs) <- levels(factor(df$val))
          gg + scale_color_manual(values = clrs)
        }
      })
      
      output[[paste0(sc, "_gene")]] <- renderPlot({
        req(input[[paste0(sc, "_gene_input")]])
        gene_idx <- gene[[input[[paste0(sc, "_gene_input")]]]]
        dimX <- meta[[conf[conf$UI == def$dimred[1], "ID"][[1]]]]
        dimY <- meta[[conf[conf$UI == def$dimred[2], "ID"][[1]]]]
        expr_vals <- h5[["grp"]][["data"]]$read(args = list(gene_idx, quote(expr =)))
        df <- data.frame(x = dimX, y = dimY, val = pmax(expr_vals, 0))
        ggplot(df, aes(x = x, y = y, color = val)) +
          geom_point(size = 0.6) +
          theme_minimal(base_size = sList[[input[[paste0(sc, "_font_size")]]]]) +
          labs(x = def$dimred[1], y = def$dimred[2]) +
          scale_color_viridis_c()
      })
      
      output[[paste0(sc, "_vln")]] <- renderPlot({
        gene_idx <- gene[[input[[paste0(sc, "_gene_input")]]]]
        expr_vals <- h5[["grp"]][["data"]]$read(args = list(gene_idx, quote(expr =)))
        group <- factor(meta[[conf[conf$UI == def$grp1, "ID"][[1]]]])
        df <- data.frame(group = group, expr = pmax(expr_vals, 0))
        clrs <- strsplit(conf[conf$UI == def$grp1, "fCL"][[1]], "\\|")[[1]]
        names(clrs) <- levels(group)
        ggplot(df, aes(group, expr, fill = group)) +
          geom_violin() +
          scale_fill_manual(values = clrs) +
          theme_minimal(base_size = sList[[input[[paste0(sc, "_font_size")]]]]) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      })
      
      output[[paste0(sc, "_bubble")]] <- renderPlot({
        genes <- strsplit(input[[paste0(sc, "_multi_gene")]], "[,;\\s]+")[[1]]
        genes <- intersect(genes, names(gene))
        validate(need(length(genes) > 0, "Enter at least one valid gene."))
        
        group <- factor(meta[[conf[conf$UI == def$grp1, "ID"][[1]]]])
        expr_matrix <- sapply(genes, function(g) {
          g_idx <- gene[[g]]
          h5[["grp"]][["data"]]$read(args = list(g_idx, quote(expr =)))
        })
        
        if (is.null(dim(expr_matrix))) {
          expr_matrix <- matrix(expr_matrix, ncol = 1)
          colnames(expr_matrix) <- genes
        }
        
        expr_matrix <- pmax(expr_matrix, 0)
        df <- data.frame(Group = group, expr_matrix)
        df_long <- melt(df, id.vars = "Group", variable.name = "Gene", value.name = "Expression")
        
        avg_expr <- aggregate(Expression ~ Group + Gene, data = df_long, FUN = mean)
        pct_expr <- aggregate(Expression ~ Group + Gene, data = df_long, FUN = function(x) mean(x > 0) * 100)
        names(pct_expr)[3] <- "PercentExpressing"
        merged <- merge(avg_expr, pct_expr)
        
        ggplot(merged, aes(x = Gene, y = Group)) +
          geom_point(aes(size = PercentExpressing, color = Expression)) +
          scale_color_viridis_c() +
          theme_minimal(base_size = sList[[input[[paste0(sc, "_font_size")]]]]) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = "DotPlot: Mean Expression and % Expressing")
      })
      
      output[[paste0(sc, "_download_dot")]] <- downloadHandler(
        filename = function() paste0(sc, "_dotplot.pdf"),
        content = function(file) {
          pdf(file)
          genes <- strsplit(input[[paste0(sc, "_multi_gene")]], "[,;\\s]+")[[1]]
          genes <- intersect(genes, names(gene))
          if (length(genes) == 0) {
            plot.new(); text(0.5, 0.5, "No genes selected"); dev.off(); return()
          }
          
          group <- factor(meta[[conf[conf$UI == def$grp1, "ID"][[1]]]])
          expr_matrix <- sapply(genes, function(g) {
            g_idx <- gene[[g]]
            h5[["grp"]][["data"]]$read(args = list(g_idx, quote(expr =)))
          })
          expr_matrix <- pmax(expr_matrix, 0)
          df <- data.frame(Group = group, expr_matrix)
          df_long <- melt(df, id.vars = "Group", variable.name = "Gene", value.name = "Expression")
          
          avg_expr <- aggregate(Expression ~ Group + Gene, data = df_long, FUN = mean)
          pct_expr <- aggregate(Expression ~ Group + Gene, data = df_long, FUN = function(x) mean(x > 0) * 100)
          names(pct_expr)[3] <- "PercentExpressing"
          merged <- merge(avg_expr, pct_expr)
          
          print(
            ggplot(merged, aes(x = Gene, y = Group)) +
              geom_point(aes(size = PercentExpressing, color = Expression)) +
              scale_color_viridis_c() +
              theme_minimal(base_size = 14)
          )
          dev.off()
        }
      )
      
      output[[paste0(sc, "_heatmap")]] <- renderPlot({
        genes <- strsplit(input[[paste0(sc, "_multi_gene")]], "[,; ]+")[[1]]
        genes <- intersect(genes, names(gene))
        if (length(genes) == 0) return(NULL)
        
        groups <- factor(meta[[conf[conf$UI == def$grp1, "ID"][[1]]]])
        mat <- sapply(genes, function(g) h5[["grp"]][["data"]]$read(args = list(gene[[g]], quote(expr =))))
        mat <- pmax(mat, 0)
        avg_expr <- aggregate(mat, by = list(group = groups), FUN = mean)
        rownames(avg_expr) <- avg_expr$group
        avg_expr <- avg_expr[, -1, drop = FALSE]
        
        heatmap(as.matrix(avg_expr), Rowv = NA, Colv = NA,
                col = viridis(50), scale = "none", margins = c(6, 6))
      })
    })
  }
  
  deg_data <- reactive({
    req(input$deg_contrast)
    if (grepl("\\.csv$", input$deg_contrast)) {
      read.csv(input$deg_contrast, row.names = 1)
    } else {
      read.delim(input$deg_contrast, row.names = 1)
    }
  })
  
  output$volcano_plot <- renderPlot({
    deg <- deg_data()
    deg$Gene <- rownames(deg)
    deg$logP <- -log10(deg$padj)
    deg$Significant <- ifelse(
      deg$padj < input$padj_thresh & abs(deg$log2FoldChange) > input$fc_thresh,
      "Yes", "No"
    )
    
    highlight <- character(0)
    if (!is.null(input$volcano_highlight_genes) && nzchar(input$volcano_highlight_genes)) {
      highlight <- unlist(strsplit(input$volcano_highlight_genes, ",\\s*"))
      highlight <- intersect(highlight, deg$Gene)
    }
    
    gg <- ggplot(deg, aes(x = log2FoldChange, y = logP, color = Significant)) +
      geom_point(alpha = 0.5) +
      scale_color_manual(values = c("grey", "red")) +
      theme_minimal(base_size = 16) +
      labs(x = "log2 Fold Change", y = "-log10(padj)")
    
    if (length(highlight) > 0) {
      gg <- gg + ggrepel::geom_text_repel(
        data = subset(deg, Gene %in% highlight),
        aes(label = Gene),
        size = 4,
        max.overlaps = Inf
      )
    }
    
    gg
  })
  
  output$bulk_expr_plot <- renderPlot({
    req(input$bulk_gene)
    df <- data.frame(Expression = as.numeric(bulk_rpkm[input$bulk_gene, ]), Group = bulk_meta$Group)
    ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(na.rm = TRUE) +
      theme_minimal(base_size = 16) +
      labs(title = paste("Expression of", input$bulk_gene)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$bulk_deg_table <- DT::renderDT({
    deg <- deg_data()
    deg$Gene <- rownames(deg)
    deg <- subset(deg, padj < input$padj_thresh & abs(log2FoldChange) > input$fc_thresh)
    deg[, c("Gene", "log2FoldChange", "padj")]
  })
  
  output$bulk_heatmap <- renderPlotly({
    req(input$bulk_heat_genes)
    gene_list <- unlist(strsplit(input$bulk_heat_genes, ",\\s*"))
    gene_list <- intersect(gene_list, rownames(bulk_rpkm))
    if (length(gene_list) < 2) return(NULL)
    mat <- bulk_rpkm[gene_list, , drop = FALSE]
    rownames(mat) <- gene_list
    colnames(mat) <- gsub("_", "\n", colnames(mat))
    
    heatmaply::heatmaply(
      mat,
      scale = "row",
      colors = viridisLite::viridis(256),
      xlab = "Samples", ylab = "Genes",
      main = "Bulk Expression Heatmap",
      dendrogram = "both",
      hide_colorbar = FALSE
    )
  })
  
  output$download_bulk_heatmap <- downloadHandler(
    filename = function() "bulk_heatmap.pdf",
    content = function(file) {
      gene_list <- unlist(strsplit(input$bulk_heat_genes, ",\\s*"))
      gene_list <- intersect(gene_list, rownames(bulk_rpkm))
      if (length(gene_list) < 2) {
        pdf(file); plot.new(); text(0.5, 0.5, "Not enough genes to plot"); dev.off(); return()
      }
      mat <- bulk_rpkm[gene_list, , drop = FALSE]
      pheatmap::pheatmap(mat, scale = "row", fontsize = 10, filename = file)
    }
  )
})
