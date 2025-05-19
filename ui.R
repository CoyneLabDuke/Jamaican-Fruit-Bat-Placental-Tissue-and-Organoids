library(shiny)
library(bslib)
library(viridisLite)
library(DT)        # for DTOutput()
library(plotly) 

# Panel loader for snRNA-seq datasets
load_panel <- function(prefix, label) {
  tabPanel(
    title = label,
    sidebarLayout(
      sidebarPanel(
        selectInput(paste0(prefix, "_meta"), "UMAP color by (metadata):", choices = NULL),
        selectizeInput(paste0(prefix, "_gene_input"), "Gene (UMAP + violin):", choices = NULL,
                       options = list(placeholder = "Type or paste gene name")),
        textInput(paste0(prefix, "_multi_gene"), "Genes (comma-separated for DotPlot/heatmap):", value = ""),
        radioButtons(paste0(prefix, "_plot_size"), "Plot size:", c("Small", "Medium", "Large"), selected = "Medium", inline = TRUE),
        radioButtons(paste0(prefix, "_font_size"), "Font size:", c("Small", "Medium", "Large"), selected = "Medium", inline = TRUE)
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("UMAP (metadata)", downloadButton(paste0(prefix, "_download_umap"), "Download"),
                   plotOutput(paste0(prefix, "_umap"), height = "600px")),
          tabPanel("UMAP (gene)", downloadButton(paste0(prefix, "_download_gene"), "Download"),
                   plotOutput(paste0(prefix, "_gene"), height = "600px")),
          tabPanel("Violin Plot", downloadButton(paste0(prefix, "_download_vln"), "Download"),
                   plotOutput(paste0(prefix, "_vln"), height = "600px")),
          tabPanel("DotPlot", downloadButton(paste0(prefix, "_download_dot"), "Download"),
                   plotOutput(paste0(prefix, "_bubble"), height = "600px")),
          tabPanel("Heatmap", downloadButton(paste0(prefix, "_download_heatmap"), "Download"),
                   plotOutput(paste0(prefix, "_heatmap"), height = "600px"))
        )
      )
    )
  )
}

shinyUI(
  fluidPage(
    theme = bs_theme(bootswatch = "flatly", version = 5),
    titlePanel("Jamaican Fruit Bat Placental Atlas"),
    
    navbarPage(
      title = NULL,
      
      # Home tab
      tabPanel("Home",
               fluidRow(
                 column(12,
                        h3("Explore the Jamaican Fruit Bat Placenta: A Multi-Modal Gene Expression Atlas"),
                        img(src = "bat_placenta.png", height = "300px", style = "display:block; margin:auto;"),
                        br(),
                        h4("About This App"),
                        p("This interactive atlas showcases the cellular and molecular landscape of the Jamaican fruit bat (Artibeus jamaicensis) placenta, using cutting-edge single-nucleus and bulk RNA sequencing. Explore how different cell types form, function, and interact during pregnancy. Compare bat, human, and mouse placental cells, and see how the bat placenta responds to immune challenges like viral mimicry (poly I:C). This resource offers a unique window into placental biology with implications for both reproductive health and species conservation."),
                        br(),
                        
                        h4("Available Datasets"),
                        tags$ul(
                          tags$li(strong("snRNA-seq:"), " Jfb placenta, trophoblast organoids (JfbTOs), and in vivo trophoblast lineages"),
                          tags$li(strong("Cross-species integration:"), " Bat, human, and mouse placenta-derived cell types"),
                          tags$li(strong("Bulk RNA-seq:"), " JfbTO, JfbDO, and JfbFib samples"),
                          tags$li(strong("Immune stimulation:"), " pIC vs Mock treatment in JfbTOs and human trophoblasts")
                        ),
                        br(),
                        
                        h4("How to Use the App"),
                        tags$ul(
                          tags$li(strong("UMAP (metadata):"), " Color cells by cluster, sample, or annotation"),
                          tags$li(strong("UMAP (gene):"), " Visualize individual gene expression in snRNA-seq space"),
                          tags$li(strong("Violin Plot:"), " Compare expression of a single gene across cell groups"),
                          tags$li(strong("DotPlot / Heatmap:"), " View expression of multiple genes across clusters"),
                          tags$li(strong("Bulk Boxplot:"), " Compare gene expression in TO, DO, or Fib samples"),
                          tags$li(strong("Volcano Plot:"), " Explore differentially expressed genes across 5 key comparisons"),
                          tags$li(strong("Bulk Heatmap:"), " Enter multiple genes to visualize expression trends across all bulk samples"),
                          tags$li(strong("All plots can be downloaded as PDFs"))
                        ),
                        br(),
                        
                        h4("Citation & Contact"),
                        p("For questions or feedback, please contact: ", a("carolyn.coyne@duke.edu", href = "mailto:carolyn.coyne@duke.edu")),
                        p(strong("Citation:"), " Caldwell, Yang, Casazza et al., Cellular and Immune Adaptations at the Maternal–Fetal Interface in Bats, 2025")
                 )
               )
      ),
      
      # snRNA-seq tabs
      load_panel("sc1", "Jfb Placental Tissue"),
      load_panel("sc2", "Trophoblast Lineages"),
      load_panel("sc3", "Jfb Organoid Model"),
      load_panel("sc4", "Cross-Species Atlas"),
      
      # Bulk RNA-seq tab
      tabPanel("Bulk RNA-seq",
               sidebarLayout(
                 sidebarPanel(
                   selectInput("bulk_gene", "Select gene for boxplot:", choices = NULL),
                   selectInput("deg_contrast", "DEG Comparison:",
                               choices = c(
                                 "JfbTO vs JfbDO" = "res_TO_v_DO_DESeq2.txt",
                                 "JfbTO vs JfbFib" = "res_TO_v_Fib_DESeq2.txt",
                                 "JfbDO vs JfbFib" = "res_DO_v_Fib_DESeq2.txt",
                                 "JfbTO: pIC vs Mock" = "DESeq2_results_pIC_vs_Mock_Jfb.csv",
                                 "HumanTO: pIC vs Mock" = "DESeq2_results_pIC_vs_Mock_Human.csv"
                               )
                   ),
                   sliderInput("padj_thresh", "FDR (padj) threshold:", min = 1e-10, max = 0.1, value = 0.05, step = 0.005),
                   sliderInput("fc_thresh", "Absolute log2 Fold Change threshold:", min = 0, max = 5, value = 1, step = 0.1),
                   textInput("bulk_heat_genes", "Genes for heatmap (comma-separated):", value = ""),
                   textInput("volcano_highlight_genes", "Highlight genes on volcano (comma-separated):", value = ""),
                   downloadButton("download_bulk_plot", "Download Boxplot"),
                   downloadButton("download_bulk_heatmap", "Download Heatmap (PDF)"),
                   downloadButton("download_volcano", "Download Volcano Plot")
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Boxplot", plotOutput("bulk_expr_plot", height = "400px")),
                     tabPanel("Volcano Plot", plotOutput("volcano_plot", height = "400px")),
                     tabPanel("DEG Table", DT::DTOutput("bulk_deg_table")),
                     tabPanel("Bulk Heatmap", plotlyOutput("bulk_heatmap", height = "600px"))
                   )
                 )
               )
      )
    )
  )
)
