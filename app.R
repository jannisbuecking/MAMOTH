# --- 0. Load Shiny and Other Libraries ---
# Ensure all these packages are installed first
library(shiny)
library(shinyjs)
library(openxlsx)
library(shinythemes)
library(DESeq2)
library(bslib)                
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)                
library(stringr)
library(forcats)
library(patchwork)
library(SummarizedExperiment)
library(ggpubr)
library(org.Hs.eg.db)
library(annotate)
library(AnnotationHub)
library(grid)
library(Homo.sapiens)
library(clusterProfiler)
library(MOFA2)              # For MOFA
library(tibble)               # For MOFA
library(viridis)              # For MOFA
library(GGally)               # For MOFA
library(RColorBrewer)         # For MOFA
library(rlang)                # For MOFA
library(gprofiler2)           # For MOFA GSEA
library(cowplot)              # For MOFA GSEA
library(pheatmap)             # For MOFA Heatmaps
library(arrow)

mammoth_icon_url <- "https://cdn-icons-png.flaticon.com/512/1888/1888595.png"

# --- 1. Load Data (Executed once when the app starts) ---
data_loaded_flag <- FALSE
rna_data_prepared <- FALSE
mofa_data_loaded <- FALSE
available_gois <- c("Loading...")
master_condition_order <- NULL
available_conditions <- c()

tryCatch({
  message("Loading data from fast .arrow files...")
  
  # --- Helper function to reassemble SummarizedExperiment objects ---
  reassemble_se <- function(name_prefix) {
    row_data_df <- read_feather(file.path("data_arrow", paste0(name_prefix, "_rowData.arrow")))
    col_data_df <- read_feather(file.path("data_arrow", paste0(name_prefix, "_colData.arrow")))
    assay_df <- read_feather(file.path("data_arrow", paste0(name_prefix, "_assay.arrow")))
    
    # Restore rownames from the saved columns
    row_data_df <- column_to_rownames(row_data_df, "feature_id")
    col_data_df <- column_to_rownames(col_data_df, "sample_id")
    assay_matrix <- as.matrix(column_to_rownames(assay_df, "feature_id"))
    
    SummarizedExperiment(
      assays = list(counts = assay_matrix),
      rowData = row_data_df,
      colData = col_data_df
    )
  }
  
  # --- Reassemble all main data objects ---
  # Result objects (for p-values and volcano plots)
  dep_8330 <- reassemble_se("dep_8330")
  dep_MGH <- reassemble_se("dep_MGH")
  Ubi_dep_8330 <- reassemble_se("ubi_dep_8330")
  Ubi_dep_MGH <- reassemble_se("ubi_dep_MGH")
  
  # --- ADDED THESE FOUR LINES TO LOAD THE _imp_ OBJECTS ---
  # Raw data objects (for gene expression plots)
  new_data_imp_8330 <- reassemble_se("new_data_imp_8330")
  new_data_imp_MGH <- reassemble_se("new_data_imp_MGH")
  Ubi_data_imp_8330 <- reassemble_se("Ubi_data_imp_8330")
  Ubi_data_imp_MGH <- reassemble_se("Ubi_data_imp_MGH")
  
  # --- Reassemble the main DESeqDataSet object ---
  counts_df_rna <- read_feather("data_arrow/dds_RNA_assay.arrow")
  col_data_df_rna <- read_feather("data_arrow/dds_RNA_colData.arrow")
  row_data_df_rna <- read_feather("data_arrow/dds_RNA_rowData.arrow")
  
  counts_matrix_rna <- as.matrix(column_to_rownames(counts_df_rna, "feature_id"))
  col_data_rna <- column_to_rownames(col_data_df_rna, "sample_id")
  row_data_rna <- column_to_rownames(row_data_df_rna, "feature_id")
  
  dds_RNA <- DESeqDataSetFromMatrix(
    countData = counts_matrix_rna, colData = col_data_rna, design = ~ 1
  )
  
  # --- THIS WAS THE BUGGY LINE - NOW CORRECTED ---
  # We use rowData()<- instead of the incorrect mcols()<- approach
  rowData(dds_RNA) <- DataFrame(row_data_rna)
  
  # --- Load smaller list and MOFA objects from consolidated file ---
  load("secondary_data.RData")
  data_loaded_flag <- TRUE
  
  if (exists("dds_RNA") && is(dds_RNA, "DESeqDataSet")) {
    available_gois <- sort(unique(rownames(dds_RNA)))
  } else { available_gois <- c("Transcriptome Data Error"); data_loaded_flag <- FALSE }
  
  if (exists("dds_RNA") && is(dds_RNA, "DESeqDataSet") && "background" %in% colnames(colData(dds_RNA))) {
    dds_RNA_8330 <- dds_RNA[, dds_RNA$background == "8330"]; dds_RNA_MGH <- dds_RNA[, dds_RNA$background == "MGH"]
    colData(dds_RNA_8330)$condition <- colData(dds_RNA_8330)$genotype; colData(dds_RNA_MGH)$condition <- colData(dds_RNA_MGH)$genotype
    master_condition_order <- levels(colData(dds_RNA_8330)$condition)
    available_conditions <- master_condition_order
    rna_data_prepared <- TRUE
  }
  
  # MOFA metadata wrangling
  MOFA_out@samples_metadata <- MOFA_out@samples_metadata %>%
    mutate(
      Background = str_extract(sample, "(?<=_)(8330|MGH)"), Background = factor(Background),
      MutationPrefix = str_extract(sample, "^[A-Za-z]+"),
      Mutation = case_when(
        MutationPrefix == "WT" ~ "WT", MutationPrefix == "DelC" ~ "DelC",
        MutationPrefix == "DupC" ~ "DupC", MutationPrefix == "HetDel" ~ "HetDel",
        MutationPrefix == "HomoDel" ~ "HomoDel", MutationPrefix == "PWSTI" ~ "PWSTI",
        MutationPrefix == "PWSTII" ~ "PWSTII", TRUE ~ MutationPrefix
      ),
      Mutation = factor(Mutation, levels = c("DelC", "DupC", "HetDel", "HomoDel", "PWSTI", "PWSTII", "WT")),
      Mutation_BG = paste0(Mutation, "_", Background)
    ) %>% mutate(group = Mutation) %>% dplyr::select(-MutationPrefix)
  if (!identical(rownames(MOFA_out@samples_metadata), MOFA_out@samples_metadata$sample)) {
    rownames(MOFA_out@samples_metadata) <- MOFA_out@samples_metadata$sample
  }
  mofa_data_loaded <- TRUE
  message("All data loaded and reassembled successfully.")
  
}, error = function(e) {
  warning("Error during data loading: ", e$message)
})

# --- 2. Global Objects & Helper Functions ---
# MOFA Globals
N_FACTORS <- if(exists("MOFA_out")) tryCatch(get_dimensions(MOFA_out)$K, error = function(e) { 15 }) else 15
VIEWS <- if(exists("MOFA_out")) tryCatch(views_names(MOFA_out), error = function(e) { c("RNA", "protein", "Ubiq") }) else c("RNA", "protein", "Ubiq")
view_choices <- setNames(VIEWS, case_when(VIEWS == "RNA" ~ "Transcriptome", VIEWS == "protein" ~ "Proteome", VIEWS == "Ubiq" ~ "Ubiquitome", TRUE ~ VIEWS))
gsea_view_choices <- c("Transcriptome" = "RNA", "Proteome" = "protein", "Ubiquitome" = "Ubiq")
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(100)
font_color_dark <- "grey20"
view_labeller <- c(`protein` = "Proteome", `RNA` = "Transcriptome", `Ubiq` = "Ubiquitome")

# General Helper Functions
format_feature_names <- function(feature_vector) { gsub("_RNA|_protein|_ubiq|_Ubiq", "", feature_vector) }
gg_color_hue <- function(n) { hues = seq(15, 375, length = n + 1); hcl(h = hues, l = 65, c = 100)[1:n] }
get_significance_label <- function(p_adj) { if (is.na(p_adj)) { return("") } else if (p_adj <= 0.0001) { return("****") } else if (p_adj <= 0.001) { return("***") } else if (p_adj <= 0.01) { return("**") } else if (p_adj <= 0.05) { return("*") } else { return("ns") } }
custom_density <- function(data, mapping, ...) { ggplot(data = data, mapping = mapping) + geom_density(aes(color = !!mapping$colour, fill = !!mapping$colour), alpha = 0.5) + theme(panel.grid = element_blank()) }
plotAbundance <- function(data, GOI, dep_results, ref_group = "WT",
                          plot_color = "grey80", data_source_label = "",
                          master_order,
                          star_vjust = -0.5, star_size = 6) {
  
  if (!GOI %in% rownames(data)) {
    return(ggplot() + theme_void() + annotate("text", x=0.5, y=0.5, label="Gene not found\nin this dataset") + ggtitle(GOI) + labs(subtitle=data_source_label))
  }
  
  if (inherits(data, "DESeqDataSet")) {
    if (is.null(sizeFactors(data))) data <- estimateSizeFactors(data)
    expr_values <- counts(data, normalized = TRUE)[GOI, ]
    y_axis_label <- "Normalized Counts"
  } else if (is(data, "SummarizedExperiment")) {
    expr_values <- assay(data)[GOI, ]
    y_axis_label <- "Normalized log2 Intensity"
  }
  
  plot_df <- data.frame(expression = expr_values, condition = colData(data)$condition)
  plot_df$condition <- factor(plot_df$condition, levels = master_order)
  other_conditions <- setdiff(master_order, ref_group)
  
  annotation_info <- list()
  y_offset <- diff(range(plot_df$expression, na.rm = TRUE)) * 0.05
  if (!is.finite(y_offset) || y_offset == 0) y_offset <- max(abs(plot_df$expression), na.rm=TRUE) * 0.05
  
  if (is(dep_results, "SummarizedExperiment")) {
    results_metadata <- as.data.frame(elementMetadata(dep_results))
    goi_idx <- which(results_metadata$name == GOI)[1]
    if (length(goi_idx) > 0 && !is.na(goi_idx)) {
      for (cond in other_conditions) {
        p_col_name <- paste0(cond, "_vs_", ref_group, "_p.adj")
        if (p_col_name %in% names(results_metadata)) {
          p_adj <- results_metadata[goi_idx, p_col_name]
          signif_label <- get_significance_label(p_adj)
          if (signif_label != "ns" && signif_label != "") {
            max_y <- max(plot_df$expression[plot_df$condition == cond], na.rm = TRUE)
            annotation_info[[cond]] <- data.frame(condition = cond, y_pos = max_y + y_offset, label = signif_label)
          }
        }
      }
    }
  } else if (is(dep_results, "list")) {
    for (cond in other_conditions) {
      res_name_fwd <- paste0(cond, "_vs_", ref_group)
      res_name_rev <- paste0(ref_group, "_vs_", cond)
      res_name <- NULL
      if (res_name_fwd %in% names(dep_results)) {
        res_name <- res_name_fwd
      } else if (res_name_rev %in% names(dep_results)) {
        res_name <- res_name_rev
      }
      if (!is.null(res_name)) {
        res_obj <- dep_results[[res_name]]
        if (is(res_obj, "DESeqResults") && GOI %in% rownames(res_obj)) {
          p_adj <- res_obj[GOI, "padj"]
          signif_label <- get_significance_label(p_adj)
          if (signif_label != "ns" && signif_label != "") {
            max_y <- max(plot_df$expression[plot_df$condition == cond], na.rm = TRUE)
            annotation_info[[cond]] <- data.frame(condition = cond, y_pos = max_y + y_offset, label = signif_label)
          }
        }
      }
    }
  }
  annotation_df <- bind_rows(annotation_info)
  
  # --- PLOTTING LOGIC UPDATED HERE ---
  p <- ggplot(plot_df, aes(x = condition, y = expression)) +
    # Layer 1: Individual points
    geom_jitter(color = plot_color, width = 0.2, size = 3, alpha = 0.7) +
    
    # Layer 2: Mean and Standard Deviation error bars
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                 geom = "errorbar", width = 0.2, color = plot_color) + # Changed color
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = plot_color) + # Changed color
    
    # --- END OF PLOTTING LOGIC UPDATE ---
    
    scale_x_discrete(drop = FALSE) +
    labs(title = GOI, subtitle = data_source_label, y = y_axis_label, x = NULL) +
    theme_bw(base_size = 16) +
    theme(
      text = element_text(color = "#545252"), axis.text = element_text(color = "#545252"),
      axis.ticks = element_line(color = "#545252"), axis.title.y = element_text(size = 12),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(t = 20, r = 5, b = 5, l = 5, "pt")
    )
  
  if (nrow(annotation_df) > 0) {
    p <- p + geom_text(data = annotation_df, aes(x = condition, y = y_pos, label = label), inherit.aes = FALSE, size = star_size, vjust = star_vjust)
  }
  
  min_y <- min(plot_df$expression, na.rm = TRUE)
  max_y_data <- max(plot_df$expression, na.rm = TRUE)
  max_y_annot <- if (nrow(annotation_df) > 0) max(annotation_df$y_pos, na.rm = TRUE) else -Inf
  highest_point <- max(max_y_data, max_y_annot, na.rm = TRUE)
  final_upper_limit <- highest_point + (diff(range(plot_df$expression, na.rm=TRUE)) * 0.15)
  if(is.finite(min_y) && is.finite(final_upper_limit)){
    p <- p + coord_cartesian(ylim = c(min_y, final_upper_limit), clip = "off")
  }
  
  return(p)
}

# --- 3. UI DEFINITION ---
ui <- fluidPage(
  useShinyjs(), # Initialize shinyjs
  
  # This div contains the login screen UI
  div(
    id = "login_page",
    style = "width: 500px; max-width: 100%; margin: 0 auto; padding-top: 100px;",
    wellPanel(
      h2("MAMOTH Login", class = "text-center", style = "padding-top: 0;"),
      p("Please enter the password to access the application.", class = "text-center"),
      passwordInput("password", "Password:"),
      div(
        class = "text-center",
        actionButton("login_button", "Log in", class = "btn-primary")
      ),
      textOutput("login_message")
    )
  ),
  
  # This div will contain the main app UI after successful login
  # It is hidden by default
  shinyjs::hidden(
    div(
      id = "main_app",
      # Your original navbarPage UI goes here. It is now wrapped in this hidden div.
      navbarPage(
        title = div(img(src = mammoth_icon_url, height = "25px", style = "margin-top: -5px; padding-right: 10px;"), "MAMOTH"),
        theme = shinytheme("sandstone"),
        id = "main_nav",
        
        tags$head(
          tags$style(HTML("
                    body { color: #545252; }
                    .navbar-default {
                      background-color: #F8F5F0 !important; 
                      /* This rule removes the navbar's own border */
                      border: none !important; 
                    }
                    .navbar-default .navbar-brand, .navbar-default .navbar-nav > li > a {
                      color: #545252 !important; font-size: 16px !important;
                    }
                    .navbar-default .navbar-nav > .active > a,
                    .navbar-default .navbar-nav > .active > a:hover,
                    .navbar-default .navbar-nav > .active > a:focus {
                      color: #545252 !important; background-color: #FFFFFF !important;
                    }
                    #nav_gene_viewer, #nav_multi_omics, #nav_mofa, #nav_about {
                      background-color: #F8F5F0 !important; border-color: #E7E7E7 !important; color: #545252 !important;
                    }
                    #nav_gene_viewer:hover, #nav_multi_omics:hover, #nav_mofa:hover, #nav_about:hover {
                      background-color: #E7E7E7 !important;
                    }

                    /* --- THIS RULE REMOVES THE TOP BORDER OF THE PAGE CONTENT --- */
                    .navbar + .container-fluid {
                        border-top: none !important;
                    }

                    /* --- CSS RULE FOR bslib SIDEBARS (MOFA TABS) --- */
                    .bslib-sidebar-layout > aside {
                        background-color: #F8F5F0 !important;
                    }
                    
                    /* --- CSS RULES FOR ABOUT PAGE FONT SIZE --- */
                    #about-section p {
                        font-size: 16px !important; 
                    }
                    #about-section h4 {
                        font-size: 22px !important; 
                    }

                  "))
        ),
        
        tabPanel("Home",
                 fluidPage(
                   fluidRow(
                     column(12, align = "center", style = "padding-top: 50px;",
                            tags$h1(
                              img(src = mammoth_icon_url, height="68px", style = "vertical-align: -4px; margin-right: 20px;"),
                              "MAMOTH", 
                              style = "font-size: 72px; font-weight: bold;"
                            ),
                            tags$h3("MAGEL2 Multi-Omics Targeting Hubs")
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(10, offset = 1, align = "center",
                            p("Welcome to the MAMOTH interactive data portal. This application provides tools to explore multi-omics datasets (transcriptome, proteome, and ubiquitome) related to MAGEL2 research. Use the navigation panels below or the tabs at the top to access the different viewers.", style = "font-size: 16px;")
                     )
                   ),
                   br(),
                   fluidRow(
                     column(3, div(class = "panel", style = "padding: 15px; text-align: center;", actionButton("nav_gene_viewer", label = tagList(h4(icon("dna"), "Gene Expression Viewer")), width = "100%"), p("Visualize the expression of individual genes."))),
                     column(3, div(class = "panel", style = "padding: 15px; text-align: center;", actionButton("nav_multi_omics", label = tagList(h4(icon("sitemap"), "Multi-Omics Browser")), width = "100%"), p("Explore differential expression and enrichment results."))),
                     column(3, div(class = "panel", style = "padding: 15px; text-align: center;", actionButton("nav_mofa", label = tagList(h4(icon("project-diagram"), "MOFA Browser")), width = "100%"), p("Deconvolve multi-omic variability."))),
                     column(3, div(class = "panel", style = "padding: 15px; text-align: center;", actionButton("nav_about", label = tagList(h4(icon("info-circle"), "About")), width = "100%"), p("Learn more about the project and data.")))
                   )
                 )
        ),
        
        tabPanel("Gene expression viewer",
                 fluidPage(
                   sidebarLayout(
                     sidebarPanel(
                       width = 3,
                       h4("Select Gene"),
                       selectizeInput("goi", "Gene Symbol:", choices = NULL, options = list(placeholder = 'Type or select gene...', maxOptions = 10000)),
                       hr(),
                       h4("Display Options"),
                       checkboxInput("show_trans", "Show Transcriptome", value = TRUE),
                       checkboxInput("show_prot", "Show Proteome", value = TRUE),
                       checkboxInput("show_ubi", "Show Ubiquitome", value = TRUE),
                       hr(),
                       helpText("Select a gene to view its expression. Use checkboxes to toggle data types. Significant stars refer to the adjusted p-value compared to WT")
                     ),
                     mainPanel(
                       width = 9,
                       conditionalPanel("input.show_trans == true", fluidRow(column(6, h5("Transcriptome - 8330", align="center"), plotOutput("plotRna8330", height = "350px")), column(6, h5("Transcriptome - MGH", align="center"), plotOutput("plotRnaMGH", height = "350px")))),
                       conditionalPanel("input.show_prot == true", fluidRow(column(6, h5("Proteome - 8330", align="center"), plotOutput("plotProt8330", height = "350px")), column(6, h5("Proteome - MGH", align="center"), plotOutput("plotProtMGH", height = "350px")))),
                       conditionalPanel("input.show_ubi == true", fluidRow(column(6, h5("Ubiquitome - 8330", align="center"), plotOutput("plotUbi8330", height = "350px")), column(6, h5("Ubiquitome - MGH", align="center"), plotOutput("plotUbiMGH", height = "350px"))))
                     )
                   )
                 )
        ),
        
        tabPanel("Multi-Omics Browser",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Volcano Plot", 
                            fluidPage(
                              br(),
                              # --- Reverted to classic sidebarLayout ---
                              sidebarLayout(
                                sidebarPanel(
                                  width = 3,
                                  h4("Controls"),
                                  selectInput("volcano_group1", "Select Group 1:", choices = NULL),
                                  selectInput("volcano_group2", "Select Group 2 (Contrast):", choices = NULL),
                                  hr(),
                                  checkboxGroupInput("volcano_omics", "Select Omics Layers:", choices = c("Transcriptome", "Proteome", "Ubiquitome"), selected = c("Transcriptome", "Proteome", "Ubiquitome")),
                                  hr(),
                                  helpText("Significant hits (p.adj < 0.05 and |LFC| > 1) are colored red."),
                                  hr(),
                                  downloadButton("download_volcano_results", "Download Volcano Results (.xlsx)")
                                ),
                                mainPanel(
                                  width = 9,
                                  conditionalPanel("input.volcano_omics.includes('Transcriptome')", fluidRow(column(6, plotOutput("volcano_rna_8330")), column(6, plotOutput("volcano_rna_MGH")))),
                                  conditionalPanel("input.volcano_omics.includes('Proteome')", fluidRow(column(6, plotOutput("volcano_prot_8330")), column(6, plotOutput("volcano_prot_MGH")))),
                                  conditionalPanel("input.volcano_omics.includes('Ubiquitome')", fluidRow(column(6, plotOutput("volcano_ubi_8330")), column(6, plotOutput("volcano_ubi_MGH"))))
                                )
                              )
                            )
                   ),
                   tabPanel("Gene Ontology Enrichment",
                            fluidPage(
                              br(),
                              # --- Reverted to classic sidebarLayout ---
                              sidebarLayout(
                                sidebarPanel(
                                  width = 3,
                                  h4("Analysis Controls"),
                                  selectInput("go_omics_type", "Select Omics Layer:", choices = c("Transcriptome", "Proteome", "Ubiquitome")),
                                  selectInput("go_ontology", "Select GO Ontology:", choices = c("Biological Process" = "BP", "Cellular Component" = "CC", "Molecular Function" = "MF")),
                                  selectInput("go_group1", "Select Group 1:", choices = NULL),
                                  selectInput("go_group2", "Select Group 2 (Contrast):", choices = NULL),
                                  sliderInput("go_n_terms", "Number of Top Terms to Display:", min = 3, max = 20, value = 10, step = 1),
                                  actionButton("go_run_analysis", "Run Analysis", icon = icon("rocket"), class = "btn-primary"),
                                  hr(),
                                  helpText("This tool finds genes significantly changed in both backgrounds and performs GO enrichment."),
                                  hr(),
                                  downloadButton("download_go_results", "Download GO Results (.xlsx)")
                                ),
                                mainPanel(
                                  width = 9,
                                  plotOutput("go_plot", height = "700px")
                                )
                              )
                            )
                   )
                 )
        ),
        
        tabPanel("MOFA Browser",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Overview",
                            fluidPage(
                              titlePanel("Overview of trained Multi-Omics Factor Analysis (MOFA) Model"),
                              h3("Factor Overview"),
                              plotOutput("overview_factor_plot", height="800px"),
                              hr(),
                              h3("Total Variance Explained"),
                              plotOutput("plot_variance_group_total")
                            )
                   ),
                   tabPanel("Factor Exploration",
                            fluidPage(
                              br(), # Added for spacing
                              sidebarLayout(
                                sidebarPanel(
                                  h4("Factor Plot Options"),
                                  checkboxGroupInput("factors_to_plot", "Select Factors to Plot:", choices = 1:N_FACTORS, selected = 1:5, inline = TRUE),
                                  selectInput("factors_color_by", "Color Samples By:", choices = c("Background", "Mutation"), selected = "Background"),
                                  helpText("Select at least two factors to generate plots.")
                                ),
                                mainPanel(
                                  h3("Factor Scatter Plot Matrix"),
                                  plotOutput("factor_plot_matrix", height = "600px"),
                                  hr(),
                                  h3("Factor Correlations"),
                                  plotOutput("factor_correlation_heatmap")
                                )
                              )
                            )
                   ),
                   tabPanel("Feature Weights",
                            fluidPage(
                              br(), # Added for spacing
                              sidebarLayout(
                                sidebarPanel(h4("Weight Plot Options"), selectInput("weights_view", "Select Omics View:", choices = view_choices), selectInput("weights_factor", "Select Factor:", choices = 1:N_FACTORS, selected = 1), sliderInput("weights_nfeatures", "Number of Top Features:", min = 5, max = 250, value = 10, step = 1)),
                                mainPanel(h3("Top Feature Weights"), uiOutput("plot_top_weights_selected_ui"))
                              )
                            )
                   ),
                   tabPanel("Data Heatmaps",
                            fluidPage(
                              br(), # Added for spacing
                              # --- Reverted to classic sidebarLayout ---
                              sidebarLayout(
                                sidebarPanel(
                                  width = 3,
                                  h4("Heatmap Options"),
                                  selectInput("heatmap_view", "Select Omics View:", choices = view_choices),
                                  selectInput("heatmap_factor", "Select Factor:", choices = 1:N_FACTORS, selected = 1),
                                  sliderInput("heatmap_nfeatures", "Number of Top Features:", min = 10, max = 100, value = 25, step = 5),
                                  checkboxInput("heatmap_cluster_rows", "Cluster Rows (Features)", value = TRUE),
                                  checkboxInput("heatmap_cluster_cols", "Cluster Columns (Samples)", value = TRUE),
                                  checkboxInput("heatmap_show_rownames", "Show Feature Names", value = TRUE),
                                  checkboxInput("heatmap_show_colnames", "Show Sample Names", value = TRUE),
                                  selectInput("heatmap_scale", "Scale Data:", choices = c("Row" = "row", "Column" = "column", "None" = "none"), selected = "row")
                                ),
                                mainPanel(
                                  width = 9,
                                  h3("Heatmap of Top Features for Selected Factor"), 
                                  uiOutput("plot_data_heatmap_selected_ui")
                                )
                              )
                            )
                   ),
                   tabPanel("GSEA Enrichment",
                            fluidPage(
                              br(), # Added for spacing
                              # --- Reverted to classic sidebarLayout ---
                              sidebarLayout(
                                sidebarPanel(
                                  width = 3,
                                  h4("GSEA Options (Live Analysis)"),
                                  selectInput("gsea_view", "Select Omics View:", choices = gsea_view_choices),
                                  selectInput("gsea_factor", "Select Factor:", choices = 1:N_FACTORS, selected = 1),
                                  selectInput("gsea_database", "Select Gene Set Database:",
                                              choices = c("Reactome" = "REAC", "Gene Ontology (BP)" = "GO:BP", "Gene Ontology (MF)" = "GO:MF", "Gene Ontology (CC)" = "GO:CC", "Human Phenotype Ontology" = "HP", "MicroRNAs (miRTarBase)" = "MIRNA", "Transcription Factors (TRANSFAC)" = "TF"),
                                              selected = "REAC"),
                                  sliderInput("gsea_n_pathways", "Number of Top Pathways:", min = 5, max = 30, value = 15, step = 1),
                                  hr(),
                                  actionButton("run_gsea", "Run Analysis", icon = icon("rocket")),
                                  helpText("Click 'Run Analysis' to fetch results from g:Profiler. This may take a moment."),
                                  hr(),
                                  downloadButton("download_gsea_results", "Download GSEA Results (.xlsx)")
                                ),
                                mainPanel(
                                  width = 9,
                                  h3("Live Gene Set Enrichment Analysis (GSEA) using g:Profiler"),
                                  plotOutput("gsea_combined_plot", height = "600px")
                                )
                              )
                            )
                   )
                 )
        ),
        
        tabPanel("About",
                 fluidPage(
                   id = "about-section", 
                   titlePanel("About the MAMOTH Application"),
                   hr(),
                   fluidRow(
                     column(8,
                            h4("General Information"),
                            p("The MAGEL2 Multi-Omics Targeting Hubs (MAMOTH) application is an interactive, open-access web portal developed to empower the research community and facilitate the exploration of this complex dataset. It provides complete and user-friendly access to the entire transcriptome, proteome, and ubiquitinome datasets presented in the study, enabling researchers to independently visualize findings, test novel hypotheses, and accelerate progress in understanding the molecular basis of MAGEL2-related neurodevelopmental disorders."),
                            p("The data originates from a comprehensive multi-omics analysis of CRISPR/Cas9-engineered isogenic human pluripotent stem cell (hiPSC)-derived cortical neurons published in in Buecking et al., 2025 preprint."),
                            
                            h4("Background"),
                            p("Variants in the gene MAGEL2 are associated with the severe neurodevelopmental disorders Prader-Willi Syndrome (PWS) and the more severe Schaaf-Yang Syndrome (SYS), yet the underlying molecular pathophysiology in the human brain remains poorly understood. Despite known links to processes like endosomal trafficking and protein ubiquitination, the specific role of MAGEL2 in human cortical neurons during critical developmental stages has been unclear."),
                            p("This project addresses this gap by presenting the first integrative multi-omics analysis for MAGEL2-related disorders."),
                            
                            h4("Feedback and Contact"),
                            p("This application is under active development. We welcome any feedback, bug reports, or suggestions for new features. For inquiries, please contact the Laugsch Lab at Heidelberg University."),
                            p(HTML("<strong>Contact:</strong> <a href='mailto:jannis.buecking@gmail.com'>jannis.buecking@gmail.com</a>")),
                            p(HTML("<strong>Lab Website:</strong> <a href='https://www.klinikum.uni-heidelberg.de/humangenetik/forschung/ag-laugsch' target='_blank'>Laugsch Lab</a>")),
                            
                            h4("Funding"),
                            p("This project was funded by the Foundation for Prader-Willi Research (FPWR)."),
                            
                            h4("Credits"),
                            p(HTML("Mammoth icon created by <a href='https://www.flaticon.com/free-icons/mammoth' title='mammoth icons'>Freepik - Flaticon</a>."))
                     )
                   )
                 )
        )
      )
    )
  )
)



# --- 4. SERVER LOGIC ---
server <- function(input, output, session) {
  
  USER <- reactiveValues(authenticated = FALSE)
  
  observeEvent(input$login_button, {
    if (input$password == "Mammuthus_primigenius") {
      USER$authenticated <- TRUE
      shinyjs::hide("login_page")
      shinyjs::show("main_app")
    } else {
      output$login_message <- renderText("Invalid password.")
    }
  })
  
  
  # --- Home Page Navigation ---
  observeEvent(input$nav_gene_viewer, { updateNavbarPage(session, "main_nav", selected = "Gene expression viewer") })
  observeEvent(input$nav_multi_omics, { updateNavbarPage(session, "main_nav", selected = "Multi-Omics Browser") })
  observeEvent(input$nav_mofa, { updateNavbarPage(session, "main_nav", selected = "MOFA Browser") })
  observeEvent(input$nav_about, { updateNavbarPage(session, "main_nav", selected = "About") })
  
  # --- Server Logic for Gene Expression Viewer ---
  observe({ 
    updateSelectizeInput(session, "goi", choices = available_gois, selected = "RHOA", server = TRUE) 
  })
  
  plot_colors <- tryCatch({ gg_color_hue(2) }, error = function(e) { c("grey80", "grey70") })
  my_color_8330 <- plot_colors[1]; my_color_MGH <- plot_colors[2]
  render_abundance_plot <- function(data_obj, dep_results_obj, color, label, star_sz = 6, star_v = -0.5, plot_subtitle = "") {
    renderPlot({
      validate(need(data_loaded_flag, "Waiting for data..."), need(!is.null(master_condition_order), "Condition order not set."),
               need(!is.null(data_obj), paste(plot_subtitle,"-",label,"data object not loaded.")))
      selected_goi <- input$goi; validate(need(selected_goi, "Please select a GOI."))
      plotAbundance(data=data_obj, GOI=selected_goi, dep_results=dep_results_obj, ref_group="WT",
                    plot_color=color, data_source_label=plot_subtitle, master_order=master_condition_order,
                    star_size=star_sz, star_vjust=star_v)
    })
  }
  output$plotProt8330 <- render_abundance_plot(new_data_imp_8330, dep_8330, my_color_8330, "Prot 8330", plot_subtitle="Proteome")
  output$plotProtMGH <- render_abundance_plot(new_data_imp_MGH, dep_MGH, my_color_MGH, "Prot MGH", plot_subtitle="Proteome", 5, -0.2)
  output$plotUbi8330 <- render_abundance_plot(Ubi_data_imp_8330, Ubi_dep_8330, my_color_8330, "Ubi 8330", plot_subtitle="Ubiquitome")
  output$plotUbiMGH <- render_abundance_plot(Ubi_data_imp_MGH, Ubi_dep_MGH, my_color_MGH, "Ubi MGH", plot_subtitle="Ubiquitome", 5, -0.2)
  output$plotRna8330 <- renderPlot({
    validate(need(data_loaded_flag, "Waiting for data..."), need(rna_data_prepared, "Transcriptome data prep failed."),
             need(!is.null(master_condition_order), "Condition order not set."))
    selected_goi <- input$goi; validate(need(selected_goi, "Please select a GOI."))
    plotAbundance(data=dds_RNA_8330, GOI=selected_goi, dep_results=DEG_list_8330, ref_group="WT",
                  plot_color=my_color_8330, data_source_label="Transcriptome", master_order=master_condition_order, star_size=6, star_vjust=-0.5)
  })
  output$plotRnaMGH <- renderPlot({
    validate(need(data_loaded_flag, "Waiting for data..."), need(rna_data_prepared, "Transcriptome data prep failed."),
             need(!is.null(master_condition_order), "Condition order not set."))
    selected_goi <- input$goi; validate(need(selected_goi, "Please select a GOI."))
    plotAbundance(data=dds_RNA_MGH, GOI=selected_goi, dep_results=DEG_list_MGH, ref_group="WT",
                  plot_color=my_color_MGH, data_source_label="Transcriptome", master_order=master_condition_order, star_size=5, star_vjust=-0.2)
  })
  
  # --- Server Logic for Volcano Plots ---
  observe({
    req(data_loaded_flag, !is.null(available_conditions))
    choices <- available_conditions
    selected_g1 <- if ("DupC" %in% choices) "DupC" else choices[1]
    selected_g2 <- if ("WT" %in% choices) "WT" else if (length(choices) > 1) choices[2] else choices[1]
    updateSelectInput(session, "volcano_group1", choices = choices, selected = selected_g1)
    updateSelectInput(session, "volcano_group2", choices = choices, selected = selected_g2)
  })
  create_volcano_plot <- function(results_df, plot_title) {
    validate(need(is.data.frame(results_df) && nrow(results_df) > 0, "Comparison not available."))
    validate(need(all(c("name", "logFC", "p_adj") %in% names(results_df)), "Results missing required columns."))
    df <- results_df %>% filter(!is.na(p_adj) & !is.na(logFC)) %>% mutate(significant = ifelse(p_adj < 0.05 & abs(logFC) > 1, "Yes", "No"), label = ifelse(rank(-abs(logFC)) <= 10, name, ""))
    
    ggplot(df, aes(x = logFC, y = -log10(p_adj))) +
      geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
      geom_text_repel(aes(label = label), max.overlaps = 15, size = 3.5, box.padding = 0.5, color = "#545252") +
      
      # --- CHANGE IS HERE ---
      # The legend is removed by setting guide = "none"
      scale_color_manual(values = c("Yes" = "firebrick", "No" = "grey50"), guide = "none") +
      
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
      labs(title = plot_title, x = bquote("Log"[2]*" Fold Change"), y = bquote(-log[10]~'(p.adj)')) +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        text = element_text(color = "#545252"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
      )
  }
  volcano_data <- reactive({
    req(input$volcano_group1, input$volcano_group2, data_loaded_flag)
    validate(need(input$volcano_group1 != input$volcano_group2, "Please select two different groups."))
    g1 <- input$volcano_group1; g2 <- input$volcano_group2
    process_se_data <- function(se_object) {
      meta <- as.data.frame(rowData(se_object))
      lfc_col <- paste0(g1, "_vs_", g2, "_diff"); padj_col <- paste0(g1, "_vs_", g2, "_p.adj")
      if (all(c(lfc_col, padj_col) %in% names(meta))) { return(data.frame(name = meta$name, logFC = meta[[lfc_col]], p_adj = meta[[padj_col]])) }
      lfc_col_rev <- paste0(g2, "_vs_", g1, "_diff"); padj_col_rev <- paste0(g2, "_vs_", g1, "_p.adj")
      if (all(c(lfc_col_rev, padj_col_rev) %in% names(meta))) { return(data.frame(name = meta$name, logFC = -1 * meta[[lfc_col_rev]], p_adj = meta[[padj_col_rev]])) }
      return(data.frame())
    }
    prot_8330_df <- process_se_data(dep_8330); prot_MGH_df <- process_se_data(dep_MGH)
    ubi_8330_df <- process_se_data(Ubi_dep_8330); ubi_MGH_df <- process_se_data(Ubi_dep_MGH)
    process_deg_list <- function(deg_list){
      list_name_fwd <- paste0(g1, "_vs_", g2)
      if (list_name_fwd %in% names(deg_list)) { res <- as.data.frame(deg_list[[list_name_fwd]]); return(data.frame(name = rownames(res), logFC = res$log2FoldChange, p_adj = res$padj)) }
      list_name_rev <- paste0(g2, "_vs_", g1)
      if (list_name_rev %in% names(deg_list)) { res <- as.data.frame(deg_list[[list_name_rev]]); return(data.frame(name = rownames(res), logFC = -1 * res$log2FoldChange, p_adj = res$padj)) }
      return(data.frame())
    }
    rna_8330_df <- process_deg_list(DEG_list_8330); rna_MGH_df <- process_deg_list(DEG_list_MGH)
    list(prot_8330=prot_8330_df, prot_MGH=prot_MGH_df, ubi_8330=ubi_8330_df, ubi_MGH=ubi_MGH_df, rna_8330=rna_8330_df, rna_MGH=rna_MGH_df)
  })
  output$volcano_prot_8330 <- renderPlot({ create_volcano_plot(volcano_data()$prot_8330, "Proteome - 8330") })
  output$volcano_prot_MGH  <- renderPlot({ create_volcano_plot(volcano_data()$prot_MGH,  "Proteome - MGH") })
  output$volcano_ubi_8330  <- renderPlot({ create_volcano_plot(volcano_data()$ubi_8330,  "Ubiquitome - 8330") })
  output$volcano_ubi_MGH   <- renderPlot({ create_volcano_plot(volcano_data()$ubi_MGH,   "Ubiquitome - MGH") })
  output$volcano_rna_8330  <- renderPlot({ create_volcano_plot(volcano_data()$rna_8330,  "Transcriptome - 8330") })
  output$volcano_rna_MGH   <- renderPlot({ create_volcano_plot(volcano_data()$rna_MGH,   "Transcriptome - MGH") })
  
  # --- Server Logic for GO Enrichment ---
  observe({
    req(data_loaded_flag); choices <- available_conditions
    updateSelectInput(session, "go_group1", choices = choices, selected = choices[1])
    updateSelectInput(session, "go_group2", choices = choices, selected = choices[2])
  })
  go_results <- eventReactive(input$go_run_analysis, {
    req(input$go_group1 != input$go_group2)
    withProgress(message = 'Running GO Analysis...', value = 0, {
      g1 <- input$go_group1; g2 <- input$go_group2; omics <- input$go_omics_type
      incProgress(0.1, detail = paste("Fetching", tolower(omics), "data"))
      if (omics == "Transcriptome") {
        get_deg_df <- function(deg_list, g1, g2) {
          fwd <- paste0(g1, "_vs_", g2); rev <- paste0(g2, "_vs_", g1)
          if (fwd %in% names(deg_list)) { df <- as.data.frame(deg_list[[fwd]]); df$name <- rownames(df); df$logFC <- df$log2FoldChange; return(df) }
          if (rev %in% names(deg_list)) { df <- as.data.frame(deg_list[[rev]]); df$name <- rownames(df); df$logFC <- -df$log2FoldChange; return(df) }
          return(NULL)
        }
        df_8330 <- get_deg_df(DEG_list_8330, g1, g2); df_mgh <- get_deg_df(DEG_list_MGH, g1, g2)
      } else {
        get_se_df <- function(se_obj, g1, g2) {
          meta <- as.data.frame(rowData(se_obj))
          fwd_lfc <- paste0(g1, "_vs_", g2, "_diff"); fwd_padj <- paste0(g1, "_vs_", g2, "_p.adj")
          rev_lfc <- paste0(g2, "_vs_", g1, "_diff"); rev_padj <- paste0(g2, "_vs_", g1, "_p.adj")
          if (all(c(fwd_lfc, fwd_padj) %in% names(meta))) { return(data.frame(name = meta$name, logFC = meta[[fwd_lfc]], padj = meta[[fwd_padj]])) }
          if (all(c(rev_lfc, rev_padj) %in% names(meta))) { return(data.frame(name = meta$name, logFC = -meta[[rev_lfc]], padj = meta[[rev_padj]])) }
          return(NULL)
        }
        if (omics == "Proteome") { df_8330 <- get_se_df(dep_8330, g1, g2); df_mgh <- get_se_df(dep_MGH, g1, g2)
        } else { df_8330 <- get_se_df(Ubi_dep_8330, g1, g2); df_mgh <- get_se_df(Ubi_dep_MGH, g1, g2) }
      }
      validate(need(!is.null(df_8330) && !is.null(df_mgh), "This comparison is not available."))
      incProgress(0.3, detail = "Finding overlapping genes")
      sig_up_8330 <- df_8330 %>% filter(padj < 0.05, logFC > 0) %>% pull(name)
      sig_down_8330 <- df_8330 %>% filter(padj < 0.05, logFC < 0) %>% pull(name)
      sig_up_mgh <- df_mgh %>% filter(padj < 0.05, logFC > 0) %>% pull(name)
      sig_down_mgh <- df_mgh %>% filter(padj < 0.05, logFC < 0) %>% pull(name)
      overlap_up <- intersect(sig_up_8330, sig_up_mgh)
      overlap_down <- intersect(sig_down_8330, sig_down_mgh)
      run_go <- function(gene_list, progress_msg, progress_val) {
        incProgress(progress_val, detail = progress_msg)
        if (length(gene_list) == 0) return(NULL)
        entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
        if (length(entrez_ids) == 0) return(NULL)
        enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = input$go_ontology, pAdjustMethod = "BH", readable = TRUE)
      }
      go_up <- run_go(overlap_up, "Analyzing upregulated genes", 0.2)
      go_down <- run_go(overlap_down, "Analyzing downregulated genes", 0.2)
      incProgress(0.2, detail = "Finalizing")
      list(upregulated = go_up, downregulated = go_down)
    })
  })
  output$go_plot <- renderPlot({
    results <- go_results()
    validate(need(!is.null(results), "Click 'Run Analysis' to generate results."), need(!is.null(results$upregulated) || !is.null(results$downregulated), "No enrichment results could be generated."))
    all_results_df <- bind_rows(as.data.frame(results$upregulated), as.data.frame(results$downregulated))
    global_plot_limit <- if (nrow(all_results_df) > 0) { ceiling(max(-log10(all_results_df$p.adjust), na.rm = TRUE)) } else { 10 }
    create_go_barplot <- function(results_df, direction_label) {
      if (is.null(results_df) || nrow(as.data.frame(results_df)) == 0) { return(ggplot() + labs(title = direction_label, subtitle = "No significant terms found") + theme_void() + theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))) }
      plot_data <- as.data.frame(results_df) %>% mutate(NegLog10Padj = -log10(p.adjust)) %>% slice_max(order_by = NegLog10Padj, n = input$go_n_terms) %>% mutate(Description_wrapped = fct_reorder(str_wrap(Description, 50), NegLog10Padj))
      color_palette <- if (direction_label == "Overlapping Upregulation") "Blues" else "Reds"
      ggplot(plot_data, aes(x = Description_wrapped, y = NegLog10Padj, fill = NegLog10Padj)) +
        geom_col() + coord_flip() +
        scale_fill_distiller(palette = color_palette, direction = 1, limits = c(0, global_plot_limit)) +
        scale_y_continuous(limits = c(0, global_plot_limit), expand = c(0, 0.1)) +
        labs(title = direction_label, x = NULL, y = bquote(-log[10]~'(Adjusted P-value)')) +
        theme_bw(base_size = 14) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
    }
    p_up <- create_go_barplot(results$upregulated, "Overlapping Upregulation")
    p_down <- create_go_barplot(results$downregulated, "Overlapping Downregulation")
    p_up + p_down
  })
  
  # --- Download Handler for Volcano Plot Results ---
  output$download_volcano_results <- downloadHandler(
    filename = function() {
      req(input$volcano_group1, input$volcano_group2)
      paste0("Volcano_Results_", input$volcano_group1, "_vs_", input$volcano_group2, "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      data_to_write <- volcano_data()
      
      # Filter the list to only include non-empty data frames
      sheets_list <- list()
      if(nrow(data_to_write$rna_8330) > 0) sheets_list[["Transcriptome_8330"]] <- data_to_write$rna_8330
      if(nrow(data_to_write$rna_MGH) > 0) sheets_list[["Transcriptome_MGH"]] <- data_to_write$rna_MGH
      if(nrow(data_to_write$prot_8330) > 0) sheets_list[["Proteome_8330"]] <- data_to_write$prot_8330
      if(nrow(data_to_write$prot_MGH) > 0) sheets_list[["Proteome_MGH"]] <- data_to_write$prot_MGH
      if(nrow(data_to_write$ubi_8330) > 0) sheets_list[["Ubiquitome_8330"]] <- data_to_write$ubi_8330
      if(nrow(data_to_write$ubi_MGH) > 0) sheets_list[["Ubiquitome_MGH"]] <- data_to_write$ubi_MGH
      
      validate(need(length(sheets_list) > 0, "No data available to download for the selected comparison."))
      
      openxlsx::write.xlsx(sheets_list, file = file)
    }
  )
  
  # --- Download Handler for Gene Ontology Results ---
  output$download_go_results <- downloadHandler(
    filename = function() {
      req(input$go_group1, input$go_group2, input$go_omics_type, input$go_ontology)
      paste0("GO_Enrichment_", input$go_omics_type, "_", input$go_ontology, "_", 
             input$go_group1, "_vs_", input$go_group2, "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      data_to_write <- go_results()
      validate(need(!is.null(data_to_write), "Please run the GO analysis before downloading."))
      
      sheets_list <- list()
      if (!is.null(data_to_write$upregulated) && nrow(as.data.frame(data_to_write$upregulated)) > 0) {
        sheets_list[["Overlapping_Upregulation"]] <- as.data.frame(data_to_write$upregulated)
      }
      if (!is.null(data_to_write$downregulated) && nrow(as.data.frame(data_to_write$downregulated)) > 0) {
        sheets_list[["Overlapping_Downregulation"]] <- as.data.frame(data_to_write$downregulated)
      }
      
      validate(need(length(sheets_list) > 0, "No enrichment results found to download."))
      
      openxlsx::write.xlsx(sheets_list, file = file)
    }
  )
  
  # --- SERVER LOGIC FOR MOFA BROWSER ---
  base_gg_theme <- theme_minimal(base_size = 18) + theme(text = element_text(colour = font_color_dark), axis.text = element_text(colour = font_color_dark), plot.title = element_text(colour = font_color_dark), legend.text = element_text(colour = font_color_dark))
  
  output$overview_factor_plot <- renderPlot({
    req(mofa_data_loaded)
    
    # Extract the values for factors 1 to 15
    factors_df <- get_factors(MOFA_out, factors = 1:15, as.data.frame = TRUE)
    
    # Extract the pre-processed metadata
    meta_df <- samples_metadata(MOFA_out)
    
    # Combine for plotting
    plot_df <- left_join(factors_df, meta_df, by = "sample")
    
    # Build the single combined plot
    ggplot(plot_df, aes(x = factor, y = value, color = Mutation, shape = Background)) +
      geom_vline(xintercept = 1:14 + 0.5, linetype = "dashed", color = "grey70") +
      
      # 1. Increased point size from 3 to 4
      geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 8, alpha = 0.8) +
      
      # 2. Changed shape values
      #    16 = circle, 17 = triangle pointing up. 
      #    Note: There is no default "triangle pointing down" symbol. 
      #    If you need a downward triangle, you would use shape = 6 (hollow) or shape = 25 (filled with border). 
      #    I will use 17 (up) and 16 (circle) as they are the standard solid equivalents.
      scale_shape_manual(name = "Background", values = c("8330" = 17, "MGH" = 16)) +
      
      scale_color_discrete(name = "Genotype") + 
      labs(x = "", y = "Factor Value") +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 18, color = "#545252"),
        axis.text = element_text(color = "#545252"),
        legend.position = "top"
      )
  })
  
  output$plot_variance_group_total <- renderPlot({ 
    req(mofa_data_loaded)
    mofa_tmp <- MOFA_out; mofa_tmp@samples_metadata$group <- mofa_tmp@samples_metadata$Mutation_BG
    tryCatch({
      plot_data <- plot_variance_explained(mofa_tmp, x = "group", plot_total = TRUE)[[2]]$data
      plot_data <- plot_data %>% mutate(Background = str_extract(group, "(8330|MGH)$"), view = factor(view, levels = c("RNA", "protein", "Ubiq")))
      bg_colors <- c("8330" = "grey40", "MGH" = "grey80") 
      
      ggplot(plot_data, aes(x = group, y = R2, fill = Background)) +
        geom_bar(stat = "identity", color = NA) + 
        facet_wrap(~view, labeller = labeller(view = view_labeller)) +
        base_gg_theme + 
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 16), 
          axis.text.y = element_text(hjust = 1, size = 16), 
          panel.grid = element_blank(), 
          legend.position = "none",
          # Added text color to the facet titles
          strip.text = element_text(size = 16, face = "bold", color = "#545252") 
        ) +
        labs(x = NULL, y = "Total variance explained (%)") +
        scale_fill_manual(name = "Background", values = bg_colors)
    }, error = function(e) { ggplot() + labs(title = "Error generating total variance plot", subtitle = e$message) + theme_void() })
  })
  
  output$factor_plot_matrix <- renderPlot({ 
    req(mofa_data_loaded, input$factors_to_plot, input$factors_color_by)
    factors_to_plot <- as.numeric(input$factors_to_plot)
    color_by_col <- input$factors_color_by
    
    validate(need(length(factors_to_plot) >= 2, "Please select at least two factors."))
    
    factors_df <- get_factors(MOFA_out, factors = factors_to_plot, as.data.frame = TRUE) %>%
      tidyr::pivot_wider(id_cols = "sample", names_from = "factor", values_from = "value")
    
    metadata_df <- MOFA_out@samples_metadata %>%
      tibble::rownames_to_column("sample_rn") %>%
      dplyr::select(sample, !!sym(color_by_col))
    
    plot_df <- dplyr::left_join(factors_df, metadata_df, by = "sample")
    
    p <- GGally::ggpairs(
      plot_df,
      columns = paste0("Factor", factors_to_plot),
      mapping = aes(color = .data[[color_by_col]]),
      
      # --- CHANGE IS HERE ---
      # Kept the default 'cor' function which respects color, 
      # but set the text 'size' to 3 to prevent overlap.
      upper = list(continuous = wrap("cor", size = 3)), 
      
      lower = list(continuous = wrap("points", size = 4, alpha = 0.6)),
      diag = list(continuous = custom_density)
    )
    
    # This logic correctly applies your grey colors only when "Background" is selected
    if (color_by_col == "Background") {
      p <- p + scale_color_manual(values = c("8330" = "grey40", "MGH" = "grey80")) +
        scale_fill_manual(values = c("8330" = "grey40", "MGH" = "grey80"))
    }
    
    p <- p + theme_bw(base_size = 14) +
      theme(
        text = element_text(colour = font_color_dark),
        strip.text = element_text(size = 14, colour = font_color_dark),
        axis.title = element_text(size = 16, colour = font_color_dark),
        axis.text = element_text(size = 12, colour = font_color_dark),
        legend.title = element_text(size = 14, colour = font_color_dark),
        legend.text = element_text(size = 12, colour = font_color_dark),
        panel.grid = element_blank()
      )
    
    print(p)
  })
  
  output$factor_correlation_heatmap <- renderPlot({ 
    req(mofa_data_loaded, input$factors_to_plot)
    selected_factors <- as.numeric(input$factors_to_plot)
    validate(need(length(selected_factors) >= 2, "Please select at least two factors."))
    
    # --- Get Correlation Matrix ---
    all_factors_list <- get_factors(MOFA_out, factors = "all")
    all_factors_matrix <- do.call(rbind, all_factors_list)
    full_cor_matrix <- cor(all_factors_matrix)
    cor_matrix_subset <- full_cor_matrix[selected_factors, selected_factors, drop = FALSE]
    
    # --- Melt matrix for ggplot2 ---
    plot_data <- as.data.frame(cor_matrix_subset) %>%
      # Use a unique temporary name for the new column
      rownames_to_column("Temp_Factor_Row") %>% 
      pivot_longer(
        cols = -Temp_Factor_Row, # Pivot everything except our new row name column
        names_to = "Factor2", 
        values_to = "Correlation"
      ) %>%
      # Rename the temporary column to its final name
      rename(Factor1 = Temp_Factor_Row) %>%
      mutate(
        Factor1 = factor(Factor1, levels = rownames(cor_matrix_subset)),
        Factor2 = factor(Factor2, levels = colnames(cor_matrix_subset))
      )
    
    # --- Create ggplot2 Heatmap ---
    ggplot(plot_data, aes(x = Factor1, y = Factor2, fill = Correlation)) +
      geom_tile(color = "white") + 
      scale_fill_gradientn(
        colors = heatmap_colors, 
        limits = c(-1, 1), 
        name = "Correlation"
      ) +
      geom_text(
        aes(label = round(Correlation, 2)), 
        color = "#545252", 
        size = 5
      ) +
      coord_fixed() +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 18) +
      theme(
        text = element_text(color = "#545252"),
        axis.text = element_text(color = "#545252"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(color = "#545252"),
        legend.text = element_text(color = "#545252"),
        panel.grid = element_blank() 
      )
  })
  
  output$plot_top_weights_selected_ui <- renderUI({ req(input$weights_nfeatures); plotOutput("plot_top_weights_selected", height = paste0(350 + (input$weights_nfeatures * 20), "px")) })
  output$plot_top_weights_selected <- renderPlot({
    req(mofa_data_loaded, input$weights_view, input$weights_factor, input$weights_nfeatures)
    weights_df <- get_weights(MOFA_out, as.data.frame = TRUE, views = input$weights_view, factors = as.numeric(input$weights_factor))
    plot_data <- weights_df %>% slice_max(abs(value), n = input$weights_nfeatures) %>% mutate(feature_formatted = format_feature_names(feature), feature_ordered = fct_reorder(feature_formatted, value))
    max_abs_weight <- max(abs(plot_data$value), na.rm = TRUE)
    ggplot(plot_data, aes(y = feature_ordered, x = value, fill = value)) +
      geom_bar(stat = "identity") + scale_fill_gradientn(colors = heatmap_colors, limits = c(-max_abs_weight, max_abs_weight)) +
      labs(title = paste("Top Features for Factor", input$weights_factor), x = "Weight", y = NULL, fill = "Weight") + base_gg_theme + theme(panel.grid.major.y = element_blank())
  })
  
  output$plot_data_heatmap_selected_ui <- renderUI({ req(input$heatmap_nfeatures); plotOutput("plot_data_heatmap_selected", height = paste0(400 + (input$heatmap_nfeatures * 15), "px")) })
  output$plot_data_heatmap_selected <- renderPlot({
    req(mofa_data_loaded, input$heatmap_view, input$heatmap_factor, input$heatmap_nfeatures)
    weights_matrix <- get_weights(MOFA_out, views = input$heatmap_view)[[1]]; factor_weights <- weights_matrix[, as.numeric(input$heatmap_factor), drop = FALSE]
    top_features_sorted <- factor_weights[order(abs(factor_weights[,1]), decreasing = TRUE), , drop = FALSE]
    features_to_plot <- head(rownames(top_features_sorted), n = input$heatmap_nfeatures)
    data_matrix_list <- MOFA_out@data[[input$heatmap_view]]; full_data_matrix <- do.call(cbind, data_matrix_list)
    data_matrix_subset <- full_data_matrix[features_to_plot, , drop = FALSE]
    rownames(data_matrix_subset) <- format_feature_names(rownames(data_matrix_subset))
    display_view_name <- view_labeller[input$heatmap_view]; title_text <- paste("Top", input$heatmap_nfeatures, display_view_name, "Features for Factor", input$heatmap_factor)
    pheatmap::pheatmap(data_matrix_subset, main = title_text, color = heatmap_colors, border_color = NA, scale = "row", cluster_rows = input$heatmap_cluster_rows, cluster_cols = input$heatmap_cluster_cols, show_rownames = input$heatmap_show_rownames, show_colnames = input$heatmap_show_colnames, fontsize = 12, fontsize_row = 10, fontsize_col = 10, fontfamily = "sans")
  })
  
  gsea_results_reactive <- eventReactive(input$run_gsea, {
    req(mofa_data_loaded)
    TARGET_FACTOR <- as.numeric(input$gsea_factor); SOURCE_GENESET <- input$gsea_database; SELECTED_VIEW <- input$gsea_view
    withProgress(message = 'Running GSEA...', value = 0, {
      incProgress(0.1, detail = "Fetching weights...")
      all_weights <- get_weights(MOFA_out, factors = TARGET_FACTOR, as.data.frame = TRUE)
      run_gsea_internal <- function(weights_df, suffix_to_remove) {
        features <- unique(format_feature_names(gsub(suffix_to_remove, "", weights_df$feature))); if (length(features) < 3) return(NULL)
        gost_result <- gost(query = features, organism = "hsapiens", sources = SOURCE_GENESET, user_threshold = 0.05, correction_method = "fdr")
        if (!is.null(gost_result)) { return(gost_result$result %>% group_by(term_name) %>% slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>% ungroup()) } else { return(NULL) }
      }
      incProgress(0.3, detail = paste("Processing view:", SELECTED_VIEW))
      view_suffix_map <- c("RNA" = "_RNA$", "protein" = "_protein$", "Ubiq" = "_ubiq$"); view_suffix <- view_suffix_map[SELECTED_VIEW]
      weights_pos <- all_weights %>% filter(view == SELECTED_VIEW & value > 0); incProgress(0.2, detail = "Querying g:Profiler (positive weights)..."); results_pos <- if (nrow(weights_pos) > 0) run_gsea_internal(weights_pos, view_suffix) else NULL
      weights_neg <- all_weights %>% filter(view == SELECTED_VIEW & value < 0); incProgress(0.2, detail = "Querying g:Profiler (negative weights)..."); results_neg <- if (nrow(weights_neg) > 0) run_gsea_internal(weights_neg, view_suffix) else NULL
      incProgress(0.2, detail = "Finalizing...")
      list(positive = results_pos, negative = results_neg)
    })
  })
  
  output$gsea_combined_plot <- renderPlot({
    req(mofa_data_loaded)
    results <- gsea_results_reactive()
    TOP_N_TERMS_PLOT <- input$gsea_n_pathways; TARGET_FACTOR <- as.numeric(input$gsea_factor); SOURCE_GENESET <- input$gsea_database; view_display_name <- names(gsea_view_choices)[gsea_view_choices == input$gsea_view]
    validate(need(!is.null(results), "Click 'Run Analysis' to generate GSEA results."))
    all_results_df <- bind_rows(results$positive, results$negative); global_plot_limit <- if (nrow(all_results_df) > 0) ceiling(max(-log10(all_results_df$p_value), na.rm=TRUE)) else 10
    
    create_gsea_plot_internal <- function(results_df, weight_sign) {
      if (is.null(results_df) || nrow(results_df) == 0) { return(ggplot() + labs(title = paste(weight_sign, "Weights Enrichment"), subtitle = "No significant terms found") + theme_void() + theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))) }
      plot_data <- results_df %>% mutate(NegLog10FDR = -log10(p_value)) %>% slice_max(order_by = NegLog10FDR, n = TOP_N_TERMS_PLOT) %>% mutate(term_name_ordered = fct_reorder(str_wrap(term_name, 50), NegLog10FDR))
      color_palette <- if (weight_sign == "Positive") "Blues" else "Reds"
      
      p <- ggplot(plot_data, aes(x = term_name_ordered, y = NegLog10FDR, fill = NegLog10FDR)) +
        geom_bar(stat = "identity") + coord_flip() + scale_fill_distiller(palette = color_palette, direction = 1, limits = c(0, global_plot_limit)) + scale_y_continuous(limits = c(0, global_plot_limit), expand = c(0, 0.1)) +
        labs(
          title = paste(weight_sign, "Weights Enrichment"),
          subtitle = paste("Factor", TARGET_FACTOR, "|", names(which(c("Reactome" = "REAC", "GO (BP)" = "GO:BP", "GO (MF)" = "GO:MF", "GO (CC)" = "GO:CC", "HPO" = "HP", "miRNA" = "MIRNA", "TF" = "TF") == SOURCE_GENESET))),
          x = NULL,
          y = expression(-log10(p.adj)) # *** LABEL UPDATED HERE ***
        ) +
        base_gg_theme + 
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none" # *** LEGEND REMOVED HERE ***
        )
      p
    }
    
    p_pos <- create_gsea_plot_internal(results$positive, "Positive"); p_neg <- create_gsea_plot_internal(results$negative, "Negative")
    title <- ggdraw() + draw_label(paste(view_display_name, "View - GSEA Results"), fontface = 'bold', x = 0.5, hjust = 0.5, size = 16, colour = font_color_dark)
    combined_plots <- plot_grid(p_pos, p_neg, ncol = 2, align = 'v')
    plot_grid(title, combined_plots, ncol = 1, rel_heights = c(0.1, 1))
  })
  
  # --- Download Handler for GSEA Results ---
  output$download_gsea_results <- downloadHandler(
    filename = function() {
      req(input$gsea_view, input$gsea_factor, input$gsea_database)
      paste0("GSEA_Results_", input$gsea_view, "_Factor", input$gsea_factor, "_", 
             input$gsea_database, "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      data_to_write <- gsea_results_reactive()
      validate(
        need(!is.null(data_to_write), "Please run the GSEA analysis before downloading."),
        need(!is.null(data_to_write$positive) || !is.null(data_to_write$negative), "No GSEA results found to download.")
      )
      
      sheets_list <- list()
      if (!is.null(data_to_write$positive) && nrow(data_to_write$positive) > 0) {
        sheets_list[["Positive_Weights_Enrichment"]] <- data_to_write$positive
      }
      if (!is.null(data_to_write$negative) && nrow(data_to_write$negative) > 0) {
        sheets_list[["Negative_Weights_Enrichment"]] <- data_to_write$negative
      }
      
      validate(need(length(sheets_list) > 0, "The analysis ran, but no significant terms were found to download."))
      
      openxlsx::write.xlsx(sheets_list, file = file)
    }
  )
  
  
}


# --- 5. Run the Application ---
shinyApp(ui = ui, server = server)
