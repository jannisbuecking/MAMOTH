# --- 0. Load Shiny and Other Libraries ---
# Ensure all these packages are installed first
library(shiny)
library(shinyjs)
library(openxlsx)
library(shinythemes)
library(DESeq2)
library(bslib)                
library(msigdbr)              # Add this line
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)                
library(stringr)
library(shinycssloaders)
library(Hmisc)
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
font_color_dark <- "#545252"
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
    return(ggplot() + theme_void() + 
             annotate("text", x=0.5, y=0.5, label="Gene not found\nin this dataset", size = 6, color = "#545252") + 
             ggtitle(GOI) + labs(subtitle=NULL) +
             theme(plot.title = element_text(hjust = 0.5, color = "#545252"))) # Add theme here too
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
  
  p <- ggplot(plot_df, aes(x = condition, y = expression)) +
    geom_jitter(color = plot_color, width = 0.2, size = 3, alpha = 0.7) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                 geom = "errorbar", width = 0.2, color = plot_color) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = plot_color) +
    scale_x_discrete(drop = FALSE) +
    
    # --- LABS UPDATED ---
    labs(title = GOI, subtitle = NULL, y = y_axis_label, x = NULL) + # Subtitle removed
    
    theme_bw(base_size = 16) +
    
    # --- THEME UPDATED ---
    theme(
      text = element_text(color = "#545252"), 
      axis.text = element_text(color = "#545252"),
      axis.ticks = element_line(color = "#545252"), 
      axis.title.y = element_text(size = 12, color = "#545252"),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "#545252"), # Color updated
      axis.line.y = element_line(color = "#545252"), # Color updated
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", color = "#545252"), # Color added
      plot.subtitle = element_text(hjust = 0.5, color = "#545252"), # Color added
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 20, r = 5, b = 5, l = 5, "pt")
    )
  # --- END OF THEME UPDATE ---
  
  if (nrow(annotation_df) > 0) {
    p <- p + geom_text(data = annotation_df, aes(x = condition, y = y_pos, label = label), 
                       inherit.aes = FALSE, size = star_size, vjust = star_vjust, color = "#545252") # Color added
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
        actionButton("login_button", "Log in")
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
                    
                    /* --- CSS RULES FOR ABOUT PAGE FONT SIZE --- */
                    #about-section p {
                        font-size: 16px !important; 
                    }
                    #about-section h4 {
                        font-size: 22px !important; 
                    }

                    /* --- CORRECTED CSS RULE TO STYLE H3 TITLES --- */
                    /* This targets h3 elements inside the main content area */
                    .col-sm-9 h3 {
                        font-size: 26px !important;
                        font-weight: bold !important;
                        color: #545252 !important;
                        text-align: center;
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
                            p("Welcome to the MAMOTH interactive data portal. This application provides tools to explore multi-omics datasets (transcriptome, proteome, and ubiquitome) from human induced pluripotent stem cell (hiPSC)-derived cortical neurons to investigate MAGEL2-associated diseases. Use the navigation panels below or the tabs at the top to access the different viewers.", style = "font-size: 16px;")
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
                       h4(""), 
                       selectizeInput("goi", "Enter Gene Symbol", choices = NULL, options = list(placeholder = 'Type or select gene...', maxOptions = 10000)),
                       hr(),
                       
                       tags$label("Select Omics View", class = "control-label"),
                       checkboxInput("show_trans", "Transcriptome", value = TRUE),
                       checkboxInput("show_prot", "Proteome", value = TRUE),
                       checkboxInput("show_ubi", "Ubiquitome", value = TRUE),
                       
                       hr(),
                       helpText("Select a gene to view its expression. Use checkboxes to toggle data types. Significant stars refer to the adjusted p-value compared to WT")
                     ),
                     
                     # --- UI UPDATED HERE ---
                     mainPanel(
                       width = 9,
                       # A single spinner now wraps a single UI output
                       withSpinner(uiOutput("gene_plots_ui"), type = 6, color = "#545252", size = 2)
                     )
                     # --- END OF UPDATE ---
                   )
                 )
        ),
        
        tabPanel("Multi-Omics Browser",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Volcano Plot", 
                            fluidPage(
                              br(),
                              sidebarLayout(
                                sidebarPanel(
                                  width = 3,
                                  selectInput("volcano_group1", "Select Genotype 1:", choices = NULL),
                                  selectInput("volcano_group2", "Select Genotype 2 (Contrast):", choices = NULL),
                                  hr(),
                                  checkboxGroupInput("volcano_omics", "Select Omics Layers:", choices = c("Transcriptome", "Proteome", "Ubiquitome"), selected = c("Transcriptome", "Proteome", "Ubiquitome")),
                                  hr(),
                                  helpText("Significant hits (p.adj < 0.05 and |LFC| > 1) are colored red."),
                                  hr(),
                                  downloadButton("download_volcano_results", "Download Volcano Results (.xlsx)")
                                ),
                                
                                # --- UI UPDATED HERE ---
                                mainPanel(
                                  width = 9,
                                  # A single spinner now wraps a single UI output
                                  withSpinner(uiOutput("volcano_plots_ui"), type = 6, color = "#545252", size = 2)
                                )
                              )
                            )
                   ),
                   tabPanel("Enrichment Analysis",
                            fluidPage(
                              br(),
                              sidebarLayout(
                                sidebarPanel(
                                  width = 3,
                                  h4(""),
                                  selectInput("go_omics_type", "Select Omics Layer:", 
                                              choices = c("Transcriptome", "Proteome", "Ubiquitome")),
                                  
                                  # --- UPDATED DROPDOWN ---
                                  selectInput("go_database", "Select Database:",
                                              choices = c("GO: Combined (All)" = "ALL",
                                                          "GO: Biological Process" = "BP",
                                                          "GO: Cellular Component" = "CC",
                                                          "GO: Molecular Function" = "MF",
                                                          "Reactome" = "REAC",
                                                          "HPO (Human Phenotype)" = "HP",
                                                          "TF (MSigDB C3)" = "TF",
                                                          "miRNA (miRTarBase)" = "MIRNA",
                                                          "Cell Types (MSigDB C8)" = "C8"),
                                              selected = "BP"),
                                  
                                  selectInput("go_background_select", "Select Genetic Background:",
                                              choices = c("Overlapping Genes (Both)" = "both",
                                                          "8330 Background" = "8330",
                                                          "MGH Background" = "MGH")),
                                  selectInput("go_group1", "Select Genotype 1:", choices = NULL),
                                  selectInput("go_group2", "Select Genotype 2 (Contrast):", choices = NULL),
                                  sliderInput("go_n_terms", "Number of Top Terms to Display:", 
                                              min = 10, max = 100, value = 10, step = 10),
                                  actionButton("go_run_analysis", "Run Analysis", icon = icon("rocket")),
                                  hr(),
                                  helpText("Select an omics layer, database, and comparison to run enrichment analysis. Running the enrichment may take a moment."),
                                  hr(),
                                  downloadButton("download_go_results", "Download Enrichment Results (.xlsx)")
                                ),
                                mainPanel(
                                  width = 9,
                                  fluidRow(
                                    column(6, withSpinner(uiOutput("go_plot_down_ui"), type = 6, color = "#545252", size = 2)),
                                    column(6, withSpinner(uiOutput("go_plot_up_ui"), type = 6, color = "#545252", size = 2))
                                  )
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
                              br(), # Added for spacing
                              h3("Overview of trained Multi-Omics Factor Analysis (MOFA) Model", align = "center"),
                              
                              withSpinner(plotOutput("overview_factor_plot", height="800px"), type = 6, color = "#545252", size = 2),
                              hr(),
                              h3("Total Variance Explained", align = "center"), # Also styled this title
                              withSpinner(plotOutput("plot_variance_group_total"), type = 6, color = "#545252", size = 2)
                            )
                   ),
                   tabPanel("Factor Exploration",
                            fluidPage(
                              br(), # Added for spacing
                              sidebarLayout(
                                sidebarPanel(
                                  checkboxGroupInput("factors_to_plot", "Select Factors to Plot:", choices = 1:N_FACTORS, selected = 1:5, inline = TRUE),
                                  
                                  # --- CHOICES RENAMED HERE ---
                                  selectInput("factors_color_by", "Color Samples By:",
                                              choices = c("Genetic Background" = "Background", 
                                                          "Genotype" = "Mutation"),
                                              selected = "Background"),
                                  
                                  helpText("Select at least two factors to generate plots.")
                                ),
                                mainPanel(
                                  h3("Factor Scatter Plot Matrix"),
                                  withSpinner(plotOutput("factor_plot_matrix", height = "600px"), type = 6, color = "#545252", size = 2),
                                  hr(),
                                  h3("Factor Correlations"),
                                  withSpinner(plotOutput("factor_correlation_heatmap"), type = 6, color = "#545252", size = 2)
                                )
                              )
                            )
                   ),
                   tabPanel("Feature Weights",
                            fluidPage(
                              br(), # Added for spacing
                              sidebarLayout(
                                sidebarPanel(selectInput("weights_view", "Select Omics View:", choices = gsea_view_choices), selectInput("weights_factor", "Select Factor:", choices = 1:N_FACTORS, selected = 1), sliderInput("weights_nfeatures", "Number of Top Features:", min = 5, max = 250, value = 10, step = 1)),
                                mainPanel(withSpinner(uiOutput("plot_top_weights_selected_ui"), type = 6, color = "#545252", size = 2))
                              )
                            )
                   ),
                   tabPanel("Data Heatmaps",
                            fluidPage(
                              br(), # Added for spacing
                              sidebarLayout(
                                sidebarPanel(
                                  width = 3,
                                  selectInput("heatmap_view", "Select Omics View:", choices = gsea_view_choices),
                                  selectInput("heatmap_factor", "Select Factor:", choices = 1:N_FACTORS, selected = 1),
                                  sliderInput("heatmap_nfeatures", "Number of Top Features:", min = 10, max = 100, value = 25, step = 5),
                                  
                                  # --- CLUSTERING CHECKBOXES REMOVED ---
                                  
                                  checkboxInput("heatmap_show_rownames", "Show Feature Names", value = TRUE),
                                  checkboxInput("heatmap_show_colnames", "Show Sample Names", value = TRUE)
                                ),
                                mainPanel(
                                  width = 9,
                                  withSpinner(uiOutput("plot_data_heatmap_selected_ui"), type = 6, color = "#545252", size = 2)
                                )
                              )
                            )
                   ),
                   tabPanel("GSEA Enrichment",
                            fluidPage(
                              br(), # Added for spacing
                              sidebarLayout(
                                sidebarPanel(
                                  width = 3,
                                  selectInput("gsea_view", "Select Omics View:", choices = gsea_view_choices),
                                  selectInput("gsea_factor", "Select Factor:", choices = 1:N_FACTORS, selected = 1),
                                  selectInput("gsea_database", "Select Gene Set Database:",
                                              choices = c("Reactome" = "REAC", "Gene Ontology (BP)" = "GO:BP", "Gene Ontology (MF)" = "GO:MF", "Gene Ontology (CC)" = "GO:CC", "Human Phenotype Ontology" = "HP", "MicroRNAs (miRTarBase)" = "MIRNA", "Transcription Factors (TRANSFAC)" = "TF"),
                                              selected = "REAC"),
                                  sliderInput("gsea_n_pathways", "Number of Top Pathways:", min = 10, max = 100, value = 10, step = 10),
                                  hr(),
                                  actionButton("run_gsea", "Run Analysis", icon = icon("rocket")),
                                  hr(),
                                  helpText("Click 'Run Analysis' to fetch results from g:Profiler. This may take a moment."),
                                  hr(),
                                  downloadButton("download_gsea_results", "Download GSEA Results (.xlsx)")
                                ),
                                mainPanel(
                                  width = 9,

                                  # --- UI UPDATED HERE ---
                                  withSpinner(uiOutput("gsea_combined_plot_ui"), type = 6, color = "#545252", size = 2)
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
  
  
  output$gene_plots_ui <- renderUI({
    tagList(
      conditionalPanel("input.show_trans == true", 
                       fluidRow(
                         column(6, h4("Transcriptome - 8330", align="center"), plotOutput("plotRna8330", height = "350px")), 
                         column(6, h4("Transcriptome - MGH", align="center"), plotOutput("plotRnaMGH", height = "350px"))
                       )),
      conditionalPanel("input.show_prot == true", 
                       fluidRow(
                         column(6, h4("Proteome - 8330", align="center"), plotOutput("plotProt8330", height = "350px")), 
                         column(6, h4("Proteome - MGH", align="center"), plotOutput("plotProtMGH", height = "350px"))
                       )),
      conditionalPanel("input.show_ubi == true", 
                       fluidRow(
                         column(6, h4("Ubiquitome - 8330", align="center"), plotOutput("plotUbi8330", height = "350px")), 
                         column(6, h4("Ubiquitome - MGH", align="center"), plotOutput("plotUbiMGH", height = "350px"))
                       ))
    )
  })
  
  # --- ADD THIS NEW BLOCK TO YOUR SERVER ---
  
  output$volcano_plots_ui <- renderUI({
    # This tagList contains all the conditional panels that were
    # previously in the UI. The spinner will now watch this
    # entire block.
    tagList(
      conditionalPanel("input.volcano_omics.includes('Transcriptome')", 
                       fluidRow(
                         column(6, plotOutput("volcano_rna_8330")), 
                         column(6, plotOutput("volcano_rna_MGH"))
                       )),
      conditionalPanel("input.volcano_omics.includes('Proteome')", 
                       fluidRow(
                         column(6, plotOutput("volcano_prot_8330")), 
                         column(6, plotOutput("volcano_prot_MGH"))
                       )),
      conditionalPanel("input.volcano_omics.includes('Ubiquitome')", 
                       fluidRow(
                         column(6, plotOutput("volcano_ubi_8330")), 
                         column(6, plotOutput("volcano_ubi_MGH"))
                       ))
    )
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
  # --- Server Logic for GO Enrichment ---
  observe({
    req(data_loaded_flag, !is.null(available_conditions))
    choices <- available_conditions
    
    # --- DEFAULTS UPDATED HERE ---
    selected_g1 <- if ("DupC" %in% choices) "DupC" else choices[1]
    selected_g2 <- if ("WT" %in% choices) "WT" else if (length(choices) > 1) choices[2] else choices[1]
    
    updateSelectInput(session, "go_group1", choices = choices, selected = selected_g1)
    updateSelectInput(session, "go_group2", choices = choices, selected = selected_g2)
  })
  go_results <- eventReactive(input$go_run_analysis, {
    req(input$go_group1 != input$go_group2)
    
    withProgress(message = 'Running Enrichment Analysis...', value = 0, {
      
      g1 <- input$go_group1; g2 <- input$go_group2; omics <- input$go_omics_type
      bg_select <- input$go_background_select; db_select <- input$go_database
      
      # --- 1. Get Gene Lists ---
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
      
      incProgress(0.2, detail = "Finding significant genes")
      sig_up_8330 <- df_8330 %>% filter(padj < 0.05, logFC > 0) %>% pull(name)
      sig_down_8330 <- df_8330 %>% filter(padj < 0.05, logFC < 0) %>% pull(name)
      sig_up_mgh <- df_mgh %>% filter(padj < 0.05, logFC > 0) %>% pull(name)
      sig_down_mgh <- df_mgh %>% filter(padj < 0.05, logFC < 0) %>% pull(name)
      
      if (bg_select == "both") {
        final_up_list <- intersect(sig_up_8330, sig_up_mgh)
        final_down_list <- intersect(sig_down_8330, sig_down_mgh)
      } else if (bg_select == "8330") {
        final_up_list <- sig_up_8330
        final_down_list <- sig_down_8330
      } else if (bg_select == "MGH") {
        final_up_list <- sig_up_mgh
        final_down_list <- sig_down_mgh
      }
      
      # --- 2. Define Analysis Functions ---
      
      # UPDATED: This function now takes the ontology (BP, CC, MF, or ALL) as an argument
      run_go_analysis <- function(gene_list, ont_selection) {
        if (length(gene_list) == 0) return(NULL)
        entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
        if (length(entrez_ids) == 0) return(NULL)
        
        universe_ids <- NULL
        if (omics == "Proteome" || omics == "Ubiquitome" || (omics == "Transcriptome" && bg_select == "MGH")) {
          if (omics == "Transcriptome") all_genes <- rownames(DEG_list_MGH[[1]])
          if (omics == "Proteome") all_genes <- rowData(dep_MGH)$name
          if (omics == "Ubiquitome") all_genes <- rowData(Ubi_dep_MGH)$name
          universe_ids <- bitr(all_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
        } 
        
        enrichGO(gene = entrez_ids, universe = universe_ids, OrgDb = org.Hs.eg.db, 
                 ont = ont_selection, # Use the selection from the UI
                 pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
      }
      
      run_gost_analysis <- function(gene_list, database) {
        if (length(gene_list) < 3) return(NULL)
        gost(query = gene_list, organism = "hsapiens", sources = database, user_threshold = 0.05, correction_method = "fdr")$result
      }
      
      run_enricher_analysis <- function(gene_list, category) {
        if (length(gene_list) == 0) return(NULL)
        msigdb_cat <- if(category == "TF") "C3" else "C8"
        db <- msigdbr(species = "Homo sapiens", category = msigdb_cat)
        if (category == "TF") {
          db <- db %>% dplyr::filter(str_starts(gs_subcat, "TFT"))
        }
        term2gene <- db %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
        enricher(gene = gene_list, TERM2GENE = term2gene, pAdjustMethod = "BH", pvalueCutoff = 0.05)
      }
      
      # --- 3. Execute Selected Analysis ---
      incProgress(0.4, detail = "Running enrichment...")
      
      # UPDATED: This block now handles the separate GO categories
      if (db_select %in% c("ALL", "BP", "CC", "MF")) {
        go_up <- run_go_analysis(final_up_list, db_select)
        go_down <- run_go_analysis(final_down_list, db_select)
      } else if (db_select %in% c("HP", "REAC", "MIRNA")) {
        go_up <- run_gost_analysis(final_up_list, db_select)
        go_down <- run_gost_analysis(final_down_list, db_select)
      } else if (db_select == "TF") {
        go_up <- run_enricher_analysis(final_up_list, "TF")
        go_down <- run_enricher_analysis(final_down_list, "TF")
      } else if (db_select == "C8") {
        go_up <- run_enricher_analysis(final_up_list, "C8")
        go_down <- run_enricher_analysis(final_down_list, "C8")
      }
      
      incProgress(0.3, detail = "Finalizing")
      list(upregulated = go_up, downregulated = go_down)
    })
  })
  # --- This entire block REPLACES the old "output$go_plot" ---
  
  # 1. New plot output for UPREGULATED terms
  # --- This entire block REPLACES the old "output$go_plot_up" ---
  
  # 1. New plot output for UPREGULATED terms (now in Red)
  # --- This entire block REPLACES the old "output$go_plot_combined" ---
  
  # --- Helper function for GO plots (to avoid repeating code) ---
  # --- Helper function for GO plots (to avoid repeating code) ---
  # --- Helper function for GO plots (AESTHETICS UPDATED) ---
  # --- Helper function for GO/Enrichment plots ---
  # --- Helper function for GO/Enrichment plots ---
  # --- Helper function for GO/Enrichment plots ---
  # --- Helper function for GO/Enrichment plots ---
  # --- Helper function for GO/Enrichment plots ---
  create_go_barplot <- function(results_df, direction_label, global_plot_limit) {
    
    color_palette <- if (direction_label == "Overlapping Downregulation") "Blues" else "Reds"
    
    if (is(results_df, "enrichResult")) {
      results_df <- as.data.frame(results_df@result)
    } else {
      results_df <- as.data.frame(results_df) # Assumes gost result
    }
    
    if (is.null(results_df) || nrow(results_df) == 0) { 
      return(ggplot() + 
               labs(title = direction_label, subtitle = "No significant terms found") + 
               theme_void() + 
               theme(plot.title = element_text(hjust=0.5, face="bold", color="#545252"), 
                     plot.subtitle = element_text(hjust=0.5, color="#545252"))) 
    }
    
    if ("p.adjust" %in% colnames(results_df)) {
      plot_data <- results_df %>% 
        mutate(NegLog10Padj = -log10(p.adjust)) %>% 
        slice_max(order_by = NegLog10Padj, n = input$go_n_terms) %>% 
        mutate(Term_Wrapped = fct_reorder(str_wrap(gsub("_", " ", Description), 30), NegLog10Padj))
    } else {
      plot_data <- results_df %>% 
        mutate(NegLog10Padj = -log10(p_value)) %>% 
        slice_max(order_by = NegLog10Padj, n = input$go_n_terms) %>% 
        mutate(Term_Wrapped = fct_reorder(str_wrap(gsub("_", " ", term_name), 30), NegLog10Padj))
    }
    
    ggplot(plot_data, aes(x = Term_Wrapped, y = NegLog10Padj, fill = NegLog10Padj)) +
      geom_col(width = 0.8) + # Keep the bar width relative but full
      coord_flip() +
      scale_fill_distiller(palette = color_palette, direction = 1, limits = c(0, global_plot_limit)) +
      scale_y_continuous(limits = c(0, global_plot_limit), expand = c(0, 0.1)) +
      labs(title = direction_label, x = NULL, y = expression(-log10(p.adj))) +
      theme_minimal(base_size = 18) + 
      theme(
        text = element_text(color = "#545252"),
        # --- TEXT OVERLAP FIX ---
        axis.text.y = element_text(color = "#545252", lineheight = 0.8), # Adjust line height
        axis.text.x = element_text(color = "#545252"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold", color = "#545252"), 
        legend.position = "none"
      )
  }
  
  
  
  # --- Render UI for DOWNREGULATED plot (for dynamic height) ---
  output$go_plot_down_ui <- renderUI({
    req(input$go_n_terms)
    
    # --- NEW HEIGHT CALCULATION ---
    # 150px base + 45px for each term. This gives ample space.
    plot_height <- 150 + (input$go_n_terms * 45)
    
    withSpinner(plotOutput("go_plot_down", height = paste0(plot_height, "px")), type = 6, color = "#545252", size = 2)
  })
  
  # --- Render the DOWNREGULATED (Blue) plot ---
  output$go_plot_down <- renderPlot({
    results <- go_results()
    validate(need(!is.null(results), "Click 'Run Analysis' to generate results."), 
             need(!is.null(results$downregulated), "No overlapping downregulated genes found."))
    
    all_results_df <- bind_rows(as.data.frame(results$upregulated), as.data.frame(results$downregulated))
    global_plot_limit <- if (nrow(all_results_df) > 0) { ceiling(max(-log10(all_results_df$p.adjust), na.rm = TRUE)) } else { 10 }
    
    create_go_barplot(results$downregulated, "Overlapping Downregulation", global_plot_limit)
  })
  
  # --- Render UI for UPREGULATED plot (for dynamic height) ---
  output$go_plot_up_ui <- renderUI({
    req(input$go_n_terms)
    
    # --- NEW HEIGHT CALCULATION ---
    # 150px base + 45px for each term
    plot_height <- 150 + (input$go_n_terms * 45)
    
    plotOutput("go_plot_up", height = paste0(plot_height, "px"))
  })
  
  # --- Render the UPREGULATED (Red) plot ---
  output$go_plot_up <- renderPlot({
    results <- go_results()
    validate(need(!is.null(results), "Click 'Run Analysis' to generate results."), 
             need(!is.null(results$upregulated), "No overlapping upregulated genes found."))
    
    all_results_df <- bind_rows(as.data.frame(results$upregulated), as.data.frame(results$downregulated))
    global_plot_limit <- if (nrow(all_results_df) > 0) { ceiling(max(-log10(all_results_df$p.adjust), na.rm = TRUE)) } else { 10 }
    
    create_go_barplot(results$upregulated, "Overlapping Upregulation", global_plot_limit)
  })
  
  # --- 1. Render the DOWNREGULATED (Blue) plot ---
  # --- 1. Render the DOWNREGULATED (Blue) plot ---
  output$go_plot_down <- renderPlot({
    results <- go_results()
    validate(need(!is.null(results), "Click 'Run Analysis' to generate results."), 
             need(!is.null(results$downregulated), "No significant terms found for downregulated set."))
    
    # Handle different result types for global limit
    df_up <- if(is(results$upregulated, "enrichResult")) as.data.frame(results$upregulated@result) else as.data.frame(results$upregulated)
    df_down <- if(is(results$downregulated, "enrichResult")) as.data.frame(results$downregulated@result) else as.data.frame(results$downregulated)
    
    all_results_df <- bind_rows(df_up, df_down)
    
    # Handle different p-value column names
    pval_col <- if ("p.adjust" %in% colnames(all_results_df)) "p.adjust" else "p_value"
    
    global_plot_limit <- if (nrow(all_results_df) > 0) { 
      ceiling(max(-log10(all_results_df[[pval_col]]), na.rm = TRUE)) 
    } else { 10 }
    
    create_go_barplot(results$downregulated, "Overlapping Downregulation", global_plot_limit)
  })
  
  # --- 2. Render the UPREGULATED (Red) plot ---
  output$go_plot_up <- renderPlot({
    results <- go_results()
    validate(need(!is.null(results), "Click 'Run Analysis' to generate results."), 
             need(!is.null(results$upregulated), "No significant terms found for upregulated set."))
    
    df_up <- if(is(results$upregulated, "enrichResult")) as.data.frame(results$upregulated@result) else as.data.frame(results$upregulated)
    df_down <- if(is(results$downregulated, "enrichResult")) as.data.frame(results$downregulated@result) else as.data.frame(results$downregulated)
    
    all_results_df <- bind_rows(df_up, df_down)
    
    pval_col <- if ("p.adjust" %in% colnames(all_results_df)) "p.adjust" else "p_value"
    
    global_plot_limit <- if (nrow(all_results_df) > 0) { 
      ceiling(max(-log10(all_results_df[[pval_col]]), na.rm = TRUE)) 
    } else { 10 }
    
    create_go_barplot(results$upregulated, "Overlapping Upregulation", global_plot_limit)
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
  # --- Download Handler for Gene Ontology Results ---
  output$download_go_results <- downloadHandler(
    filename = function() {
      req(input$go_group1, input$go_group2, input$go_omics_type, input$go_database, input$go_background_select)
      
      bg_label <- switch(input$go_background_select,
                         "both" = "Overlapping",
                         "8330" = "8330_Only",
                         "MGH" = "MGH_Only")
      
      paste0("Enrichment_", input$go_omics_type, "_", input$go_database, "_", 
             bg_label, "_", input$go_group1, "_vs_", input$go_group2, "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      data_to_write <- go_results()
      validate(need(!is.null(data_to_write), "Please run the enrichment analysis before downloading."))
      
      # Helper to correctly extract data.frame from enrichResult or gost result
      extract_df <- function(obj) {
        if (is(obj, "enrichResult")) {
          return(as.data.frame(obj@result))
        } else if (is.data.frame(obj)) {
          return(obj)
        }
        return(NULL)
      }
      
      df_up <- extract_df(data_to_write$upregulated)
      df_down <- extract_df(data_to_write$downregulated)
      
      sheets_list <- list()
      if (!is.null(df_up) && nrow(df_up) > 0) {
        sheets_list[["Upregulated_Enrichment"]] <- df_up
      }
      if (!is.null(df_down) && nrow(df_down) > 0) {
        sheets_list[["Downregulated_Enrichment"]] <- df_down
      }
      
      validate(need(length(sheets_list) > 0, "No enrichment results found to download."))
      
      openxlsx::write.xlsx(sheets_list, file = file)
    }
  )
  
  # --- SERVER LOGIC FOR MOFA BROWSER ---
  base_gg_theme <- theme_minimal(base_size = 18) + 
    theme(
      text = element_text(colour = font_color_dark), 
      axis.text = element_text(colour = font_color_dark), 
      plot.title = element_text(colour = font_color_dark), 
      legend.text = element_text(colour = font_color_dark)
    )
  
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

    # --- DEFINE COLORS AND WRAPPER FUNCTION HERE ---
    
    # 1. Define the correct color palette first
    my_colors <- if (color_by_col == "Background") {
      c("8330" = "grey40", "MGH" = "grey80")
    } else {
      # Get the default ggplot colors for the 7 mutation types
      gg_color_hue(n_distinct(plot_df[[color_by_col]]))
    }

    # 2. Create a simple wrapper for the correlation function to set the size
    #    ggally_cor is designed to respect the aes(color=...) mapping
    custom_cor <- function(data, mapping, ...) {
      ggally_cor(data = data, mapping = mapping, size = 3, ...)
    }
    
    # --- END NEW LOGIC ---

    p <- GGally::ggpairs(
      plot_df,
      columns = paste0("Factor", factors_to_plot),
      mapping = aes(color = .data[[color_by_col]]),
      
      # Use the new custom function
      upper = list(continuous = custom_cor), 
      
      lower = list(continuous = wrap("points", size = 4, alpha = 0.6)),
      diag = list(continuous = custom_density)
    )

    # 3. Apply the correct color palette to all plot components
    p <- p + scale_color_manual(values = my_colors) +
             scale_fill_manual(values = my_colors)

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
      geom_bar(stat = "identity") + 
      scale_fill_gradientn(colors = heatmap_colors, limits = c(-max_abs_weight, max_abs_weight)) +
      labs(title = paste("Top Features for Factor", input$weights_factor), x = "Weight", y = NULL, fill = "Weight") +
      
      # --- AESTHETICS UPDATED HERE ---
      base_gg_theme + # Apply the corrected base theme
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        # No manual color overrides are needed anymore
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  })
  
  output$plot_data_heatmap_selected_ui <- renderUI({ 
    req(input$heatmap_nfeatures)
    # 250px base + 15px for each feature
    plot_height <- 250 + (input$heatmap_nfeatures * 15) 
    plotOutput("plot_data_heatmap_selected", height = paste0(plot_height, "px"))
  })
  
  output$plot_data_heatmap_selected <- renderPlot({
    req(mofa_data_loaded, input$heatmap_view, input$heatmap_factor, input$heatmap_nfeatures)
    
    # --- 1. Data Preparation (same as before) ---
    weights_matrix <- get_weights(MOFA_out, views = input$heatmap_view)[[1]]
    factor_weights <- weights_matrix[, as.numeric(input$heatmap_factor), drop = FALSE]
    top_features_sorted <- factor_weights[order(abs(factor_weights[,1]), decreasing = TRUE), , drop = FALSE]
    features_to_plot <- head(rownames(top_features_sorted), n = input$heatmap_nfeatures)
    
    data_matrix_list <- MOFA_out@data[[input$heatmap_view]]
    full_data_matrix <- do.call(cbind, data_matrix_list)
    
    features_in_data <- features_to_plot[features_to_plot %in% rownames(full_data_matrix)]
    if(length(features_in_data) == 0) stop("No matching features found between weights and data.")
    data_matrix_subset <- full_data_matrix[features_in_data, , drop = FALSE]
    
    # --- 2. Z-Score Calculation (same as before) ---
    col_data <- as.data.frame(MOFA_out@samples_metadata)
    mgh_cols <- col_data$sample[col_data$Background == "MGH"]
    b8330_cols <- col_data$sample[col_data$Background == "8330"]
    mgh_cols_present <- intersect(mgh_cols, colnames(data_matrix_subset))
    b8330_cols_present <- intersect(b8330_cols, colnames(data_matrix_subset))
    
    mgh_matrix_unscaled <- data_matrix_subset[, mgh_cols_present, drop = FALSE]
    b8330_matrix_unscaled <- data_matrix_subset[, b8330_cols_present, drop = FALSE]
    
    calculate_z_score <- function(mat) {
      mat_scaled <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      variances <- apply(mat, 1, var, na.rm=TRUE)
      has_variance <- (variances > 1e-6) & !is.na(variances)
      if(any(has_variance)) {
        mat_subset <- mat[has_variance, , drop = FALSE] 
        scaled_subset <- t(scale(t(mat_subset), center = TRUE, scale = TRUE))
        mat_scaled[has_variance, ] <- scaled_subset
      }
      return(mat_scaled)
    }
    
    z_mgh <- calculate_z_score(mgh_matrix_unscaled)
    z_8330 <- calculate_z_score(b8330_matrix_unscaled)
    z_score_matrix_combined <- cbind(z_mgh, z_8330)
    z_score_matrix_final <- z_score_matrix_combined[, colnames(data_matrix_subset)] 
    
    z_score_matrix_final[z_score_matrix_final > 2] <- 2
    z_score_matrix_final[z_score_matrix_final < -2] <- -2
    
    # --- 3. Convert to long data frame for ggplot2 ---
    plot_data <- as.data.frame(z_score_matrix_final) %>%
      rownames_to_column("Feature") %>%
      pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Z_Score") %>%
      mutate(
        # --- CHANGE IS HERE ---
        # Apply the name formatting function to the Feature column
        Feature = factor(format_feature_names(Feature), levels = rev(format_feature_names(features_in_data))), 
        Sample = factor(Sample, levels = colnames(z_score_matrix_final)) 
      )
    
    display_view_name <- view_labeller[input$heatmap_view]
    title_text <- paste("Top", input$heatmap_nfeatures, display_view_name, "Features for Factor", input$heatmap_factor)
    
    # --- 4. Create ggplot2 Heatmap ---
    p <- ggplot(plot_data, aes(x = Sample, y = Feature, fill = Z_Score)) +
      geom_tile() +
      scale_fill_gradientn(
        colors = heatmap_colors, 
        limits = c(-2, 2), 
        name = "Z-Score"
      ) +
      labs(title = title_text, x = NULL, y = NULL) +
      base_gg_theme + # Apply the consistent theme
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    # --- 5. Apply Show/Hide Names Logic ---
    if (input$heatmap_show_rownames) {
      p <- p + theme(axis.text.y = element_text(size = 12))
    } else {
      p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    
    if (input$heatmap_show_colnames) {
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
    } else {
      p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    
    print(p)
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
  
  # --- Render UI for GSEA plot (for dynamic height) ---
  output$gsea_combined_plot_ui <- renderUI({
    req(input$gsea_n_pathways)
    # 150px base + 45px for each term
    plot_height <- 150 + (input$gsea_n_pathways * 45)
    plotOutput("gsea_combined_plot", height = paste0(plot_height, "px"))
  })
  
  output$gsea_combined_plot <- renderPlot({
    req(mofa_data_loaded)
    results <- gsea_results_reactive()
    TOP_N_TERMS_PLOT <- input$gsea_n_pathways; TARGET_FACTOR <- as.numeric(input$gsea_factor); SOURCE_GENESET <- input$gsea_database; view_display_name <- names(gsea_view_choices)[gsea_view_choices == input$gsea_view]
    validate(need(!is.null(results), "Click 'Run Analysis' to generate GSEA results."))
    all_results_df <- bind_rows(results$positive, results$negative); global_plot_limit <- if (nrow(all_results_df) > 0) ceiling(max(-log10(all_results_df$p_value), na.rm=TRUE)) else 10
    
    create_gsea_plot_internal <- function(results_df, weight_sign) {
      if (is.null(results_df) || nrow(results_df) == 0) { return(ggplot() + labs(title = paste(weight_sign, "Weights Enrichment"), subtitle = "No significant terms found") + theme_void() + theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))) }
      
      # --- UPDATED PLOT DATA (TEXT WRAPPING) ---
      plot_data <- results_df %>% 
        mutate(NegLog10FDR = -log10(p_value)) %>% 
        slice_max(order_by = NegLog10FDR, n = TOP_N_TERMS_PLOT) %>% 
        mutate(term_name_ordered = fct_reorder(str_wrap(gsub("_", " ", term_name), 30), NegLog10FDR)) # Wrap at 30 chars
      
      color_palette <- if (weight_sign == "Positive") "Blues" else "Reds"
      
      p <- ggplot(plot_data, aes(x = term_name_ordered, y = NegLog10FDR, fill = NegLog10FDR)) +
        
        # --- UPDATED BAR WIDTH ---
        geom_bar(stat = "identity", width = 0.8) + # Set bar width
        
        coord_flip() + scale_fill_distiller(palette = color_palette, direction = 1, limits = c(0, global_plot_limit)) + scale_y_continuous(limits = c(0, global_plot_limit), expand = c(0, 0.1)) +
        labs(
          title = paste(weight_sign, "Weights"),
          subtitle = paste("Factor", TARGET_FACTOR, "|", names(which(c("Reactome" = "REAC", "GO (BP)" = "GO:BP", "GO (MF)" = "GO:MF", "GO (CC)" = "GO:CC", "HPO" = "HP", "miRNA" = "MIRNA", "TF" = "TF") == SOURCE_GENESET))),
          x = NULL,
          y = expression(-log10(p.adj))
        ) +
        base_gg_theme + 
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none",
          
          # --- UPDATED TEXT SPACING ---
          axis.text.y = element_text(lineheight = 0.8) # Adjust line height for wrapped text
        )
      p
    }
    
    p_pos <- create_gsea_plot_internal(results$positive, "Positive"); p_neg <- create_gsea_plot_internal(results$negative, "Negative")
    title <- ggdraw() + draw_label(paste(view_display_name, "View - GSEA Results"), fontface = 'bold', x = 0.5, hjust = 0.5, size = 16, colour = font_color_dark)
    combined_plots <- plot_grid(p_pos, p_neg, ncol = 2, align = 'v')
    plot_grid(combined_plots, ncol = 1, rel_heights = c(0.1, 1))
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
