# CAD_app.R - CAD Translational Analytics Dashboard

library(DT)
library(shiny)
library(tidyverse)
library(ggpubr)
library(viridis)
library(enrichplot)
library(shinythemes)

load("precalculated_enrichment.RData")

df <- read.table("dataset.txt", header = TRUE)

df$cad.class <- as.factor(df$cad.class)
df$gender <- as.factor(df$gender)
df$diabetes <- as.factor(df$diabetes)
df$hypertension <- as.factor(df$hypertension)
df$hyperlipid <- as.factor(df$hyperlipid)

# Extract gene columns
clinical_cols <- c("cad.class", "gender", "bmi", "diabetes", "hyperlipid", "hypertension", "cxcl5.rank")
gene_cols = setdiff(names(df), clinical_cols)

# Key Metrics
total_patients <- nrow(df)
male_count <- sum(df$gender == "M")
female_count <- sum(df$gender == "F")
male_pct <- round((male_count / total_patients) * 100, 1)
female_pct <- round((female_count / total_patients) * 100, 1)

# Gene Expression for heatmap
gene_profile_matrix <- df %>%
  group_by(gender, cad.class) %>%
  summarise(across(all_of(gene_cols), mean, na.rm = TRUE), .groups = 'drop')

top_var_genes <- names(head(sort(apply(df[, gene_cols], 2, var), decreasing = TRUE), 50))

df_melted_genes <- gene_profile_matrix %>%
  pivot_longer(cols = all_of(top_var_genes), names_to = "Gene", values_to = "Mean_Expression")

# CXCL5 statistics
cor_cad <- cor.test(as.numeric(df$cad.class), df$cxcl5.rank, method = "spearman")
cor_bmi <- cor.test(df$bmi, df$cxcl5.rank, method = "spearman")
wilcox_gender <- wilcox.test(cxcl5.rank ~ gender, data = df)


# User Interface (UI) design
ui <- fluidPage(
  theme = shinytheme("cosmo") %>% tryCatch(error = function(e) NULL),
  
  titlePanel(
    div(style = "padding: 10px 0px;",
        h2("Coronary Artery Disease (CAD) Translational Analytics", style = "font-weight: bold; color: #2c3e50;"),
        h5("An Integrative Dashboard Linking Clinical Covariates, Biomarkers, and Functional Genomics", style = "color: #7f8c8d;")
    )
  ),
  
  # Top row
  fluidRow(
    column(3, wellPanel(style = "background-color: #f8f9fa; border-left: 5px solid #3498db; padding: 10px 15px;",
                        h4("Total Cohort Size", style = "margin:0; color:#7f8c8d; font-size:12px; text-transform: uppercase;"),
                        h2(total_patients, style = "margin:5px 0 0 0; font-weight:bold; color:#2c3e50;"))),
    column(3, wellPanel(style = "background-color: #f8f9fa; border-left: 5px solid #2ecc71; padding: 10px 15px;",
                        h4("Gender Demographics", style = "margin:0; color:#7f8c8d; font-size:12px; text-transform: uppercase;"),
                        h3(sprintf("M: %s%% | F: %s%%", male_pct, female_pct), style = "margin:7px 0 0 0; font-weight:bold; color:#2c3e50; font-size:20px;"))),
    column(3, wellPanel(style = "background-color: #f8f9fa; border-left: 5px solid #e67e22; padding: 10px 15px;",
                        h4("CXCL5 vs CAD Severity", style = "margin:0; color:#7f8c8d; font-size:12px; text-transform: uppercase;"),
                        h3(sprintf("Rho = %s (p = %s)", round(cor_cad$estimate, 2), format.pval(cor_cad$p.value, digits = 2)), 
                           style = "margin:7px 0 0 0; font-weight:bold; color:#2c3e50; font-size:18px;"))),
    column(3, wellPanel(style = "background-color: #f8f9fa; border-left: 5px solid #9b59b6; padding: 10px 15px;",
                        h4("Enriched GO Terms", style = "margin:0; color:#7f8c8d; font-size:12px; text-transform: uppercase;"),
                        h2(nrow(ego_df), style = "margin:5px 0 0 0; font-weight:bold; color:#2c3e50;")))
  ),
  
  # Navigation panel split into tabs
  tabsetPanel(
    type = "tabs",
    
    tabPanel("Patient Cohort & BMI", 
             fluidRow(style = "padding-top: 25px; padding-bottom: 25px;",
                      column(7, 
                             div(style = "padding-right: 10px;",
                                 plotOutput("bmiPlot", height = "500px")
                             )
                      ),
                      column(5, 
                             wellPanel(style = "background-color: #ffffff; border: 1px solid #e3e6f0; box-shadow: 0 0.15rem 1.75rem 0 rgba(58, 59, 69, 0.05); height: 500px; display: flex; flex-direction: column; justify-content: space-between;",
                                       div(
                                         h4("Interesting Cohort Insight", style = "font-weight: bold; color: #2c3e50; margin-top: 0;"),
                                         p("This view analyzes body mass index (BMI) dynamics stratified across the 5 clinical stages of coronary artery disease severity (Classes 0-4) and separated by biological sex.", style = "color: #7f8c8d; font-size: 13px; margin-bottom: 0;")
                                       ),
                                       hr(style = "margin-top: 10px; margin-bottom: 10px;"),
                                       div(style = "flex-grow: 1;",
                                           plotOutput("diabetesBmiPlot", height = "310px")
                                       )
                             )
                      )
             )
    ),
    
    tabPanel("Transcriptomic Profiling", 
             fluidRow(style = "padding-top: 20px; padding-bottom: 10px;",
                      column(12, 
                             wellPanel(style = "background-color: #f8f9fa; border-left: 5px solid #2c3e50; margin-bottom: 15px;",
                                       h4("Expression Topography of Top 50 Highly Variable Structural and Functional Genes", style = "font-weight: bold; color: #2c3e50; margin-top: 0;"),
                                       p("Visualizing transcriptomic markers across severity stages and genders. Notice distinct genes like ", 
                                         tags$b("NPC2"), " (implicated in intracellular cholesterol trafficking) maintaining high baseline expression while displaying significant variance across subgroups.", style = "margin-bottom: 0;")
                             )
                      )
             ),
             fluidRow(
               column(12, plotOutput("heatmapPlot", height = "850px"))
             )
    ),
    
    tabPanel("Functional Genomics (GO ORA)", 
             fluidRow(style = "padding-top: 20px;",
                      column(6, h4("Top Enriched Biological Processes", style="font-weight:bold;"), plotOutput("goDotplot", height = "450px")),
                      column(6, h4("Gene-Concept Network Linkages", style="font-weight:bold;"), plotOutput("goCnetplot", height = "450px"))
             ),
             fluidRow(style = "padding-top: 15px;",
                      column(12, h4("Enriched Pathways Data Table", style="font-weight:bold;"), dataTableOutput("goTable"))
             )
    ),
    
    tabPanel("Comorbidity Risk Analysis", 
             sidebarLayout(
               sidebarPanel(style = "margin-top: 20px;",
                            selectInput("comorb_var", "Select Clinical Comorbidity:",
                                        choices = c("Hyperlipidemia" = "hyperlipid", "Hypertension" = "hypertension")),
                            hr(),
                            h4("Statistical Contingency Summary:"),
                            verbatimTextOutput("contingencyText")
               ),
               mainPanel(style = "margin-top: 20px;",
                         plotOutput("comorbPlot", height = "400px")
               )
             )
    ),
    
    tabPanel("Biomarker Focus: CXCL5", 
             fluidRow(style = "padding-top: 20px;",
                      column(7, plotOutput("cxcl5Comprehensive", height = "400px")),
                      column(5, plotOutput("cxcl5BmiCorr", height = "400px"))
             ),
             fluidRow(
               column(12, wellPanel(
                 h4("CXCL5 Biomarker Discussion", style = "font-weight:bold; color:#2c3e50;"),
                 p("This analysis demonstrates a statistically robust, highly significant negative correlation between ", 
                   tags$code("CXCL5.rank"), " and clinical CAD progression. While it might be suggested that inflammation markers rise during disease progression, ", 
                   tags$b("CXCL5 acts as a known atheroprotective agent and endothelial survival factor."), 
                   " The observed decline signifies exhaustion of vascular compensatory angiogenesis and plaque stabilization systems in late-stage (Class 3-4) multi-vessel coronary disease. Conversely, the marker shows zero significant correlation with BMI, confirming its role as a specific vascular stress indicator rather than an obesity surrogate.")
               ))
             )
    )
  )
)

# Server logic (Plots & Axes - Sentence Case)
server <- function(input, output, session) {
  
  # Tab1.1 BMI distribution plot
  output$bmiPlot <- renderPlot({
    ggplot(df, aes(x = cad.class, y = bmi, fill = gender)) +
      geom_violin(alpha = 0.65, position = position_dodge(0.6), color = "#d1d3e2") +
      geom_boxplot(width = 0.12, position = position_dodge(0.6), color = "#2e3440", outlier.shape = NA, fatten = 2.5) +
      scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.7, labels = c("Female (F)", "Male (M)")) +
      labs(title = "BMI profiles stratified by CAD severity and gender",
           x = "Clinical CAD severity class (0: Control → 4: Severe)", 
           y = "Body mass index (BMI)", 
           fill = "Gender group") +
      theme_minimal(base_size = 13) + 
      theme(
        plot.title = element_text(face = "bold", size = 15, color = "#2c3e50", margin = margin(b = 12)),
        legend.position = "top",
        legend.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "#ccd1d9"),
        axis.line.y = element_line(color = "#ccd1d9")
      )
  })
  
  # Tab1.2 Diabetes boxplot
  output$diabetesBmiPlot <- renderPlot({
    ggplot(df, aes(x = diabetes, y = bmi, fill = diabetes)) +
      geom_boxplot(alpha = 0.8, outlier.shape = 16, width = 0.5, color = "#2e3440") +
      facet_wrap(~gender, labeller = labeller(gender = c("F" = "Female cohort", "M" = "Male cohort"))) +
      scale_fill_manual(values = c("#440154FF", "#21908CFF"), labels = c("No", "Yes")) +
      labs(title = "Correlation between diabetes status and BMI", 
           x = "Diabetes diagnosis status", 
           y = "BMI index", 
           fill = "Diabetes") +
      theme_minimal(base_size = 11) + 
      theme(
        plot.title = element_text(face = "bold", size = 11, color = "#2c3e50", margin = margin(b = 8)),
        legend.position = "none",
        strip.background = element_rect(fill = "#f8f9fa", color = "#e3e6f0", size = 1),
        strip.text = element_text(face = "bold", color = "#2c3e50"),
        panel.grid.minor = element_blank()
      )
  })
  
  # Tab2. High variable genes heatmap
  output$heatmapPlot <- renderPlot({
    ggplot(df_melted_genes, aes(x = cad.class, y = Gene, fill = Mean_Expression)) +
      geom_tile(color = "white", size = 0.1) +
      facet_wrap(~gender) +
      scale_fill_viridis_c(option = "B") +
      labs(x = "CAD severity classification", y = "Annotated gene symbols", fill = "Mean log2 expression") +
      theme_minimal(base_size = 14) +
      theme(
        panel.spacing = unit(3, "lines"), 
        axis.text.y = element_text(size = 10, face = "bold", font = "mono"),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold")
      )
  })
  
  # Tab3. Functional enrichment
  output$goDotplot <- renderPlot({
    req(nrow(ego_df) > 0)
    dotplot(ego, showCategory = 8, title = "GO: Biological processes enrichment") + theme_classic(base_size = 12)
  })
  
  output$goCnetplot <- renderPlot({
    req(nrow(ego_df) > 0)
    cnetplot(ego, layout = "nicely")
  })
  
  output$goTable <- renderDataTable({
    req(nrow(ego_df) > 0)
    ego_df[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust", "geneID")]
  }, options = list(pageLength = 5, scrollX = TRUE))
  
  # Tab4. Comorbidities
  reactive_table <- reactive({
    sub_df <- df %>% 
      filter(cad.class %in% c("0", "3", "4")) %>%
      mutate(Group = ifelse(cad.class == "0", "Class 0 (Control)", "Class 3-4 (Severe CAD)"))
    table(sub_df$Group, sub_df[[input$comorb_var]])
  })
  
  output$contingencyText <- renderPrint({
    tbl <- reactive_table()
    cat("--- Contingency matrix ---\n")
    print(tbl)
    cat("\n--- Inferential testing ---\n")
    if (any(tbl < 5)) {
      cat("Executing Exact Test due to small local cell count (<5):\n")
      print(fisher.test(tbl))
    } else {
      cat("Executing Pearson's Chi-Square test:\n")
      print(chisq.test(tbl))
    }
  })
  
  output$comorbPlot <- renderPlot({
    sub_df <- df %>% 
      filter(cad.class %in% c("0", "3", "4")) %>%
      mutate(Group = ifelse(cad.class == "0", "Class 0 (Control)", "Class 3-4 (Severe CAD)"))
    
    title_label <- ifelse(input$comorb_var == "hyperlipid", "Hyperlipidemia", "Hypertension")
    
    ggplot(sub_df, aes(x = Group, fill = .data[[input$comorb_var]])) +
      geom_bar(position = "fill", width = 0.6) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.6, labels = c("Absence (No)", "Presence (Yes)")) +
      labs(title = paste(title_label, "proportional ratios across CAD cohorts"),
           x = "Phenotypic cohort group", y = "Relative percentage ratio", fill = title_label) +
      theme_classic(base_size = 13)
  })
  
  # Tab5. CXCL5 evaluation
  output$cxcl5Comprehensive <- renderPlot({
    pA <- ggplot(df, aes(x = cad.class, y = cxcl5.rank, fill = cad.class)) +
      geom_boxplot(alpha = 0.8, outlier.shape = 16) +
      scale_fill_viridis_d(option = "E") +
      labs(subtitle = "a) Trend across CAD classes", x = "CAD class", y = "CXCL5 rank score") +
      theme_classic(base_size = 11) + theme(legend.position = "none")
    
    pB <- ggplot(df, aes(x = gender, y = cxcl5.rank, fill = gender)) +
      geom_boxplot(alpha = 0.8, width = 0.5) +
      scale_fill_manual(values = c("#35b719FF", "#4419FF")) +
      labs(subtitle = "b) Variance across genders", x = "Gender", y = "CXCL5 rank score") +
      theme_classic(base_size = 11) + theme(legend.position = "none")
    
    p_comp <- ggarrange(pA, pB, ncol = 2, nrow = 1)
    annotate_figure(p_comp, top = text_grob("CXCL5 kinetics evaluation", face = "bold", size = 14))
  })
  
  output$cxcl5BmiCorr <- renderPlot({
    ggplot(df, aes(x = bmi, y = cxcl5.rank)) +
      geom_point(aes(color = factor(cad.class)), alpha = 0.7, size = 2.5) +
      geom_smooth(method = "lm", color = "darkred", fill = "red", alpha = 0.15, linetype = "dashed") +
      scale_color_viridis_d(option = "D", name = "CAD class") +
      labs(
        title = "Independence check: BMI vs CXCL5 kinetics",
        subtitle = paste0("Spearman's Rho = ", round(cor_bmi$estimate, 3), 
                          " (p = ", format.pval(cor_bmi$p.value, digits = 3), ")"),
        x = "Body mass index (BMI)", y = "CXCL5 rank score"
      ) +
      theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  })
}

# Execution
shinyApp(ui = ui, server = server)