# =============================================================================
# Shiny app: Volcano explorer (Script_10)
# =============================================================================
#
# Interactive volcano plot for one comparison of the multiple-factor
# differential-expression analysis (Script_6, Method 4 results): 8 age,
# 8 sex and 24 region comparisons. Points are colored by the estimate and
# -log10(FDR) cutoffs set with the sliders (defaults: |estimate| > log2(1.5),
# FDR < 0.05); changing a cutoff changes only the highlighting, not the
# precomputed statistics.
#
# Data: helper volcano plot.RData (made by Script-10_R_Shiny.R):
#   age.vs.young, female.vs.male, region.vs.region - result tables named by
#   their display label. Comparisons are looked up by these names.
#   (Earlier versions matched labels to tables by list position, which
#   mislabelled the age x sex group of all 24 region comparisons.)
#
# Young = 7 weeks; Aged = 7-11 months.
# Run: shiny::runApp("apps/volcano_explorer")
# =============================================================================

library(shiny)
library(ggplot2)
library(DT)
library(dplyr)

load("helper volcano plot.RData")

# Comparison labels, taken from the names of the result tables.
comparisons <- list(
  age = names(age.vs.young),
  sex = names(female.vs.male),
  region = names(region.vs.region)
)

# All result tables in one list, keyed by "<type>_<label>".
all_deg_data <- c(
  setNames(age.vs.young, paste("age", names(age.vs.young), sep = "_")),
  setNames(female.vs.male, paste("sex", names(female.vs.male), sep = "_")),
  setNames(region.vs.region, paste("region", names(region.vs.region), sep = "_"))
)

# UI section
ui <- fluidPage(
  titlePanel(div("Differential Expression Analysis",
                 tags$br(),
                 tags$small("Young = 7 weeks; Aged = 7-11 months"),
                 style = "color: #333; padding: 10px")),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$head(
        tags$style(HTML("
          .control-label {
            font-size: 14px;
            font-weight: bold;
          }
          .irs-min, .irs-max, .irs-single, .irs-from, .irs-to {
            font-size: 12px !important;
          }
          .selectize-dropdown, .selectize-input {
            font-size: 14px;
          }
          .js-irs-0 .irs-bar,
          .js-irs-1 .irs-bar {
            border-top: 1px solid #ff0000 !important;
            border-bottom: 1px solid #ff0000 !important;
            background: #ff0000 !important;
          }
          .js-irs-0 .irs-bar-edge,
          .js-irs-1 .irs-bar-edge {
            border: 1px solid #ff0000 !important;
            background: #ff0000 !important;
          }
          .js-irs-0 .irs-single,
          .js-irs-1 .irs-single,
          .js-irs-0 .irs-from,
          .js-irs-1 .irs-from,
          .js-irs-0 .irs-to,
          .js-irs-1 .irs-to {
            background: #ff0000 !important;
          }
          .js-irs-0 .irs-single:after,
          .js-irs-1 .irs-single:after,
          .js-irs-0 .irs-from:after,
          .js-irs-1 .irs-from:after,
          .js-irs-0 .irs-to:after,
          .js-irs-1 .irs-to:after {
            border-top-color: #ff0000 !important;
          }
          .js-irs-0 .irs-slider,
          .js-irs-1 .irs-slider {
            box-shadow: 1px 1px 3px rgba(255, 0, 0, 0.3) !important;
            border: 2px solid #ff0000 !important;
            background: #fff !important;
          }
          .js-irs-0 .irs-line,
          .js-irs-1 .irs-line {
            background: #f9f9f9 !important;
            border: 1px solid #ff0000 !important;
          }
          .js-irs-0 .irs-grid-text,
          .js-irs-1 .irs-grid-text {
            color: #ff0000 !important;
          }
          .selectize-input:hover {
            border-color: #ff0000;
          }
          .well {
            background-color: #f8f9fa;
            border: 1px solid #ddd;
            border-radius: 5px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
          }
        "))
      ),
      
      downloadButton("downloadPlot", "Download Plot", 
                     style = "margin-bottom: 15px; width: 100%;"),
      
      selectInput("comparison_type",
                  "Select Comparison Type:",
                  choices = c("Age", "Sex", "Region"),
                  selected = "Age"),
      
      uiOutput("comparison_selector"),
      
      sliderInput("fc_cutoff",
                  "Estimate Change Threshold:",
                  min = 0,
                  max = 4,
                  value = log2(1.5),
                  step = 0.1),
      
      sliderInput("p_cutoff",
                  "- Log10 (FDR)",
                  min = 0,
                  max = 20,
                  value = -log10(0.05),
                  step = 0.1),
      
      sliderInput("point_size",
                  "Point Size:",
                  min = 1,
                  max = 5,
                  value = 2.5,
                  step = 0.5)
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("Volcano Plot", 
                 div(style = "position: relative;",
                     plotOutput("volcanoPlot",
                                click = "plot_click",
                                brush = "plot_brush",
                                height = "800px"),
                     verbatimTextOutput("plot_stats"),
                     verbatimTextOutput("click_info")
                 )),
        tabPanel("Protein Table", 
                 div(style = "margin-top: 20px",
                     downloadButton("downloadData", "Download Data"),
                     DTOutput("geneTable")))
      )
    )
  )
)

# Server section
server <- function(input, output, session) {
  
  output$comparison_selector <- renderUI({
    comp_type <- tolower(input$comparison_type)
    selectInput("selected_comparison",
                "Select Comparison:",
                choices = comparisons[[comp_type]])
  })
  
  selected_data <- reactiveVal()
  
  get_comparison_data <- reactive({
    req(input$comparison_type, input$selected_comparison)
    key <- paste(tolower(input$comparison_type), input$selected_comparison, sep = "_")
    all_deg_data[[key]]
  })
  
  create_volcano_plot <- function(gene_data) {
    # Create color categories
    gene_data$color_category <- "Not Significant"
    gene_data$color_category[gene_data$estimate > input$fc_cutoff & 
                               gene_data$`-Log10(FDR)` > input$p_cutoff] <- "Up"
    gene_data$color_category[gene_data$estimate < -input$fc_cutoff & 
                               gene_data$`-Log10(FDR)` > input$p_cutoff] <- "Down"
    
    # Convert to factor with specific order
    gene_data$color_category <- factor(gene_data$color_category, 
                                       levels = c("Not Significant", "Up", "Down"))
    
    ggplot(gene_data, aes(x = estimate, y = `-Log10(FDR)`)) +
      geom_point(aes(color = color_category),
                 size = input$point_size) +
      geom_hline(yintercept = input$p_cutoff, linetype = "dashed") +
      geom_vline(xintercept = c(-input$fc_cutoff, input$fc_cutoff), 
                 linetype = "dashed") +
      scale_color_manual(values = c("grey", "red", "green3")) +
      labs(x = "Estimate Change",
           y = "-Log10(FDR)",
           title = paste("Volcano Plot -", input$selected_comparison,
                         "\nPlease select a rectangle covering interested data points by left-click and drag.")) +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank()
      )
  }
  
  output$volcanoPlot <- renderPlot({
    req(input$selected_comparison)
    gene_data <- get_comparison_data()
    gene_data[,"-Log10(FDR)"] <- -log10(gene_data$p.value.adjusted)
    create_volcano_plot(gene_data)
  })
  
  output$plot_stats <- renderPrint({
    gene_data <- get_comparison_data()
    gene_data[,"-Log10(FDR)"] <- -log10(gene_data$p.value.adjusted)
    
    # Calculate up-regulated proteins
    up_proteins <- sum(gene_data$estimate > input$fc_cutoff & 
                         gene_data$`-Log10(FDR)` > input$p_cutoff)
    
    # Calculate down-regulated proteins
    down_proteins <- sum(gene_data$estimate < -input$fc_cutoff & 
                           gene_data$`-Log10(FDR)` > input$p_cutoff)
    
    # Calculate total significant proteins
    total_sig <- up_proteins + down_proteins
    
    cat(sprintf("\nTotal Significant Proteins: %d\nUp-regulated Proteins: %d\nDown-regulated Proteins: %d\n",
                total_sig, up_proteins, down_proteins))
  })
  
  output$click_info <- renderPrint({
    if (!is.null(input$plot_click)) {
      gene_data <- get_comparison_data()
      gene_data$neglog10p <- -log10(gene_data$p.value.adjusted)
      
      nearest_point <- nearPoints(gene_data, 
                                  input$plot_click, 
                                  xvar = "estimate",
                                  yvar = "neglog10p",
                                  threshold = 10, 
                                  maxpoints = 1)
      
      if (nrow(nearest_point) > 0) {
        nearest_point[,"-Log10(FDR)"] <- nearest_point$neglog10p
        nearest_point$neglog10p <- NULL
        cat("Selected Protein:\n")
        print(nearest_point)
      }
    }
  })
  
  observe({
    brush <- input$plot_brush
    if (!is.null(brush)) {
      gene_data <- get_comparison_data()
      gene_data$neglog10p <- -log10(gene_data$p.value.adjusted)
      
      selected <- subset(gene_data,
                         estimate >= brush$xmin &
                           estimate <= brush$xmax &
                           neglog10p >= brush$ymin &
                           neglog10p <= brush$ymax)
      
      selected[,"-Log10(FDR)"] <- selected$neglog10p
      selected$neglog10p <- NULL
      selected_data(selected)
    }
  })
  
  output$geneTable <- renderDT({
    if (!is.null(selected_data())) {
      data <- selected_data()
    } else {
      data <- get_comparison_data()
    }
    
    data[,"-Log10(FDR)"] <- -log10(data$p.value.adjusted)
    
    datatable(data,
              options = list(pageLength = 10,
                             scrollX = TRUE),
              rownames = FALSE) %>%
      formatSignif(columns = c("estimate", "p.value.adjusted", "-Log10(FDR)"), digits = 3)
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("volcano_plot_", input$selected_comparison, ".png", sep = "")
    },
    content = function(file) {
      gene_data <- get_comparison_data()
      gene_data[,"-Log10(FDR)"] <- -log10(gene_data$p.value.adjusted)
      plot <- create_volcano_plot(gene_data)
      ggsave(file, plot, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("protein_data_", input$selected_comparison, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(get_comparison_data(), file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
