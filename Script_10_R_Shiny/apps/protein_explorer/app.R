# =============================================================================
# Shiny app: Protein explorer (Script_10)
# =============================================================================
#
# Look up one or more proteins (by gene symbol) and show their estimates and
# FDRs across the multiple-factor comparisons of Script_6: 8 age (Aged versus
# Young), 8 sex (Female versus Male) or 24 region comparisons, as a table and
# a bubble plot. FDR is NA ("Not tested") when the protein was removed by the
# Method 4 observed-value filter in that comparison.
#
# Data (made by Script-10_R_Shiny.R): helper age.RData, helper region.RData,
#   helper sex.RData - df4: estimates (first half of the rows) and FDRs
#   (second half) per gene symbol, in the comparison order listed below;
#   df1 / df6: protein annotations.
#
# Young = 7 weeks; Aged = 7-11 months.
# Run: shiny::runApp("apps/protein_explorer")
# =============================================================================

library(shiny)
library(ggplot2)
library(DT)
library(dplyr)
library(tidyr)
library(ggtext)

# Load data (each helper file contains df1, df4 and df6).
load("helper age.RData")
df4_age <- df4
load("helper region.RData")
df4_region <- df4
load("helper sex.RData")
df4_sex <- df4

# Get available gene symbols
available_gene_symbols <- sort(df1$Gene.symbol)

# Define UI
ui <- fluidPage(
  titlePanel(div("Protein Analysis Dashboard", tags$br(), tags$small("Young = 7 weeks; Aged = 7-11 months"))),
  
  tags$head(
    tags$style(
      HTML("
        th { text-align: center !important; }
        table, th, td { font-size: 16px; }
        .download-btn { margin-top: 10px; }
      ")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 10,
      radioButtons("analysis_type", 
                   "Select Analysis Type:",
                   choices = c("Age Analysis" = "age",
                               "Region Analysis" = "region",
                               "Sex Analysis" = "sex")),
      
      selectizeInput(
        "gene_symbols",
        "Enter Gene Symbols:", 
        choices = available_gene_symbols,
        multiple = TRUE,
        options = list(placeholder = 'Type to search for a gene symbol')
      ),
      downloadButton("downloadData", "Download Selected Protein Data", class = "download-btn"),
      DT::dataTableOutput("geneTable")
    ),
    mainPanel(
      uiOutput("outputArea")
    )
  )
)

# Define server
server <- function(input, output) {
  
  # Render the main search table
  output$geneTable <- DT::renderDataTable({
    datatable(
      df6,
      options = list(
        pageLength = 20, 
        dom = 'ftip',
        scrollY = '300px'
      ),
      selection = 'multiple'
    )
  })
  
  # Reactive expression for data calculation
  calculatedData <- reactive({
    req(input$gene_symbols)
    
    if(input$analysis_type == "region") {
      # Region analysis (24 comparisons)
      region_index_all <- c(
        # Young Male comparisons
        "Young Male: Bottom versus Central",
        "Young Male: Bottom versus Peripheral",
        "Young Male: Bottom versus Top",
        "Young Male: Central versus Peripheral",
        "Young Male: Central versus Top",
        "Young Male: Peripheral versus Top",
        
        # Young Female comparisons
        "Young Female: Bottom versus Central",
        "Young Female: Bottom versus Peripheral",
        "Young Female: Bottom versus Top",
        "Young Female: Central versus Peripheral",
        "Young Female: Central versus Top",
        "Young Female: Peripheral versus Top",
        
        # Aged Male comparisons
        "Aged Male: Bottom versus Central",
        "Aged Male: Bottom versus Peripheral",
        "Aged Male: Bottom versus Top",
        "Aged Male: Central versus Peripheral",
        "Aged Male: Central versus Top",
        "Aged Male: Peripheral versus Top",
        
        # Aged Female comparisons
        "Aged Female: Bottom versus Central",
        "Aged Female: Bottom versus Peripheral",
        "Aged Female: Bottom versus Top",
        "Aged Female: Central versus Peripheral",
        "Aged Female: Central versus Top",
        "Aged Female: Peripheral versus Top"
      )
      
      df_list <- lapply(input$gene_symbols, function(gene_symbol) {
        # Get correct number of estimates and FDR values
        estimates <- df4_region[[gene_symbol]][1:24]
        fdr_values <- df4_region[[gene_symbol]][25:48]
        
        df6 <- data.frame(
          Description = region_index_all,
          Estimate = estimates,
          stringsAsFactors = FALSE
        )
        
        # Add additional columns
        df6 <- df6 %>%
          tidyr::separate(Description, 
                          into = c("Age and Sex group", "Comparison"),
                          sep = ": ", remove = FALSE) %>%
          tidyr::separate(Comparison,
                          into = c("First region", "Second region"),
                          sep = " versus ", remove = FALSE) %>%
          mutate(
            Direction = ifelse(Estimate > 0, 
                               paste(`First region`, ">", `Second region`),
                               paste(`First region`, "<", `Second region`)),
            FDR = fdr_values,
            `Significant with FDR < 0.05` = ifelse(is.na(FDR), "Not tested", ifelse(FDR < 0.05, "Yes", "No")),
            `-Log10(FDR)` = signif((-1) * log10(FDR), 4)
          )
        
        return(list(gene = gene_symbol, data = df6))
      })
      
    } else if(input$analysis_type == "age") {
      # Age analysis (8 comparisons)
      age_vs_young_index_all <- c("Bottom Male", "Bottom Female", "Top Male", "Top Female", 
                                  "Central Male", "Central Female", "Peripheral Male", "Peripheral Female")
      
      df_list <- lapply(input$gene_symbols, function(gene_symbol) {
        df5 <- data.frame(Value = df4_age[[gene_symbol]])
        df6 <- data.frame(
          Description = age_vs_young_index_all, 
          Estimate = df5[1:8, ]
        ) %>%
          mutate(
            Direction = ifelse(Estimate > 0, "Aged > Young", "Aged < Young"),
            FDR = df5[9:16, ],
            `Significant with FDR < 0.05` = ifelse(is.na(FDR), "Not tested", ifelse(FDR < 0.05, "Yes", "No")),
            `-Log10(FDR)` = signif((-1) * log10(FDR), 4)
          )
        
        return(list(gene = gene_symbol, data = df6))
      })
      
    } else {
      # Sex analysis (8 comparisons)
      female_vs_male_index_all <- c("Bottom Aged", "Bottom Young", "Top Aged", "Top Young", 
                                    "Central Aged", "Central Young", "Peripheral Aged", "Peripheral Young")
      
      df_list <- lapply(input$gene_symbols, function(gene_symbol) {
        df5 <- data.frame(Value = df4_sex[[gene_symbol]])
        df6 <- data.frame(
          Description = female_vs_male_index_all, 
          Estimate = df5[1:8, ]
        ) %>%
          mutate(
            Direction = ifelse(Estimate > 0, "Female > Male", "Female < Male"),
            FDR = df5[9:16, ],
            `Significant with FDR < 0.05` = ifelse(is.na(FDR), "Not tested", ifelse(FDR < 0.05, "Yes", "No")),
            `-Log10(FDR)` = signif((-1) * log10(FDR), 4)
          )
        
        return(list(gene = gene_symbol, data = df6))
      })
    }
    
    return(df_list)
  })
  
  # Generate UI elements for each selected protein
  output$outputArea <- renderUI({
    df_list <- calculatedData()
    
    output_list <- lapply(seq_along(df_list), function(i) {
      df_item <- df_list[[i]]
      gene <- df_item$gene
      
      title_suffix <- switch(input$analysis_type,
                             "age" = "Aged versus Young groups",
                             "region" = "Region Comparisons",
                             "sex" = "Female versus Male groups")
      
      tagList(
        h3(HTML(paste0("Results for Protein: ", gene, "<br>", title_suffix))),
        br(),
        DT::dataTableOutput(paste0("geneData", i)),
        br(),
        plotOutput(paste0("bubblePlot", i), 
                   height = if(input$analysis_type == "region") 800 else 400),
        br()
      )
    })
    
    do.call(tagList, output_list)
  })
  
  # Update visualizations
  observe({
    df_list <- calculatedData()
    
    for (i in seq_along(df_list)) {
      local({
        my_i <- i
        gene_symbol <- df_list[[my_i]]$gene
        df6 <- df_list[[my_i]]$data
        
        # Table rendering
        output[[paste0("geneData", my_i)]] <- DT::renderDataTable({
          if(input$analysis_type == "region") {
            rowColors <- ifelse(!is.na(df6$FDR) & df6$FDR < 0.05,
                                ifelse(df6$Estimate > 0, 
                                       '#ffcdd2',  # Light red for positive significant estimates
                                       '#c8e6c9'), # Light green for negative significant estimates
                                       'white')      # White for non-significant rows
                                
            DT::datatable(
              df6,
              options = list(
                pageLength = 24,
                dom = 't',
                columnDefs = list(
                  list(className = 'dt-center', targets = '_all'),
                  list(
                    targets = '_all',
                    createdCell = JS(
                      sprintf(
                        "function(td, cellData, rowData, row, col) {
              $(td).css('background-color', %s[row]);
            }",
                        jsonlite::toJSON(rowColors)
                      )
                    )
                  )
                )
              )
            )
          } else {
            rowColors <- ifelse(!is.na(df6$FDR) & df6$FDR < 0.05,
                                ifelse(df6$Estimate > 0, 
                                       '#ffcdd2',
                                       '#c8e6c9'),
                                       'white')
                                
            datatable(df6,
                      options = list(
                        pageLength = 8,
                        dom = 't',
                        columnDefs = list(
                          list(className = 'dt-center', targets = '_all'),
                          list(
                            targets = '_all',
                            createdCell = JS(
                              sprintf(
                                "function(td, cellData, rowData, row, col) {
                          $(td).css('background-color', %s[row]);
                        }",
                                jsonlite::toJSON(rowColors)
                              )
                            )
                          )
                        )
                      )
            )
          } %>%
            formatSignif(columns = c("Estimate", "FDR", "-Log10(FDR)"), digits = 3)
        })
        
        # Bubble plot
        output[[paste0("bubblePlot", my_i)]] <- renderPlot({
          plot_title <- switch(input$analysis_type,
                               "age" = "Aged versus Young Groups",
                               "region" = "Region Comparisons",
                               "sex" = "Female versus Male Groups")
          
          if (input$analysis_type %in% c("age", "sex")) {
            current_levels <- unique(df6$Description)
            df6$Description <- factor(df6$Description, 
                                      levels = rev(current_levels))
          }
          
          # -log10(FDR) for the color scale (FDR = 0 would give Inf -> NA).
          df6$`-Log10(FDR)` <- as.numeric(as.character(df6$FDR))
          df6$`-Log10(FDR)` <- -log10(df6$`-Log10(FDR)`)
          df6$`-Log10(FDR)`[is.infinite(df6$`-Log10(FDR)`)] <- NA
          
          y_axis_label <- switch(input$analysis_type,
                                 "age" = "Region and Sex group",
                                 "sex" = "Region and Aged group",
                                 "region" = "Comparisons")
          
          ggplot(df6, aes(x = Estimate, y = Description,  
                          color = `-Log10(FDR)`,
                          size = abs(Estimate))) +
            geom_point(alpha = 0.7) +
            scale_size_continuous(range = c(5, 15)) + 
            scale_colour_gradient(low = "blue", high = "red", na.value = "grey50") +
            ggtitle(paste0("Protein: ", gene_symbol, "\n", plot_title)) +
            geom_vline(xintercept = 0, colour = "darkred", size = 0.5, linetype = 2) +
            theme_bw() +
            theme(
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              plot.title = element_text(size = 30, face = "bold"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 20),
              panel.grid.major = element_line(colour = "grey90"),
              panel.grid.minor = element_line(colour = "grey95")
            ) +
            labs(x = "Estimate",
                 y = y_axis_label,
                 size = "Estimate Magnitude")
        })
      })
    }
  })
  
  # Download handler
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("protein_data_", input$analysis_type, "_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      df_list <- calculatedData()
      all_data <- do.call(rbind, lapply(df_list, function(x) {
        cbind(Protein = x$gene, x$data)
      }))
      write.csv(all_data, file, row.names = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server) 