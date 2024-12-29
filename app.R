library(shiny)
library(tidyverse)
library(DT)
library(ggplot2)
library(colourpicker)
library(shinythemes)
library(pheatmap)
library(stringr)
library(RColorBrewer)

# Set maximum upload size for files to 30 MB
options(shiny.maxRequestSize = 30 * 1024^2)

#---------------------------------- UI ----------------------------------
ui <- fluidPage(
  theme = shinytheme("yeti"),
  titlePanel("BF591 Final Project"), # Add a title for the application
  p("Anuradha Basyal"),
  
  tabsetPanel(
    #---------------------- Samples Tab ----------------------
    tabPanel("Samples",
             p("Use this tab to get distinct values and distributions of sample information. This tab allows user to load and examine a sample information matrix."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("samples_csvfile", "Upload Sample Metadata:", accept = ".csv"),   #only allows uploading .csv file, all others wont be able to upload
                 helpText("Upload a well-formatted CSV file containing metadata."),
                 actionButton("submit_samples", "Submit", width = '100%')
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", tableOutput("samples_summary")),
                   tabPanel("Table", dataTableOutput("samples_table")),
                   tabPanel("Plots",
                            selectInput("plot_column", "Select Column for Plotting:", choices = NULL),
                            selectInput("group_column", "Select Grouping Column (Optional):", choices = NULL),
                            selectInput("plot_type", "Select Plot Type:", choices = c("Histogram", "Density", "Violin")),
                            actionButton("generate_plot", "Generate Plot", width = '100%'),
                            plotOutput("samples_plots")
                   )
                 )
               )
             )
    ),
    
    #---------------------- Counts Matrix Tab ----------------------
    tabPanel("Counts",
             p("Use this tab to explore, filter, and visualize the counts matrix based on statistical properties like variance and non-zero values. It provides tools to view diagnostic plots, generate heatmaps, and perform PCA analysis."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("counts_csvfile", "Upload Normalized Counts Matrix:", accept = ".csv"),  #only allows uploading .csv file, all others wont be able to upload
                 helpText("Upload a normalized gene expression counts matrix."),
                 sliderInput("variance_slider", "Minimum Variance Percentile:", min = 0, max = 100, value = 50),
                 sliderInput("nonzero_slider", "Minimum Non-Zero Samples Percentile:", min = 0, max = 100, value = 50),
                 actionButton("apply_filters", "Apply Filters", width = '100%')
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Filter Summary", tableOutput("filter_summary")),
                   tabPanel("Diagnostic Plots", 
                            plotOutput("variance_vs_median"), 
                            plotOutput("nonzero_vs_median")
                   ),
                   tabPanel("Heatmap",
                            p("Clustered heatmap of filtered counts matrix."),
                            plotOutput("heatmap", height = "600px", width = "800px")
                   ),
                   tabPanel("PCA",
                            selectInput("pca_x", "Select PC for X-axis:", choices = NULL),
                            selectInput("pca_y", "Select PC for Y-axis:", choices = NULL),
                            actionButton("generate_pca", "Generate PCA Plot", width = '100%'),
                            plotOutput("pca_plot")
                   )
                 )
               )
             )
    ),
    
    #---------------------- Differential Expression Tab ----------------------
    tabPanel("Differential Expression",
             p("Use this tab to upload and analyze differential expression analysis (DEA) results. It includes visualizations like a Volcano Plot to examine statistical significance and effect size and a data table for browsing the results.

"),
             sidebarLayout(
               sidebarPanel(
                 fileInput("DE_Anaylsis", "Upload DEA Results CSV:", accept = ".csv"),   #only allows uploading .csv file, all others wont be able to upload
                 helpText("Upload the results of differential expression analysis."),
                 selectInput("x_axis", "X-axis", 
                             choices = c("log2FoldChange", "baseMean", "lfcSE", "stat", "pvalue", "padj"),
                             selected = "log2FoldChange"),
                 selectInput("y_axis", "Y-axis", 
                             choices = c("log2FoldChange", "baseMean", "lfcSE", "stat", "pvalue", "padj"), selected = "padj"),
                 sliderInput("padjusted", "P-Adjusted Threshold:", min = -50, max = 0, value = -15),
                 colourInput("de_base", "Base Point Color", "#FF0000"),
                 colourInput("de_highlight", "Highlight Point Color", "#000000")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Volcano Plot", plotOutput("volcano", height = "500px")),
                   tabPanel("Table", DT::dataTableOutput("table"))
                 )
               )
             )
    ),
    
    #---------------------- GSEA Tab ----------------------
    tabPanel("GSEA",
             p("Use this tab to load and visualize Gene Set Enrichment Analysis (GSEA) results. It includes tools for creating bar plots, scatter plots, and filtering tables based on significance thresholds."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("GSEA_csvfile", "Upload GSEA Results Table:", accept = c(".csv")),  #only allows uploading .csv file, all others wont be able to upload
                 helpText("Upload a CSV containing GSEA results."),
                 downloadButton("download_NES_table", "Download Filtered Table")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("NES Bar Plot",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("GSEA_Bar", "Top pathways shown by adjusted p-value", min = 5, max = 50, value = 10)
                              ),
                              mainPanel(plotOutput("GSEA_barplot"))
                            )
                   ),
                   tabPanel("Table",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("GSEA_table_pvalue", "P-Adjusted Threshold:", min = -20, max = 0, value = -10),
                                selectInput("GSEA_pathways_choice", "Filter Pathways:", choices = c("All", "Positive", "Negative"))
                              ),
                              mainPanel(DT::dataTableOutput("GSEA_results_table"))
                            )
                   ),
                   tabPanel("NES Scatter Plot",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("GSEA_scatter_pvalue", "P-Adjusted Threshold for Scatter Plot:", min = -20, max = 0, value = -12)
                              ),
                              mainPanel(plotOutput("GSEA_scatterplot"))
                            )
                   )
                 )
               )
             )
    )
  )
)


#---------------------------------- SERVER ----------------------------------
server <- function(input, output, session) {
  
  #----------------------Samples Logic----------------------
  # Reactive function to load sample data
  samples_load_data <- reactive({
    req(input$samples_csvfile)
    file <- input$samples_csvfile
    validate(need(file, "Please upload a valid CSV file."))  #only allows uploading .csv file, all others wont be able to upload
    
    # Attempt to load the file
    tryCatch(
      {
        datafile <- read_csv(file$datapath)
        return(datafile)
      },
      error = function(e) {
        showNotification("Error reading the file. Ensure it's a well-formatted CSV.", type = "error")
        return(NULL)
      }
    )
  })
  
  # Function to generate summary table
  samples_table_summary <- function(dataf) {
    # Identify columns
    column <- colnames(dataf)
    column_type <- sapply(dataf, function(x) {
      if (is.numeric(x)) "double" else "factor"
    })
    
    # Generate numeric summaries or distinct values
    summary_stat <- sapply(dataf, function(x) {
      if (is.numeric(x)) {
        sprintf("%.2f (+/- %.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
      } else {
        NA
      }
    })
    
    # Generate distinct values for categorical columns
    distinct_values <- sapply(dataf, function(x) {
      if (!is.numeric(x)) paste(unique(x), collapse = ", ") else NA
    })
    
    # Combine into a formatted table
    df <- tibble(
      `Column Name` = column,
      `Type` = column_type,
      `Mean (sd) or Distinct Values` = ifelse(column_type == "double", summary_stat, distinct_values)
    )
    
    return(df)
  }
  
  
  # Function to generate plots dynamically based on selected type
  samples_plot <- function(dataf, plot_column, group_column = NULL, plot_type) {
    validate(need(plot_column, "Please select a column to plot."))
    validate(need(is.numeric(dataf[[plot_column]]), "Selected column is not numeric."))
    
    p <- ggplot(dataf, aes_string(x = plot_column, fill = group_column))
    
    if (plot_type == "Histogram") {
      bin_width <- (max(dataf[[plot_column]], na.rm = TRUE) - min(dataf[[plot_column]], na.rm = TRUE)) / 30
      p <- p + geom_histogram(binwidth = bin_width, alpha = 0.7, position = "identity", color = "black")
    } else if (plot_type == "Density") {
      p <- p + geom_density(alpha = 0.5)
    } else if (plot_type == "Violin") {
      validate(need(group_column, "Please select a grouping column for violin plots."))
      p <- ggplot(dataf, aes_string(x = group_column, y = plot_column)) +
        geom_violin(fill = "pink", alpha = 0.5)
    }
    
    p <- p + labs(title = paste(plot_type, "Plot:"), x = plot_column, y = plot_type) + theme_minimal()
    
    return(p)
  }
  
  # Rendering table summary
  output$samples_summary <- renderTable({
    data <- samples_load_data()
    req(data)
    samples_table_summary(data)
  })
  
  # Rendering data table
  output$samples_table <- renderDataTable({
    data <- samples_load_data()
    req(data)
    datatable(data, options = list(pageLength = 10))
  })
  
  # Dynamically update column choices for plots
  observeEvent(input$submit_samples, {
    data <- samples_load_data()
    req(data)
    numeric_cols <- colnames(data)[sapply(data, is.numeric)]
    updateSelectInput(session, "plot_column", choices = numeric_cols)
    
    group_cols <- colnames(data)[!sapply(data, is.numeric)]
    updateSelectInput(session, "group_column", choices = c("None", group_cols))
  })
  
  # Rendering plots dynamically
  output$samples_plots <- renderPlot({
    req(input$generate_plot)
    data <- samples_load_data()
    req(data)
    samples_plot(data, input$plot_column, ifelse(input$group_column == "None", NULL, input$group_column), input$plot_type)
  })
  
  #-----------------------Counts Logic------------------
  # Reactive function to load counts matrix
  counts_load_data <- reactive({
    req(input$counts_csvfile)
    file <- input$counts_csvfile
    validate(need(file, "Please upload a valid CSV file.")) #only allows uploading .csv file, all others wont be able to upload
    tryCatch(
      {
        datafile <- read_csv(file$datapath)
        return(datafile)
      },
      error = function(e) {
        showNotification("Error reading the file. Ensure it's a well-formatted CSV.", type = "error")
        return(NULL)
      }
    )
  })
  
  # Reactive function for filtered counts
  filtered_counts <- reactive({
    req(counts_load_data())
    data <- counts_load_data()
    
    # Calculate variance and non-zero values
    data$variance <- apply(data[,-1], 1, var)
    data$nonzero <- rowSums(data[,-1] > 0)
    
    # Filter based on sliders
    data <- data %>% mutate(
      filter_status = variance >= quantile(variance, input$variance_slider / 100) &
        nonzero >= input$nonzero_slider
    )
    return(data)
  })
  
  # Render filter summary
  output$filter_summary <- renderTable({
    req(filtered_counts())
    original <- counts_load_data()
    filtered <- filtered_counts() %>% filter(filter_status)
    
    tibble(
      `Total Samples` = ncol(original) - 1,
      `Total Genes` = nrow(original),
      `Genes Passing Filter` = nrow(filtered),
      `Percent Passing` = round((nrow(filtered) / nrow(original)) * 100, 2),
      `Genes Filtered Out` = nrow(original) - nrow(filtered),
      `Percent Filtered Out` = round(100 - (nrow(filtered) / nrow(original)) * 100, 2)
    )
  })
  
  # Render diagnostic scatter plots
  output$variance_vs_median <- renderPlot({
    req(filtered_counts())
    data <- filtered_counts()
    
    # Calculate median counts and log10(variance + 1)
    median_counts <- apply(data[,-1], 1, median)
    log_variance <- log10(data$variance + 1)
    
    ggplot(data, aes(x = median_counts, y = log_variance, color = filter_status)) +
      geom_point(alpha = 0.7) +
      scale_x_log10() +  # Apply log scale to x-axis (median counts)
      scale_color_manual(values = c("grey", "black"), labels = c("Filtered Out", "Passing Filter")) +
      labs(
        title = "Variance vs Median Counts",
        x = "Median Counts (Log Scale)",
        y = "Log10(Variance + 1)",
        color = "Filter Status"
      ) +
      theme_minimal()
  })
  
  output$nonzero_vs_median <- renderPlot({
    req(filtered_counts())
    data <- filtered_counts()
    
    # Calculate median counts and raw nonzero counts
    median_counts <- apply(data[,-1], 1, median)
    raw_nonzero <- data$nonzero
    
    # Plot with linear y-scale and jitter for better visualization
    ggplot(data, aes(x = median_counts, y = log10(raw_nonzero + 1), color = filter_status)) +
      geom_point(alpha = 0.7, position = position_jitter(width = 0.2, height = 0.05)) +
      scale_x_log10() +  # Log scale for median counts
      scale_color_manual(values = c("grey", "black"), labels = c("Filtered Out", "Passing Filter")) +
      labs(
        title = "Non-Zero Samples vs Median Counts",
        x = "Median Counts (Log Scale)",
        y = "Log10(Non-Zero Samples + 1)",
        color = "Filter Status"
      ) +
      theme_minimal()
  })
  
  output$heatmap <- renderPlot({
    req(filtered_counts())  # Ensure filtered counts data exists
    
    # Filter to keep only genes passing filter
    data <- filtered_counts() %>% filter(filter_status)
    
    # Exclude 'variance', 'nonzero', and 'filter_status' columns
    data <- data %>% select(-variance, -nonzero, -filter_status)
    
    # Convert the data to a matrix (exclude the gene ID column)
    mat <- as.matrix(data[, -1])  # Exclude the first column (gene names)
    
    # Set gene names as rownames
    rownames(mat) <- make.unique(as.character(data[[1]]))
    
    # Define a darker color palette with more contrast
    heatmap_color <- colorRampPalette(c("darkblue", "white", "darkred"))(200)
    
    # Define custom breakpoints for better color scaling
    breaks <- seq(-6, 6, length.out = 201)
    
    # Generate heatmap
    pheatmap(
      log10(mat + 1),  # Log-transform the matrix
      scale = "row",   # Normalize rows (genes)
      color = heatmap_color,
      breaks = breaks,  # Use custom breakpoints
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE,  # Optional: Hide row names
      show_colnames = TRUE,
      legend = TRUE
    )
  }, height = 500, width = 800)
  
  
  
  # Render PCA plot
  observeEvent(input$apply_filters, {
    data <- counts_load_data()
    pca <- prcomp(t(data[,-1]), center = TRUE, scale. = TRUE)
    updateSelectInput(session, "pca_x", choices = colnames(pca$x))
    updateSelectInput(session, "pca_y", choices = colnames(pca$x))
  })
  
  output$pca_plot <- renderPlot({
    req(filtered_counts())  # Ensure filtered counts data is available
    data <- counts_load_data()  # Load counts matrix
    
    # Perform PCA analysis
    pca <- prcomp(t(data[,-1]), center = TRUE, scale. = TRUE)
    
    # Calculate variance explained by each principal component
    variance_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100
    x_label <- paste0(input$pca_x, " (", round(variance_explained[as.numeric(gsub("PC", "", input$pca_x))], 2), "%)")
    y_label <- paste0(input$pca_y, " (", round(variance_explained[as.numeric(gsub("PC", "", input$pca_y))], 2), "%)")
    
    # Generate the PCA plot with updated labels
    ggplot(as.data.frame(pca$x), aes_string(x = input$pca_x, y = input$pca_y)) +
      geom_point(alpha = 0.7) +
      labs(
        title = "PCA Plot",
        x = x_label,
        y = y_label
      ) +
      theme_minimal()
  })
  
  
  #----------------------DEA Logic----------------------
  DE_data <- reactive({
    req(input$DE_Anaylsis)
    read_csv(input$DE_Anaylsis$datapath)   #only allows uploading .csv file, all others wont be able to upload
  })
  
  output$volcano <- renderPlot({
    ggplot(DE_data(), aes(x = !!sym(input$x_axis), y = -log10(!!sym(input$y_axis)))) +
      geom_point(aes(color = !!sym(input$y_axis) < 10^input$padjusted)) +
      scale_color_manual(values = c(input$de_base, input$de_highlight)) +
      theme_minimal()
  })
  
  output$table <- renderDataTable({
    datatable(DE_data())
  })
  
  #----------------------GSEA Logic----------------------
  
  # Reactive function to load GSEA data
  GSEA_load_data <- reactive({
    req(input$GSEA_csvfile)
    read_csv(input$GSEA_csvfile$datapath) # Only allow uploading .csv files
  })
  
  # Function to create a bar plot for top pathways
  GSEA_bar <- function(dataf, pathways_slider) {
    filtered <- dataf %>%
      arrange(padj) %>%
      slice_head(n = pathways_slider) %>%
      mutate(
        pathway = str_trunc(pathway, 50),  # Truncate pathway names to 50 characters
        status = ifelse(NES > 0, "Positive", "Negative")
      )
    
    ggplot(filtered, aes(x = reorder(pathway, NES), y = NES, fill = status)) +
      geom_col(width = 0.7) +
      coord_flip() +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.y = element_text(size = 11, hjust = 1),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_fill_manual(values = c("Positive" = "#FF6F61", "Negative" = "#6FA3FF")) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      labs(
        title = "Top Pathways by Normalized Enrichment Score (NES)",
        x = NULL,
        y = "NES",
        fill = "Status"
      )
  }
  
  # Function to filter and prepare the GSEA table
  GSEA_table <- function(dataf, pvalue_choice, type_pathways_choice) {
    filtered_data <- dataf %>%
      mutate(status = ifelse(NES > 0, "Positive", "Negative")) %>%
      filter(padj <= 10^pvalue_choice)
    
    if (type_pathways_choice != "All") {
      filtered_data <- filtered_data %>% filter(status == type_pathways_choice)
    }
    
    filtered_data
  }
  
  # Function to create a scatter plot for GSEA results
  GSEA_scatter <- function(dataf, pvalue_choice) {
    filtered <- dataf %>%
      mutate(filter_status = ifelse(padj < 10^pvalue_choice, "Passed", "Filtered"))
    
    ggplot(filtered, aes(x = NES, y = -log10(padj), color = filter_status)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(
        values = c("Passed" = "#1f77b4", "Filtered" = "grey"),
        labels = c("Below Threshold", "Above Threshold")
      ) +
      theme_light(base_size = 14) +
      labs(
        title = "Scatter Plot: NES vs -log10 Adjusted p-value",
        x = "Normalized Enrichment Score (NES)",
        y = "-log10 Adjusted p-value",
        color = "Filter Status"
      ) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "bottom"
      )
  }
  
  # Render bar plot
  output$GSEA_barplot <- renderPlot({
    req(GSEA_load_data())
    GSEA_bar(GSEA_load_data(), input$GSEA_Bar)
  }, height = reactive(50 * input$GSEA_Bar))  # Adjust height dynamically
  
  # Render results table
  output$GSEA_results_table <- renderDataTable({
    req(GSEA_load_data())
    GSEA_table(
      GSEA_load_data(),
      input$GSEA_table_pvalue,
      input$GSEA_pathways_choice
    )
  })
  
  # Render scatter plot
  output$GSEA_scatterplot <- renderPlot({
    req(GSEA_load_data())
    GSEA_scatter(GSEA_load_data(), input$GSEA_scatter_pvalue)
  })
  
  # Download handler for GSEA results
  output$download_NES_table <- downloadHandler(
    filename = function() { "NES_Results.csv" },
    content = function(file) {
      write.csv(
        GSEA_table(GSEA_load_data(), input$GSEA_table_pvalue, input$GSEA_pathways_choice),
        file,
        row.names = FALSE
      )
    }
  )
  
}

#---------------------------------- Run App ----------------------------------
shinyApp(ui = ui, server = server)








