# Modules/DrugFeature.R
# Source required modules and functions
source("Modules/DataAdapter.R")
source("Functions/FeatureLoading.R")

# UI Component
uiDrugFeature <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    fluidRow(
      column(4,
             # Drug Selection Panel
             wellPanel(
               # Direct selectizeInput for drug selection
               selectizeInput(
                 ns("select_drug"), "Drug Selection:", choices = NULL,
                 options = list(
                   placeholder = 'Please select a drug',
                   onInitialize = I('function() { this.setValue(""); }')
                 )),
               
               # Data type and tumor type filters
               selectInput(inputId = ns("filter_data_type"), 
                           "Filter by data type:", 
                           choices = c("All" = "all",
                                       "Cell Lines" = "CellLine",
                                       "Patient-Derived Cells" = "PDC",
                                       "Patient-Derived Organoids" = "PDO",
                                       "Patient-Derived Xenografts" = "PDX"
                           ), selected = "all"),
               
               selectInput(inputId = ns("filter_tumor_type"), 
                           "Filter by tumor type:", 
                           choices = c("All" = "all",
                                       "Aerodigestive Tract Cancer" = "aerodigestive tract cancer",
                                       "Bladder Cancer" = "bladder cancer",
                                       "Breast Cancer" = "breast cancer",
                                       "Cervical Cancer" = "cervical cancer",
                                       "Choriocarcinoma" = "choriocarcinoma",
                                       "Endometrial Cancer" = "endometrial cancer",
                                       "Gastrointestinal Cancer" = "gastrointestinal cancer",
                                       "Haematopoietic/Lymphoid Cancer" = "haematopoietic/lymphoid cancer",
                                       "Kidney Cancer" = "kidney cancer",
                                       "Liver Cancer" = "liver cancer",
                                       "Lung Cancer" = "lung cancer",
                                       "Nasopharyngeal Cancer" = "nasopharyngeal cancer",
                                       "Nervous System Cancer" = "nervous system cancer",
                                       "Non-Cancer" = "non-cancer",
                                       "Ovarian Cancer" = "ovarian cancer",
                                       "Pancreatic Cancer" = "pancreatic cancer",
                                       "Prostate Cancer" = "prostate cancer",
                                       "Retinoblastoma" = "retinoblastoma",
                                       "Sarcoma" = "sarcoma",
                                       "Skin Cancer" = "skin cancer",
                                       "Stomach Cancer" = "stomach cancer",
                                       "Testicular Cancer" = "testicular cancer",
                                       "Thyroid Cancer" = "thyroid cancer",
                                       "Uterine Cancer" = "uterine cancer",
                                       "Vulvar Cancer" = "vulvar cancer"
                           ), selected = "all"),
               
               # Add "Compare by which" selection
               hr(),
               h4("Comparison Options"),
               selectInput(inputId = ns("compare_by"), 
                           "Compare by which:", 
                           choices = c("TumorType", "Gender", "FullEthnicity", "Age"),
                           selected = "TumorType"),
               
               # Remove the checkbox since we're now showing both values
               
               # Download section
               hr(), # Add a horizontal rule for separation
               prettyRadioButtons(
                 inputId = ns("download_type"),
                 label = "Select download format:",
                 choices = c(
                   "Data (RDS)" = "data_rds",
                   "Data (CSV)" = "data_csv"
                 ),
                 selected = "data_csv",
                 icon = icon("check"), 
                 animation = "jelly",
                 status = "primary"
               ), 
               downloadBttn(
                 outputId = ns("download_content"),
                 label = "Download",
                 style = "gradient",
                 color = "default",
                 block = TRUE,
                 size = "sm"
               )
             )
      ),
      
      column(8,
             # Results Panels with tabs for different visualizations
             tabsetPanel(
               tabPanel("Drug Data", 
                        h4("Drug Sensitivity Data"),
                        DT::dataTableOutput(ns("drug_table"))),
               
               tabPanel("Annotated Data", 
                        h4("Drug Sensitivity with Sample Annotations"),
                        DT::dataTableOutput(ns("annotated_drug_table"))),
               
               tabPanel("Comparison", 
                        h4("Drug Sensitivity by Selected Attribute"),
                        uiOutput(ns("comparison_ui")),
                        plotOutput(ns("comparison_plot"), height = "500px")
               )
             )
      )
    )
  )
}

# Server Component
serverDrugFeature <- function(input, output, session) {
  ns <- session$ns
  
  # Update drug selection choices from DROMA database
  observe({
    drugs_available <- getAvailableFeatures("drug")

    updateSelectizeInput(session = session, inputId = 'select_drug',
                         choices = drugs_available, server = TRUE,
                         options = list(placeholder = 'Please select a drug'))
  }, once = TRUE)
  
  # Process drug data using DROMA_R
  processed_data <- reactive({
    shiny::validate(
      shiny::need(input$select_drug != "", "Please select a drug.")
    )

    # Convert filter parameters
    data_type_filter <- if(input$filter_data_type == "all") NULL else input$filter_data_type
    tumor_type_filter <- if(input$filter_tumor_type == "all") NULL else input$filter_tumor_type

    # Use DROMA_R to get drug sensitivity data
    drug_data <- getDrugData(
      drug_name = input$select_drug,
      data_type_filter = data_type_filter,
      tumor_type_filter = tumor_type_filter
    )

    return(drug_data)
  })

  # Merge with sample annotations (DROMA_R already includes annotations)
  annotated_data <- reactive({
    req(processed_data())

    # DROMA_R's getDrugSensitivityData already includes annotations
    # So we just return the processed data
    return(processed_data())
  })
  
  # Drug sensitivity table
  output$drug_table <- DT::renderDataTable({
    req(processed_data())
    
    # Use the extracted function instead of inline code
    format_drug_table(
      processed_data(),
      caption = "Drug sensitivity data - showing both raw and Z-score normalized values"
    )
  })
  
  # Annotated drug table
  output$annotated_drug_table <- DT::renderDataTable({
    req(annotated_data())
    
    # Use the extracted function instead of inline code
    format_drug_table(
      annotated_data(),
      caption = "Drug sensitivity with sample annotations - showing both raw and Z-score normalized values"
    )
  })
  
  # Dynamic UI for comparison settings
  output$comparison_ui <- renderUI({
    req(input$compare_by)
    
    # Check if the comparison variable is continuous (numeric)
    if (is.numeric(annotated_data()[[input$compare_by]])) {
      # For continuous variables
      fluidRow(
        column(6,
               sliderInput(ns("group_bins"), 
                           "Number of group groups:",
                           min = 2, max = 10, value = 4, step = 1)
        ),
        column(6,
               checkboxInput(ns("show_jitter"), "Show Groups Boxplot", value = TRUE)
        )
      )
    } else {
      # For categorical variables - no additional controls needed
      NULL
    }
  })
  
  # Comparison visualization
  output$comparison_plot <- renderPlot({
    req(annotated_data(), input$compare_by)
    data <- annotated_data()
    
    # Always use z-score value for plots
    plot_data <- data
    plot_data$value <- plot_data$zscore_value
    value_type_label <- "Z-score normalized drug sensitivity values"
    
    # Handle missing values in the comparison variable
    plot_data <- plot_data[!is.na(plot_data[[input$compare_by]]), ]
    
    if (nrow(plot_data) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data available for this comparison") + 
               theme_void())
    }
    
    # Different visualization based on the type of the comparison variable
    if (is.numeric(annotated_data()[[input$compare_by]])) {
       return(create_drug_comparison_plot(
        data = plot_data, 
        comparison_var = input$compare_by,
        value_column = "value",
        value_label = value_type_label,
        num_bins = input$group_bins,
        show_groups_boxplot = input$show_jitter
      ))
    } else {
      # For other categorical variables
      return(create_drug_comparison_plot(
        data = plot_data,
        comparison_var = input$compare_by,
        value_column = "value",
        value_label = value_type_label
      ))
    }
  })
  
  # Download handler
  output$download_content <- downloadHandler(
    filename = function() {
      type <- input$download_type
      drug_name_safe <- gsub("[^A-Za-z0-9_]", "_", input$select_drug)
      base_name <- paste0("DrugFeature_", drug_name_safe)
      
      switch(type,
             "data_rds" = paste0(base_name, "_data.rds"),
             "data_csv" = paste0(base_name, "_data.csv")
      )
    },
    content = function(filename) {
      type <- input$download_type
      
      # Determine which data to download - annotated or raw
      data_to_save <- if (exists("sample_anno")) annotated_data() else processed_data()
      
      switch(type,
             "data_rds" = saveRDS(data_to_save, filename),
             "data_csv" = write.csv(data_to_save, filename, row.names = FALSE)
      )
    }
  )
}