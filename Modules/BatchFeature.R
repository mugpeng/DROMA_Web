# Source required modules and functions
source("Modules/DataAdapter.R")
source("Functions/FeatureLoading.R")
source("Functions/ZScoreNormalization.R")
source("Functions/CommonFunctions.R")

uiBatchFeature <- function(id){
  ns <- NS(id)
  fluidPage(
    useWaiter(), # for overall loading
    useSweetAlert(), # for notifications
    
    # Add progress elements
    tags$head(
      tags$style(HTML("
        #progressBox {
          background-color: #f8f9fa;
          padding: 15px;
          border-radius: 5px;
          margin-bottom: 20px;
        }
        .progress {
          margin-bottom: 10px;
        }
      "))
    ),
    fluidRow(
      # Select Features 1 ----
      column(4,
             selectInput(inputId = ns("select_features1"), 
                         "Please select the feature type:", 
                         choices = c("Copy Number Data" = "cnv",
                                     "DNA Methylation" = "meth",
                                     "Gene Fusion" = "fusion",
                                     "Gene Mutation" = "mutation_gene",
                                     "Gene Site Mutation" = "mutation_site",
                                     "mRNA Expression" = "mRNA",
                                     "Protein RPPA Expression" = "proteinrppa",
                                     "Protein MS Expression" = "proteinms",
                                     "Drug Sensitivity" = "drug"
                         ), selected = "proteinms"
             )),
      # Select specific feature ----
      column(4,
             selectizeInput(
               ns("select_specific_feature"), "Features Selection:", choices = NULL,
               options = list(
                 placeholder = 'Please select a feature',
                 onInitialize = I('function() { this.setValue(""); }'), selected = ""
               ))),
      # Select Features 2 ----
      column(4,
             selectInput(inputId = ns("select_features2"), 
                         "Please select the second feature:", 
                         choices = c("Copy Number Data" = "cnv",
                                     "DNA Methylation" = "meth",
                                     "Gene Fusion" = "fusion",
                                     "Gene Mutation" = "mutation_gene",
                                     "Gene Site Mutation" = "mutation_site",
                                     "mRNA Expression" = "mRNA",
                                     "Protein RPPA Expression" = "proteinrppa",
                                     "Protein MS Expression" = "proteinms",
                                     "Drug Sensitivity" = "drug"
                         ), selected = "drug"
             )),
    ),
    # Add data_type and tumor_type filters
    fluidRow(
      column(6,
             selectInput(inputId = ns("data_type"), 
                         "Filter by data type:", 
                         choices = c("All" = "all",
                                     "Cell Lines" = "cell",
                                     "Patient-Derived Cells" = "PDC",
                                     "Patient-Derived Organoids" = "PDO",
                                     "Patient-Derived Xenografts" = "PDX"
                         ), selected = "all"
             )),
      column(6,
             selectInput(inputId = ns("tumor_type"), 
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
                         ), selected = "all"
             ))
    ),
    fluidRow(
      column(12,
             hidden(
               div(id = ns("progressBox"),
                   h4("Analysis Progress"),
                   progressBar(ns("progressBar"), value = 0, 
                               title = "Calculating...",
                               display_pct = TRUE) 
               )
             )
      )
    ),
    # Number of cores selection
    # fluidRow(
    #   column(4,
    #          numericInput(ns("n_cores"), 
    #                       "Number of CPU cores to use:", 
    #                       value = 1, min = 1, max = parallel::detectCores())
    #   )
    # ),
    # Output results ----
    wellPanel(
      column(12,
             tabsetPanel(
               tabPanel("Volcano Plot",
                        plotOutput(ns("volcano_plot"), height = "600px"))
             )),
      h5(".")
    ),
    # Download section
    fluidRow(
      column(6,
             prettyRadioButtons(
               inputId = ns("download_type"),
               label = "Select download format:",
               choices = c(
                 "PDF" = "pdf",
                 "CSV Results" = "csv",
                 "Plot Data" = "plot_data"
               ),
               selected = "pdf",
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
      ),
      column(6,
             h4(strong("NOTEs:")),
             h5("1. The analysis compares one feature against all features in the selected dataset type."),
             h5("2. The xaxis and yaxis are effect size and p value calculated from Meta analysis results."),
             h5("3. Multiple cores can speed up computation but uses more memory.")
      )
    )
  )
}

serverBatchFeature <- function(input, output, session){
  ns <- session$ns
  
  # Track z-score changes
  zscore_tracker <- reactiveVal(if(exists("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)) 
                               base::get("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)$timestamp 
                               else Sys.time())
  
  # Update tracker when global state changes
  observe({
    if(exists("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)) {
      zscore_tracker(base::get("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)$timestamp)
    }
    invalidateLater(1000) # Check every second
  })
  
  # Add reactive values for progress tracking
  progress_vals <- reactiveValues(
    start_time = NULL,
    features_done = 0,
    total_features = 0,
    estimated_time = "Calculating..."
  )
  
  # Select Features reactive ----
  features_search_sel <- reactiveValues()
  
  observeEvent(input$select_features1, {
    # Get available features from DROMA database
    features <- getAvailableFeatures(input$select_features1)

    features_search_sel$features <- features

    updateSelectizeInput(session = session,
                         inputId = 'select_specific_feature',
                         label = 'Features Selection:',
                         choices = features,
                         server = TRUE,
                         selected = "")
  }, ignoreNULL = FALSE)
  
  # Calculate results using DROMA_R ----
  results <- reactive({
    # React to z-score changes
    zscore_tracker()

    req(input$select_features1, input$select_features2, input$select_specific_feature)

    # Show progress box
    shinyjs::show("progressBox")
    progress_vals$start_time <- Sys.time()
    progress_vals$features_done <- 0

    # Simple progress indication
    updateProgressBar(session = session,
                      id = ns("progressBar"),
                      value = 10,
                      title = "Initializing...")

    # Convert filter parameters
    data_type_filter <- if(input$data_type == "all") NULL else input$data_type
    tumor_type_filter <- if(input$tumor_type == "all") NULL else input$tumor_type

    # Get number of cores to use
    used_core <- ifelse(parallel::detectCores()/2 > 8, 8, parallel::detectCores()/2)

    # Use DROMA_R batch feature analysis
    updateProgressBar(session = session,
                      id = ns("progressBar"),
                      value = 30,
                      title = "Running analysis...")

    results <- batchFindSignificantFeaturesWrapper(
      feature1_type = input$select_features1,
      feature2_type = input$select_features2,
      feature1_names = input$select_specific_feature,
      feature2_names = NULL,  # Analyze all features of type2
      data_type_filter = data_type_filter,
      tumor_type_filter = tumor_type_filter,
      cores = used_core
    )

    # Update progress to complete
    updateProgressBar(session = session,
                      id = ns("progressBar"),
                      value = 100,
                      title = "Complete!")

    # Hide progress box when done
    shinyjs::hide("progressBox")

    # Show completion message
    if (!is.null(results) && nrow(results) > 0) {
      sendSweetAlert(
        session = session,
        title = "Analysis Complete",
        text = sprintf("Found %d significant associations in %s",
                       nrow(results),
                       format_time(as.numeric(difftime(Sys.time(),
                                                       progress_vals$start_time,
                                                       units = "secs")))),
        type = "success"
      )
    } else {
      sendSweetAlert(
        session = session,
        title = "Analysis Complete",
        text = "No significant associations found",
        type = "info"
      )
    }

    results
  })
  
  # Progress outputs
  output$timeEstimate <- renderText({
    progress_vals$estimated_time
  })
  
  # Render volcano plot using DROMA_R ----
  p_volcano <- reactive({
    req(results())

    # Use DROMA_R plotMetaVolcano function if available
    # For now, create a basic volcano plot
    if (nrow(results()) > 0) {
      tryCatch({
        # Try to use DROMA_R's volcano plot function
        DROMA.R::plotMetaVolcano(
          results(),
          es_t = 0.2,
          P_t = 0.001,
          label = TRUE,
          top_label_each = 5,
          title = paste(input$select_features1, input$select_specific_feature,
                        "vs", input$select_features2)
        )
      }, error = function(e) {
        # Fallback to basic ggplot if DROMA_R function not available
        library(ggplot2)
        df <- results()
        ggplot(df, aes(x = effect_size, y = -log10(p_value))) +
          geom_point(alpha = 0.5) +
          geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "red") +
          geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "red") +
          theme_minimal() +
          labs(title = paste(input$select_features1, input$select_specific_feature,
                            "vs", input$select_features2),
               x = "Effect Size",
               y = "-log10(p-value)")
      })
    } else {
      # Return empty plot if no results
      library(ggplot2)
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No significant associations found") +
        theme_void()
    }
  })
  
  output$volcano_plot <- renderPlot({
    p_volcano()
  })
  
  # Render results table ----
  # output$result_table <- DT::renderDataTable({
  #   req(results())
  #   DT::datatable(results(),
  #                 options = list(scrollX = TRUE),
  #                 selection = 'single')
  # })
  
  # Download handler ----
  output$download_content <- downloadHandler(
    filename = function() {
      type <- input$download_type
      base_name <- paste0(input$select_features1, "_", 
                          input$select_specific_feature, "_vs_",
                          input$select_features2)
      
      switch(type,
             "pdf" = paste0(base_name, "_volcano.pdf"),
             "csv" = paste0(base_name, "_results.csv"),
             "plot_data" = paste0(base_name, "_data.rds"))
    },
    content = function(file) {
      type <- input$download_type
      
      switch(type,
             "pdf" = {
               ggsave(file, 
                      plot = p_volcano(),
                      width = 10, height = 8)
             },
             "csv" = {
               write.csv(results(), file, row.names = FALSE)
             },
             "plot_data" = {
               saveRDS(results(), file)
             })
    }
  )
}