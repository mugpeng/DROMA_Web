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
    features_search_sel$features <- switch(input$select_features1,
                                           "mRNA" = feas_search[feas_search$type %in% "mRNA",]$name,
                                           "meth" = feas_search[feas_search$type %in% "meth",]$name,
                                           "proteinrppa" = feas_search[feas_search$type %in% "proteinrppa",]$name,
                                           "proteinms" = feas_search[feas_search$type %in% "proteinms",]$name,
                                           "cnv" = feas_search[feas_search$type %in% "cnv",]$name,
                                           "drug" = feas_search[feas_search$type %in% "drug",]$name,
                                           "mutation_gene" = feas_search[feas_search$type %in% "mutation_gene",]$name,
                                           "mutation_site" = feas_search[feas_search$type %in% "mutation_site",]$name,
                                           "fusion" = feas_search[feas_search$type %in% "fusion",]$name)
    
    updateSelectizeInput(session = session, 
                         inputId = 'select_specific_feature',
                         label = 'Features Selection:', 
                         choices = features_search_sel$features, 
                         server = TRUE,
                         selected = "")
  })
  
  # Calculate results ----
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
                      value = 50,
                      title = "Processing...")
    
    # Get total number of features to process
    feas_search_sel <- feas_search[feas_search$type %in% input$select_features2,]
    
    used_core <- ifelse(parallel::detectCores()/2 > 8, 8, parallel::detectCores()/2)
    
    results <- BatchFindSigFeaturesPlus(
      feature1_type = input$select_features1,
      feature1_name = input$select_specific_feature,
      feature2_type = input$select_features2,
      data_type = input$data_type,
      tumor_type = input$tumor_type,
      cores = used_core,
      test_top_100 = FALSE
      # progress_callback = progress_callback 
    )
    
    # Update progress to complete
    updateProgressBar(session = session,
                      id = ns("progressBar"),
                      value = 100,
                      title = "Complete!")
    
    # Hide progress box when done
    shinyjs::hide("progressBox")
    
    # Show completion message with actual number of completed features
    if (!is.null(results)) {
      sendSweetAlert(
        session = session,
        title = "Analysis Complete",
        text = sprintf("Processed %d features in %s",
                       nrow(feas_search_sel),  # Use actual number of results
                       format_time(as.numeric(difftime(Sys.time(), 
                                                       progress_vals$start_time, 
                                                       units = "secs")))),
        type = "success"
      )
    }
    # waiter_hide()
    
    results
  })
  
  # Progress outputs
  output$timeEstimate <- renderText({
    progress_vals$estimated_time
  })
  
  # Render volcano plot ----
  p_volcano <- reactive({
    req(results())
    plotMetaVolcano(results(),
                    es_t = 0.2,
                    P_t = 0.001,
                    label = TRUE,
                    top_label_each = 5,
                    title = paste(input$select_features1, input$select_specific_feature,
                                  "vs", input$select_features2))
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