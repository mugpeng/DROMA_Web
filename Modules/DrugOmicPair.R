uiDrugOmicPair <- function(id){
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(4,
             selectInput(inputId = ns("select_omics"), 
                         "Please select the molecular type:", 
                         choices = c("Copy Number Data" = "cnv",
                                     "DNA Methylation" = "meth",
                                     "mRNA Expression" = "mRNA",
                                     "Protein RPPA Expression" = "proteinrppa",
                                     "Protein MS Expression" = "proteinms",
                                     "Gene Fusion" = "fusion",
                                     "Gene Mutation" = "mutation_gene",
                                     "Gene Site Mutation" = "mutation_site"
                         ), selected = "mRNA"
             )),
      # Select omics ----
      column(4,
             selectizeInput(
               ns("select_specific_omic"), "Molecule Selection:", choices = NULL,
               options = list(
                 placeholder = 'Please select a molecular feature',
                 onInitialize = I('function() { this.setValue(""); }'), selected = "ABCC3"
               )),
             # DOPselectOmicsUI("DOPselectOmics"),
      ),
      # Select drugs ----
      column(4,
             selectizeInput(
               ns("select_specific_drug"), "Drug Selection:", choices = NULL,
               options = list(
                 placeholder = 'Please select a drug',
                 onInitialize = I('function() { this.setValue(""); }'), selected = "Sepantronium bromide"
               ))
      ),
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
    # plot results ----
    wellPanel(
      # textOutput("total")
      column(12, plotOutput(ns("patchPlot"), height="20cm")),
      column(12, plotOutput(ns("metaPlot"), height="8cm")),
      h5("."),
    ),
    # Add the new download section
    fluidRow(
      column(6,
             h4(strong("NOTEs:")),
             h5("Given that there are overlapping cells in different drug and omics datasets, we have utilized the common data to assess correlations, thereby maximizing the utilization of existing information. For instance, the designation 'gdsc_ctrp1' indicates that the omics data is sourced from the GDSC project, while the drug sensitivity data is derived from the CTRP1 project.")
      ),
      column(6,
             prettyRadioButtons(
               inputId = ns("download_type"),
               label = "Select download format:",
               choices = c(
                 "PDF" = "pdf",
                 "ggplot Obj" = "plot_obj",
                 "data Rds" = "plot_data",
                 "Meta Forest" = "meta"
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
             ),
             h5("PDF: saving above plot.pdf."),
             h5("ggplot Obj: patchwork ggplot object, you can adjust it by yourself."),
             h5("data Rds: the list object used to plot."),
             h5("Meta Forest: saveing meta plot.pdf")
      ),
    )
  )
}

serverDrugOmicPair <- function(input, output, session){
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
  
  # Select ----
  ## Omics ----
  omics_search_sel <- reactiveValues()
  observeEvent(input$select_omics, {
    omics_search_sel$omics <- switch(input$select_omics, 
                                     "mRNA" = omics_search[omics_search$type %in% "mRNA",]$omics,
                                     "meth" = omics_search[omics_search$type %in% "meth",]$omics,
                                     "proteinrppa" = omics_search[omics_search$type %in% "proteinrppa",]$omics,
                                     "proteinms" = omics_search[omics_search$type %in% "proteinms",]$omics,
                                     "cnv" = omics_search[omics_search$type %in% "cnv",]$omics,
                                     "mutation_gene" = omics_search[omics_search$type %in% "mutation_gene",]$omics,
                                     "mutation_site" = omics_search[omics_search$type %in% "mutation_site",]$omics,
                                     "fusion" = omics_search[omics_search$type %in% "fusion",]$omics)
    updateSelectizeInput(session = session, inputId = 'select_specific_omic',
                         label = 'Molecule Selection:', choices = omics_search_sel$omics, server = TRUE,
                         options = list(placeholder = 'Please select a molecular feature', onInitialize = I('function() { this.setValue(""); }')),
                         selected = "ABCC3"
    )
  })
  
  ## Drugs ----
  updateSelectizeInput(session = session, inputId = 'select_specific_drug',
                       label = 'Drug Selection:', choices = drugs_search$drugs, server = TRUE,
                       options = list(placeholder = 'Please select a drug', onInitialize = I('function() { this.setValue(""); }')),
                       selected = "Sepantronium bromide"
  )
  
  # Calculate pair result and plot
  selected_obj <- reactive({
    # React to z-score changes
    zscore_tracker()
    
    # Check if drug is selected
    shiny::validate(
      shiny::need(input$select_specific_drug != "", "Please select a drug.")
    )

    # Check if omic is selected
    shiny::validate(
      shiny::need(input$select_specific_omic != "", "You are not chosen omic yet.")
    )

    # Check if z-score normalization is enabled
    merged_enabled <- FALSE
    if(exists("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)) {
      merged_enabled <- isTRUE(base::get("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)$enabled)
    }
    
    # Use the oneDrugOmicPair function which now uses our refactored functions internally
    oneDrugOmicPair(input$select_omics, input$select_specific_omic,
                   input$select_specific_drug,
                   data_type = input$data_type, 
                   tumor_type = input$tumor_type,
                   merged_enabled = merged_enabled)
  })
  
  output$patchPlot <- renderPlot({
    # modidy plot
    plot_with_axis <- create_plot_with_common_axes(selected_obj()$plot, 
                                                   x_title = "Molecular State(mRNA expression or Mutation status)", 
                                                   y_title = "drug sensitivity (higher indicates resistance)")
    plot_with_axis()
  })
  
  output$metaPlot <- renderPlot({
    req(selected_obj()$meta)  
    create_forest_plot(selected_obj()$meta)
  })
  
  # Add download handler
  output$download_content <- downloadHandler(
    filename = function() {
      type <- input$download_type
      omic_drug_name <- paste0(input$select_omics, "_", input$select_specific_omic, "_", input$select_specific_drug)
      switch(type,
             "pdf" = paste0("plot-", omic_drug_name, ".pdf"),
             "plot_obj" = paste0("plot_obj-", omic_drug_name, ".rds"),
             "plot_data" = paste0("plot_data-", omic_drug_name, ".rds"),
             "meta" = paste0("meta-", omic_drug_name, ".pdf")
      )
    },
    content = function(filename) {
      req(selected_obj())
      type <- input$download_type
      switch(type,
             "pdf" = {
               # Save as PDF
               ggplot2::ggsave(
                 filename = filename, plot = selected_obj()$plot, device = "pdf",
                 units = "cm", width = 15, height = 20, dpi = 600
               )
             },
             "plot_obj" = {
               # Save plot object
               saveRDS(selected_obj()$plot, filename)
             },
             "plot_data" = {
               # Save data
               saveRDS(selected_obj()$data, filename)
             },
             "meta" = {
               pdf(file = filename, width = 10, height = 6)
               create_forest_plot(selected_obj()$meta)
               dev.off()
             }
      )
    }
  )
}
