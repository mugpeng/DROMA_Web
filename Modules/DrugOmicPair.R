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
                                     "Patient-Derived Organoids" = "PDO"
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
    # Plot results ----
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
             h5("The y axis is the drug sensitivity metric."),
             h5("The x-axis represents the molecular state or numerical data, including the boolean state for gene mutations or mRNA expression."),
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
  # Produce plot ----
  # Select drug and omic
  selected_drug <- reactive({
    shiny::validate(
      shiny::need(input$select_specific_drug != "", "You are not chosen drug yet.")
    )
    selFeatures("drug", input$select_specific_drug, 
                data_type = input$data_type, 
                tumor_type = input$tumor_type)
  })
  selected_omic <- reactive({
    shiny::validate(
      shiny::need(input$select_specific_omic != "", "You are not chosen omic yet.")
    )
    selFeatures(input$select_omics, input$select_specific_omic,
                data_type = input$data_type,
                tumor_type = input$tumor_type) 
  })
  # Calculate pair result and plot
  selected_obj <- reactive({
    if (input$select_omics %in% c("mRNA", "meth", "cnv",
                                  "proteinrppa", "proteinms")) {
      selected_pair <- pairDrugOmic(selected_omic(), selected_drug())
      re <- plotDrugOmicPair_con(selected_pair)    
    } else{
      selected_pair <- pairDrugOmic2(selected_omic(), selected_drug()) 
      re <- plotDrugOmicPair_dis(selected_pair)     
    }
    if(length(re)>1){
      return(
        list(plot = re[[1]], meta = re[[2]], data = selected_pair)
      )
    } else {
      return(list(plot = re[[1]], data = selected_pair))
    }
  })
  output$patchPlot <- renderPlot({
    selected_obj()$plot
  })
  output$metaPlot <- renderPlot({
    req(selected_obj()$meta)  
    p_val <- selected_obj()$meta$pval.random
    p_text <- if(p_val < 0.001) {
      paste("Random-Effects Model (p =", format(p_val, scientific = TRUE, digits = 3), ")")
    } else {
      paste("Random-Effects Model (p =", round(p_val, 3), ")")
    }
    meta::forest(selected_obj()$meta, 
                 xlab = "Effect Size (95% CI)", 
                 slab = "study", 
                 print.pval.common = T,
                 boxsize = 0.2, 
                 lineheight = "auto",
                 print.pval.Q=FALSE,
                 print.I2 = F,
                 # resid.hetstat = F,
                 print.tau2	= F,
                 common = FALSE,
                 text.random = p_text
    )
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
               p_val <- selected_obj()$meta$pval.random
               p_text <- if(p_val < 0.001) {
                 paste("Random-Effects Model (p =", format(p_val, scientific = TRUE, digits = 3), ")")
               } else {
                 paste("Random-Effects Model (p =", round(p_val, 3), ")")
               }
               pdf(file = filename, width = 10, height = 6)
               meta::forest(selected_obj()$meta, 
                            xlab = "Effect Size (95% CI)", 
                            slab = "study", 
                            print.pval.common = T,
                            boxsize = 0.2, 
                            lineheight = "auto",
                            print.pval.Q=FALSE,
                            print.I2 = F,
                            # resid.hetstat = F,
                            print.tau2	= F,
                            common = FALSE,
                            text.random = p_text
               )
               dev.off()
             }
      )
    }
  )
}
