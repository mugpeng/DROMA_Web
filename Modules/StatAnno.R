uiStatAnno <- function(id){
  ns <- NS(id)
  fluidPage(
    column(12,
           navlistPanel(
             tabPanel("Overall Drug Information",
                      tabsetPanel(
                        tabPanel("Drug and Sample Counts",
                                 plotOutput(ns("p_count_drugandsample_sum_with_gap2")),
                                 plotOutput(ns("p_count_drugandsample_facet_with_gap"))
                        ),
                        tabPanel("Molecular Characteristics",
                                 plotOutput(ns("p_mol_character"))
                        ),
                        tabPanel("Drug and Sample Overlap Counts",
                                 plotOutput(ns("p_overlap_drug")),
                                 plotOutput(ns("p_overlap_sample"))
                        ),
                        tabPanel("Drug and Sample Annotation visulization",
                                 plotOutput(ns("p_tumor_bubble")),
                                 plotOutput(ns("p_drug_moa"))
                        )
                      )
             ),
             tabPanel("Annotation",
                      tabsetPanel(
                        tabPanel("Sample",
                                 DT::dataTableOutput(ns("sample_anno")),
                                 br(),
                                 wellPanel(
                                   fluidRow(
                                     column(6,
                                            prettyRadioButtons(
                                              inputId = ns("sample_download_type"),
                                              label = "Select download format:",
                                              choices = c(
                                                "Data (RDS)" = "data_rds",
                                                "Data (CSV)" = "data_csv"
                                              ),
                                              selected = "data_csv",
                                              icon = icon("check"), 
                                              animation = "jelly",
                                              status = "primary"
                                            )
                                     ),
                                     column(6,
                                            downloadBttn(
                                              outputId = ns("download_sample"),
                                              label = "Download",
                                              style = "gradient",
                                              color = "default",
                                              block = TRUE,
                                              size = "sm"
                                            )
                                     )
                                   )
                                 )
                        ),
                        tabPanel("Drug",
                                 DT::dataTableOutput(ns("drug_anno")),
                                 br(),
                                 wellPanel(
                                   fluidRow(
                                     column(6,
                                            prettyRadioButtons(
                                              inputId = ns("drug_download_type"),
                                              label = "Select download format:",
                                              choices = c(
                                                "Data (RDS)" = "data_rds",
                                                "Data (CSV)" = "data_csv"
                                              ),
                                              selected = "data_csv",
                                              icon = icon("check"), 
                                              animation = "jelly",
                                              status = "primary"
                                            )
                                     ),
                                     column(6,
                                            downloadBttn(
                                              outputId = ns("download_drug"),
                                              label = "Download",
                                              style = "gradient",
                                              color = "default",
                                              block = TRUE,
                                              size = "sm"
                                            )
                                     )
                                   )
                                 )
                        ),
                      )),
           )),
  )
}

serverStatAnno <- function(input, output, session){
  ns <- session$ns
  # Plot ----
  output$p_count_drugandsample_facet_with_gap <- renderPlot({
    p_count_drugandsample_facet_with_gap
  })
  output$p_count_drugandsample_sum_with_gap2 <- renderPlot({
    p_count_drugandsample_sum_with_gap2
  })
  output$p_mol_character <- renderPlot({
    p_mol_character
  })
  output$p_overlap_drug <- renderPlot({
    p_overlap_drug
  })
  output$p_overlap_sample <- renderPlot({
    p_overlap_sample
  })
  output$p_tumor_bubble <- renderPlot({
    p_tumor_bubble
  })
  output$p_drug_moa <- renderPlot({
    p_drug_moa
  })
  # Get data from DROMA database
  sample_metadata <- reactive({
    getSampleMetadata()
  })

  drug_metadata <- reactive({
    getTreatmentMetadata()
  })

  # Generate statistical plots using DROMA_R
  observe({
    if (.dromadata$initialized) {
      tryCatch({
        # Use DROMA_R to generate statistical plots
        plots <- DROMA.R::generateStatisticalPlots(
          projects = "all",
          connection = .dromadata$connection,
          plot_types = "all"
        )

        # Store plots for output
        p_count_drugandsample_facet_with_gap <<- plots$sample_drug_counts_facet
        p_count_drugandsample_sum_with_gap2 <<- plots$sample_drug_counts_summary
        p_mol_character <<- plots$molecular_characteristics
        p_overlap_drug <<- plots$drug_overlap
        p_overlap_sample <<- plots$sample_overlap
        p_tumor_bubble <<- plots$tumor_type_bubble
        p_drug_moa <<- plots$drug_moa

      }, error = function(e) {
        # Create placeholder plots if DROMA_R function fails
        library(ggplot2)
        empty_plot <- ggplot() + theme_void() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "Statistical plots will be available\nafter full DROMA integration")

        p_count_drugandsample_facet_with_gap <<- empty_plot
        p_count_drugandsample_sum_with_gap2 <<- empty_plot
        p_mol_character <<- empty_plot
        p_overlap_drug <<- empty_plot
        p_overlap_sample <<- empty_plot
        p_tumor_bubble <<- empty_plot
        p_drug_moa <<- empty_plot
      })
    }
  })

  # Table outputs using DROMA data
  output$drug_anno <- DT::renderDataTable({
    DT::datatable(
      drug_metadata(),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left; color: black; font-size: 14px;',
        htmltools::strong("Drug Annotation Data")
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv')
      ),
      extensions = 'Buttons',
      rownames = FALSE,
      filter = 'top'
    )
  })

  output$sample_anno <- DT::renderDataTable({
    DT::datatable(
      sample_metadata(),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left; color: black; font-size: 14px;',
        htmltools::strong("Sample Annotation Data")
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv')
      ),
      extensions = 'Buttons',
      rownames = FALSE,
      filter = 'top'
    )
  })
  
  # Download handlers for sample annotations
  output$download_sample <- downloadHandler(
    filename = function() {
      type <- input$sample_download_type
      base_name <- "DROMA_Sample_Annotations"
      
      switch(type,
             "data_rds" = paste0(base_name, ".rds"),
             "data_csv" = paste0(base_name, ".csv")
      )
    },
    content = function(filename) {
      type <- input$sample_download_type
      
      switch(type,
             "data_rds" = saveRDS(sample_metadata(), filename),
             "data_csv" = write.csv(sample_metadata(), filename, row.names = FALSE)
      )
    }
  )
  
  # Download handlers for drug annotations
  output$download_drug <- downloadHandler(
    filename = function() {
      type <- input$drug_download_type
      base_name <- "DROMA_Drug_Annotations"
      
      switch(type,
             "data_rds" = paste0(base_name, ".rds"),
             "data_csv" = paste0(base_name, ".csv")
      )
    },
    content = function(filename) {
      type <- input$drug_download_type
      
      switch(type,
             "data_rds" = saveRDS(drug_metadata(), filename),
             "data_csv" = write.csv(drug_metadata(), filename, row.names = FALSE)
      )
    }
  )
}
