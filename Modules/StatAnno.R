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
  # Table ----
  output$drug_anno <- DT::renderDataTable({ 
    DT::datatable(
      drug_anno,
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
      sample_anno,
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
             "data_rds" = saveRDS(sample_anno, filename),
             "data_csv" = write.csv(sample_anno, filename, row.names = FALSE)
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
             "data_rds" = saveRDS(drug_anno, filename),
             "data_csv" = write.csv(drug_anno, filename, row.names = FALSE)
      )
    }
  )
}
