uiStatAnno <- function(id){
  ns <- NS(id)
  fluidPage(
    column(12,
           navlistPanel(
             tabPanel("Overall Drug Information",
                      tabsetPanel(
                        tabPanel("Drug and Sample Counts",
                                 plotOutput(ns("p_count_combined"))
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
                                 DT::dataTableOutput(ns("sample_anno"))),
                        tabPanel("Drug",
                                 DT::dataTableOutput(ns("drug_anno"))),
                      )),
           )),
  )
}

serverStatAnno <- function(input, output, session){
  ns <- session$ns
  # Plot ----
  output$p_count_combined <- renderPlot({
    p_count_combined
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
    drug_anno
  }, options = list(scrollX = TRUE), selection = 'single')
  output$sample_anno <- DT::renderDataTable({ 
    sample_anno
  }, options = list(scrollX = TRUE), selection = 'single')
}
