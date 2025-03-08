# Preparation ----
library(shiny)
# library(rsconnect)
# UI
library(shinyWidgets)
library(shinyjs)
library(waiter) # wait while running
library(DT)
# library(shinydashboard)
library(config)

# Manipulate data
library(dplyr)
library(data.table)

# Meta analysis
library(meta)
library(metafor)
library(effsize)

# Plot
library(UpSetR)
library(ggpubr)
library(plotly)
# library(ggrepel)
# library(cowplot)
library(patchwork)

# Multithreads
library(snowfall)
library(parallel)

## Debug
# library(reactlog)

# Load ----
config_list <- config::get(
  config = "test"
  # Default is production mode
)

## Data----
source("Modules/LoadData.R")

## Preprocess ----
source("Modules/Preprocess.R")
# Welcome notification
str1 <- "Nice to meet you."
str2 <- "Very welcome to my version(0.1) —25/01/19"
str3 <- "You can visit https://github.com/mugpeng/DROMA_DB to reach the toturial."
modal_notification <- modalDialog(
  # p("Nice to meet you. \n, test"),
  HTML(paste(str1, str2, str3, sep = '<br/>')),
  title = "Update Notification",
  footer = tagList(
    actionButton("close_modal", "Close")
  )
)

## Function----
source("Package_Function/FuncGetData.R")
source("Package_Function/FuncDrugOmicPair.R")
source("Package_Function/FuncBatchFeature.R")

## Modules----
source("Modules/DrugOmicPair.R")
source("Modules/BatchFeature.R")
source("Modules/StatAnno.R")

# UI ----
ui <- tagList(
  tags$head(
    tags$title("DROMA(Drug Response Omics association MAp)"),
  ),
  autoWaiter(html = spin_loader(), color = transparent(0.5)),
  navbarPage("DROMA(Drug Response Omics association MAp)",
             ## Drugs-omics pairs analysis ----
             tabPanel("Drugs-omics Pairs Analysis",
                      uiDrugOmicPair("DrugOmicPair")
             ),
             ## Features database significant analysis ----
             tabPanel("Batch Features Associations Analysis",
                      uiBatchFeature("BatchFeature")
             ),
             ## Statistics and Annotations ----
             tabPanel("Statistics and Annotations",
                      uiStatAnno("StatAnno")
             ),
             ## Contact ----
             tabPanel("Contact",
                      fluidPage(
                        strong("Feel free to talk with me if you find any bugs or have any suggestions. :)"),
                        p(""),
                        p("Email: mugpeng@outlook.com"),
                        p("Email: yc47680@um.edu.mo"),
                        p("github: https://github.com/mugpeng"),
                        # p("You can visit https://github.com/mugpeng/DROMA_DB to reach the toturial.")
                      ))
  )
)
  
  
# Server ----
server <- function(input, output, session) {
  # Some setup ----
  showModal(modal_notification) # notification
  observeEvent(input$close_modal, {
    removeModal()
  })
  # stop warn
  storeWarn <- getOption("warn")
  options(warn = -1) 
  # Drugs-omics pairs analysis ----
  callModule(serverDrugOmicPair, "DrugOmicPair")
  # Features database significant analysis ----
  callModule(serverBatchFeature, "BatchFeature")
  # Statistics and Annotations ----
  callModule(serverStatAnno, "StatAnno")
}
  
  
# Run ----
shinyApp(ui = ui, server = server)
