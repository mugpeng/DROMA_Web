# UI component for Global Settings
uiGlobalSetting <- function(id) {
  ns <- NS(id)
  tagList(
    # Custom CSS for floating widget
    tags$head(
      tags$style(HTML(sprintf("
        #%s {
          position: fixed;
          left: -100px;  /* Hide most of the button by default */
          top: 50%%;
          transform: translateY(-50%%);
          background-color: rgba(200, 200, 200, 0.3);
          padding: 0;
          border-radius: 5px;
          border: 1px solid #ddd;
          cursor: pointer;
          transition: all 0.3s ease;
          z-index: 1000;
          box-shadow: 0 2px 5px rgba(0,0,0,0.1);
          opacity: 0.3;
          display: flex;
          align-items: center;
          justify-content: center;
          width: 150px;  /* Fixed width for smooth animation */
        }
        #%s:hover {
          left: 20px;  /* Slide out to show full button */
          opacity: 1;
          box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }
        #%s .btn {
          width: 100%%;
          height: 100%%;
          padding: 10px 15px;
          margin: 0;
          border: none;
          background: none;
          color: #333;
          cursor: pointer;
          text-align: center;
          white-space: nowrap;
        }
        #%s .btn:hover,
        #%s .btn:focus,
        #%s .btn:active {
          background: none !important;
          box-shadow: none !important;
          outline: none !important;
        }
      ", ns("global_settings_btn"), ns("global_settings_btn"), 
         ns("global_settings_btn"), ns("global_settings_btn"),
         ns("global_settings_btn"), ns("global_settings_btn")))
    )),
    # Floating Global Settings Button
    div(
      id = ns("global_settings_btn"),
      actionButton(
        ns("global_settings"),
        "Global Settings",
        class = "btn-block"
      )
    )
  )
}

# Server component for Global Settings
serverGlobalSetting <- function(input, output, session) {
  ns <- session$ns
  
  # Initialize reactive values for storing states
  rv <- reactiveValues(
    zscore_enabled = TRUE,  # Set default to TRUE
    settings_changed = FALSE,  # Track if settings have changed
    last_update = Sys.time()  # Timestamp for last update
  )
  
  # Apply z-score on initialization
  observe({
    # Only run once on initialization
    if (!exists("normalization_state", envir = .GlobalEnv) || 
        !isTRUE(base::get("normalization_state", envir = .GlobalEnv))) {
      # Source z-score module
      if (!exists("apply_zscore_normalization", envir = .GlobalEnv)) {
        source("Package_Function/FuncZscoreWhole.R", local = FALSE)
      }
      # Apply z-score normalization
      apply_zscore_normalization()
      
      # Set global z-score state for other modules to check
      assign("GLOBAL_ZSCORE_STATE", list(
        enabled = TRUE,
        timestamp = Sys.time()
      ), envir = .GlobalEnv)
    }
  }, priority = 1000)
  
  # Global Settings Modal
  observeEvent(input$global_settings, {
    showModal(modalDialog(
      title = "Global Settings",
      div(
        checkboxInput(ns("zscore_checkbox"), "Apply Z-score normalization to data", value = rv$zscore_enabled),
        tags$div(
          style = "margin-left: 20px; color: #666; font-size: 0.9em;",
          tags$p("When enabled:"),
          tags$ul(
            tags$li("Omics data (mRNA, CNV, methylation, protein): Z-score normalization is applied across samples for each gene"),
            tags$li("Drug data: Z-score normalization is applied row-by-row, where each drug is independently normalized across all cell lines by subtracting its mean and dividing by its standard deviation"),
            tags$li("Drug-Omic Pair Analysis: A merged dataset plot will be created, combining data from all datasets for a comprehensive view (this merged dataset is excluded from meta-analysis)")
          ),
          tags$p("This helps standardize data and make comparisons more meaningful.")
        )
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton(ns("save_settings"), "Save Settings")
      ),
      size = "m",
      easyClose = TRUE
    ))
  })
  
  # Handle save settings
  observeEvent(input$save_settings, {
    # Check if z-score setting has changed
    if (rv$zscore_enabled != input$zscore_checkbox) {
      rv$zscore_enabled <- input$zscore_checkbox
      rv$settings_changed <- TRUE
      
      # Show loading message
      showModal(modalDialog(
        "Reloading data with new normalization settings...",
        footer = NULL
      ))
      
      # Source z-score module if needed
      if (!exists("apply_zscore_normalization", envir = .GlobalEnv) || 
          !exists("reset_to_original", envir = .GlobalEnv)) {
        source("Package_Function/FuncZscoreWhole.R", local = FALSE)
      }
      
      # Apply appropriate normalization
      if (rv$zscore_enabled) {
        apply_zscore_normalization()
      } else {
        reset_to_original()
      }
      
      # Update timestamp to trigger reactivity in other modules
      rv$last_update <- Sys.time()
      
      # Store the normalization state in a global variable for other modules to check
      assign("GLOBAL_ZSCORE_STATE", list(
        enabled = rv$zscore_enabled,
        timestamp = rv$last_update
      ), envir = .GlobalEnv)
      
      # Remove loading message
      removeModal()
      
      # Show confirmation message with instructions
      showNotification(
        HTML(paste(
          "Data has been", 
          if(rv$zscore_enabled) "normalized" else "reset to original values"
          # "<br><b>Please re-run your analysis to see the changes.</b>"
        )),
        type = "message",
        duration = 10
      )
    }
    removeModal()
  })
  
  # Return reactive values to allow access to state from main app
  return(rv)
} 