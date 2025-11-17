#' ==========================================
#' Common Utility Functions
#' ==========================================
#' These are shared utility functions used across modules

#' Format time duration into human-readable string
#'
#' @param seconds Time in seconds
#' @return Formatted time string
#'
#' @examples
#' format_time(65)  # Returns "1 minute 5 seconds"
#' format_time(3661)  # Returns "1 hour 1 minute 1 second"
format_time <- function(seconds) {
  if (seconds < 60) {
    return(paste0(round(seconds, 1), " seconds"))
  } else if (seconds < 3600) {
    mins <- floor(seconds / 60)
    secs <- round(seconds %% 60)
    if (secs == 0) {
      return(paste0(mins, " minute", if(mins != 1) "s" else ""))
    } else {
      return(paste0(mins, " minute", if(mins != 1) "s" else "", " ", secs, " seconds"))
    }
  } else {
    hours <- floor(seconds / 3600)
    mins <- floor((seconds %% 3600) / 60)
    secs <- round(seconds %% 60)

    parts <- character(0)
    if (hours > 0) {
      parts <- c(parts, paste0(hours, " hour", if(hours != 1) "s" else ""))
    }
    if (mins > 0) {
      parts <- c(parts, paste0(mins, " minute", if(mins != 1) "s" else ""))
    }
    if (secs > 0 && length(parts) < 2) {
      parts <- c(parts, paste0(secs, " seconds"))
    }

    return(paste(parts, collapse = " "))
  }
}

#' Format p-values for display
#'
#' @param p Numeric p-values
#' @param digits Number of decimal places
#' @return Formatted p-values as strings
#'
#' @examples
#' format_pvalue(c(0.0001, 0.045, 0.5))
format_pvalue <- function(p, digits = 3) {
  sapply(p, function(x) {
    if (x < 0.001) {
      return("< 0.001")
    } else if (x < 0.01) {
      return(paste0("=", round(x, digits)))
    } else if (x < 0.05) {
      return(paste0("=", round(x, 2)))
    } else {
      return(paste0("=", round(x, digits)))
    }
  })
}

#' Create a safe data frame from list input
#'
#' @param data_list List containing data
#' @param default_col_names Default column names if needed
#' @return A data frame
#'
#' @examples
#' df <- safe_dataframe(list(data))
safe_dataframe <- function(data_list, default_col_names = c("ID", "Value")) {
  if (is.null(data_list)) {
    return(data.frame())
  }

  if (is.data.frame(data_list)) {
    return(data_list)
  }

  if (is.list(data_list) && !is.null(data_list$data)) {
    if (is.data.frame(data_list$data)) {
      return(data_list$data)
    }
  }

  # Try to convert to data frame
  tryCatch({
    return(as.data.frame(data_list))
  }, error = function(e) {
    # Return empty data frame with default columns
    df <- data.frame(matrix(NA, nrow = 0, ncol = length(default_col_names)))
    colnames(df) <- default_col_names
    return(df)
  })
}

#' Validate input parameters
#'
#' @param param Parameter value
#' @param param_name Name of the parameter
#' @param allowed_values Allowed values (optional)
#' @param class Expected class (optional)
#' @return TRUE if valid, FALSE otherwise
#'
#' @examples
#' validate_input("mRNA", "omic_type", allowed_values = c("mRNA", "cnv", "meth"))
validate_input <- function(param, param_name, allowed_values = NULL, class = NULL) {
  if (is.null(param)) {
    message(param_name, " cannot be NULL")
    return(FALSE)
  }

  if (!is.null(class) && !inherits(param, class)) {
    message(param_name, " must be of class ", class)
    return(FALSE)
  }

  if (!is.null(allowed_values) && !param %in% allowed_values) {
    message(param_name, " must be one of: ", paste(allowed_values, collapse = ", "))
    return(FALSE)
  }

  return(TRUE)
}

#' Show loading notification
#'
#' @param session Shiny session object
#' @param message Loading message
#' @param duration Duration in seconds (NULL for manual removal)
#' @return Notification ID
#'
#' @examples
#' show_loading(session, "Processing data...")
show_loading <- function(session, message = "Loading...", duration = NULL) {
  id <- paste0("loading_", sample(10000:99999, 1))

  if (is.null(duration)) {
    showNotification(
      div(id = id,
          icon("spinner fa-pulse"),
          span(message),
          style = "padding: 10px; background-color: #f8f9fa; border-left: 4px solid #007bff;"
      ),
      type = "default",
      duration = NULL,
      closeButton = FALSE
    )
  } else {
    showNotification(
      div(
        icon("spinner fa-pulse"),
        span(message)
      ),
      type = "default",
      duration = duration * 1000
    )
  }

  return(id)
}

#' Hide loading notification
#'
#' @param id Notification ID returned by show_loading
#'
#' @examples
#' hide_loading("loading_12345")
hide_loading <- function(id) {
  # Shiny doesn't have a direct way to hide notifications
  # This is a placeholder for future implementation
  # Could use custom HTML/JS for more control
}

#' Get file extension
#'
#' @param filename File name or path
#' @return File extension without the dot
#'
#' @examples
#' get_file_extension("plot.pdf")  # Returns "pdf"
get_file_extension <- function(filename) {
  ext <- tools::file_ext(filename)
  return(ext)
}

#' Generate unique ID
#'
#' @param prefix Prefix for the ID
#' @return Unique ID string
#'
#' @examples
#' generate_id("plot")  # Returns something like "plot_1704123456789"
generate_id <- function(prefix = "") {
  timestamp <- as.numeric(Sys.time()) * 1000
  random <- sample(1000:9999, 1)
  id <- paste0(prefix, "_", timestamp, "_", random)
  return(id)
}

#' Log activity to console or file
#'
#' @param message Message to log
#' @param level Log level (INFO, WARNING, ERROR)
#' @param log_file Optional log file path
#'
#' @examples
#' log_message("Analysis completed", "INFO")
log_message <- function(message, level = "INFO", log_file = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste0("[", timestamp, "] ", level, ": ", message)

  # Print to console
  message(log_entry)

  # Optionally write to file
  if (!is.null(log_file)) {
    write(log_entry, file = log_file, append = TRUE)
  }
}