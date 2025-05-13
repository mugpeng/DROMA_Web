
# DrugFeature ----
# Create a comparison plot for continuous variables (like Age)
plot_continuous_comparison <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity") {
  # Create scatter plot with correlation information
  p <- ggscatter(data, x = cont_column, y = value_column, alpha = 0.2) +
    stat_cor(size = 6, method = "spearman") + 
    stat_smooth(formula = y ~ x, method = "lm") + 
    theme_bw() +
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12)
    ) + 
    ggtitle(paste(value_label, "vs", cont_column)) 
    
  
  return(p)
}

# Create boxplots for continuous variable groups
plot_continuous_groups <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity", num_bins = 4) {
  # Create bins for the continuous variable
  cont_values <- data[[cont_column]]
  
  # Create bins
  cont_bins <- cut(cont_values, 
                   breaks = num_bins,
                   include.lowest = TRUE,
                   labels = FALSE)
  
  # Create labels for groups
  cont_range <- range(cont_values, na.rm = TRUE)
  bin_width <- diff(cont_range) / num_bins
  group_labels <- sapply(1:num_bins, function(i) {
    lower <- cont_range[1] + (i-1) * bin_width
    upper <- cont_range[1] + i * bin_width
    paste0(round(lower), "-", round(upper))
  })
  
  # Add group information to data
  group_column <- paste0(cont_column, "_group")
  label_column <- paste0(cont_column, "_group_label")
  
  data[[group_column]] <- cont_bins
  data[[label_column]] <- group_labels[data[[group_column]]]
  
  # Create boxplot with statistical test
  p <- ggboxplot(data, x = label_column, y = value_column,
                 fill = label_column, palette = "jco",
                 add = "jitter", add.params = list(alpha = 0.15)) +
    stat_compare_means(size = 6, label.x = 0.8,
                       label.y = (max(data[[value_column]]) - max(data[[value_column]])/8),
                       label = "p.format") + 
    theme_bw() + 
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggtitle(paste(value_label, "by", cont_column, "Group")) 
    
  
  return(p)
}

# Create a categorical comparison plot
plot_category_comparison <- function(data, category_column, value_column = "value", value_label = "Drug Sensitivity") {
  # Count observations per category and filter out categories with too few samples
  category_counts <- table(data[[category_column]])
  valid_categories <- names(category_counts)[category_counts >= 3]
  
  if (length(valid_categories) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "Not enough samples per category for comparison") + 
             theme_void())
  }
  
  data_filtered <- data[data[[category_column]] %in% valid_categories, ]
  
  # Create improved boxplot with consistent styling
  p <- ggboxplot(data_filtered, x = category_column, y = value_column,
                 fill = category_column, 
                 palette = bright_palette_26,
                 add = "jitter", 
                 add.params = list(alpha = 0.15)) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggtitle(paste(value_label, "by", category_column))
    
  
  # Add statistical comparison with appropriate method
  if (length(valid_categories) >= 2) {
    max_y <- max(data_filtered[[value_column]], na.rm = TRUE)
    label_x_pos <- length(valid_categories) / 2.5
    if (length(valid_categories) == 2) {
      # For two groups, use wilcoxon with clear label positioning
      p <- p + stat_compare_means(size = 6, 
                                  label.x = label_x_pos,
                                  label.y = (max_y - max_y/8),
                                  label = "p.format")
    } else {
      # For more than two groups, use Kruskal-Wallis
      # Add global p-value at top
      p <- p + stat_compare_means(method = "kruskal.test", 
                                  size = 6,
                                  label.x = label_x_pos,
                                  label.y = max_y + (max_y * 0.15),
                                  label = "p.format")
    }
  }
  
  return(p)
}

# Create a drug feature comparison plot
create_drug_comparison_plot <- function(data, comparison_var, value_column = "value", value_label = "Drug Sensitivity", 
                                        num_bins = 4, show_groups_boxplot = TRUE) {
  # Handle missing values in the comparison variable
  data <- data[!is.na(data[[comparison_var]]), ]
  
  if (nrow(data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No data available for this comparison") + 
             theme_void())
  }
  
  # Check if the comparison variable is numeric/continuous
  if (is.numeric(data[[comparison_var]])) {
    # For continuous variables
    p1 <- plot_continuous_comparison(data, cont_column = comparison_var, 
                                     value_column = value_column, value_label = value_label)
    
    # Also create a boxplot with bins if requested
    if (show_groups_boxplot) {
      # Create grouped boxplot
      p2 <- plot_continuous_groups(data, cont_column = comparison_var, 
                                   value_column = value_column, value_label = value_label,
                                   num_bins = num_bins)
      
      # Return grid of both plots
      return(grid.arrange(p1, p2, ncol = 2))
    }
    
    return(p1)
  } else {
    # For categorical variables
    return(plot_category_comparison(data, category_column = comparison_var, 
                                    value_column = value_column, value_label = value_label))
  }
}

# DrugOmicPair ----
create_plot_with_common_axes <- function(p, x_title = "Common X-Axis Title", 
                                         y_title = "Common Y-Axis Title") {
  # Create a function that will generate the plot when called
  function() {
    # Convert patchwork to a grob
    p_grob <- patchworkGrob(p)
    
    # Create a new plotting area
    grid.newpage()
    
    # Draw the patchwork
    grid.draw(p_grob)
    
    # Add common x-axis title
    grid.text(x_title,
              x = 0.5, y = 0.02,
              gp = gpar(fontsize = 18, fontface = "bold"))
    
    # Add common y-axis title (rotated)
    grid.text(y_title,
              x = 0.01, y = 0.5,
              rot = 90,
              gp = gpar(fontsize = 18, fontface = "bold"))
    
    # Return the grob for potential further use
    invisible(p_grob)
  }
}

# Create a standardized forest plot for meta-analysis results
create_forest_plot <- function(meta_obj, 
                              xlab = "Effect Size (95% CI)", 
                              show_common = FALSE) {
  # Validate input
  if (!inherits(meta_obj, "meta") && !inherits(meta_obj, "metagen")) {
    stop("Input must be a meta-analysis object from the 'meta' package")
  }
  
  # Format p-value text for random effects model
  p_val <- meta_obj$pval.random
  p_text <- if(p_val < 0.001) {
    paste("Random-Effects Model (p =", format(p_val, scientific = TRUE, digits = 3), ")")
  } else {
    paste("Random-Effects Model (p =", round(p_val, 3), ")")
  }
  
  # Create forest plot
  meta::forest(meta_obj, 
               xlab = xlab, 
               slab = "study", 
               print.pval.common = show_common,
               boxsize = 0.2, 
               lineheight = "auto",
               print.pval.Q = FALSE,
               print.I2 = FALSE,
               print.tau2 = FALSE,
               common = show_common,
               text.random = p_text
  )
}

# Plot continuous drug-omic correlation for a single study
plot_continuous_drugomic <- function(omic_values, drug_values, study_name) {
  # Combine data into dataframe
  cor_df <- data.frame(
    genes = omic_values,
    drugs = drug_values
  )
  
  # Create scatter plot with correlation statistics
  ggscatter(cor_df, x = "genes", y = "drugs", alpha = 0.2) +
    stat_cor(size = 6, method = "spearman") + 
    stat_smooth(formula = y ~ x, method = "lm") + 
    theme_bw() +
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12)
    ) + 
    ggtitle(study_name)
}

# Create plots for all continuous drug-omic pairs
plot_all_continuous_drugomic <- function(pairs_list) {
  # Initialize list to store plots
  p_list <- list()
  
  # Create plot for each pair
  for (i in seq_along(pairs_list)) {
    
    # Try to create the plot, continue if error
    tryCatch({
      omic_sel <- pairs_list[[i]]$omic
      drug_sel <- pairs_list[[i]]$drug
      
      # Ensure adequate data for plotting
      if (length(omic_sel) < 3 || length(drug_sel) < 3) next
      
      # Create plot and add to list
      p_list[[i]] <- plot_continuous_drugomic(omic_sel, drug_sel, names(pairs_list)[i])
    }, error = function(e) {
      # Continue to next pair on error
    })
  }
  
  # Remove NULL entries from list
  p_list <- p_list[!sapply(p_list, is.null)]
  
  # Combine plots using patchwork if plots exist
  if (length(p_list) > 0) {
    return(wrap_plots(p_list, ncol = 3))
  } else {
    return(NULL)
  }
}

# Plot discrete drug-omic comparison for a single study
plot_discrete_drugomic <- function(yes_values, no_values, study_name) {
  # Combine data into dataframe
  box_df <- data.frame(
    drugs = c(no_values, yes_values),
    events = rep(c("no", "yes"), times = c(length(no_values), length(yes_values)))
  )
  
  # Create boxplot with statistical test
  ggboxplot(data = box_df, x = "events", y = "drugs",
            fill = "events", palette = c("#BEBADAFF", "#FB8072FF"),
            add = "jitter", add.params = list(alpha = 0.15)) + 
    stat_compare_means(size = 6, label.x = 0.8,
                       label.y = (max(box_df$drugs) - max(box_df$drugs)/8),
                       label = "p.format") + 
    theme_bw() + 
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "none"
    ) + 
    coord_cartesian(ylim = c(NA, max(box_df$drugs) + max(box_df$drugs)/20)) + 
    ggtitle(study_name)
}

# Create plots for all discrete drug-omic pairs
plot_all_discrete_drugomic <- function(pairs_list) {
  # Initialize list to store plots
  p_list <- list()
  
  # Create plot for each pair
  for (i in seq_along(pairs_list)) {
    
    # Try to create the plot, continue if error
    tryCatch({
      yes_drugs <- pairs_list[[i]]$yes
      no_drugs <- pairs_list[[i]]$no
      
      # Ensure adequate data for plotting
      if (length(yes_drugs) < 3 || length(no_drugs) < 3) next
      
      # Create plot and add to list
      p_list[[i]] <- plot_discrete_drugomic(yes_drugs, no_drugs, names(pairs_list)[i])
    }, error = function(e) {
      # Continue to next pair on error
    })
  }
  
  # Remove NULL entries from list
  p_list <- p_list[!sapply(p_list, is.null)]
  
  # Combine plots using patchwork if plots exist
  if (length(p_list) > 0) {
    return(wrap_plots(p_list, ncol = 3))
  } else {
    return(NULL)
  }
}

# BatchFeature ----
#' Create a volcano plot from meta-analysis results
#' 
#' @param meta_df Data frame containing meta-analysis results with columns: 
#'                effect_size, p_value, and name
#' @param es_t Effect size threshold to consider significant
#' @param P_t P-value threshold to consider significant
#' @param label Whether to add labels to top points (TRUE/FALSE)
#' @param top_label_each Number of top points in each direction to label
#' @param label_size Size of text labels
#' @param point_size Size of points
#' @param point_alpha Alpha transparency of points
#' @param title Plot title (NULL for no title)
#' @param p_adj_method Method for p-value adjustment ("none", "BH", "bonferroni")
#' @param custom_colors Custom color vector for Up, NS, Down (NULL for defaults)
#' @return ggplot object with volcano plot
plotMetaVolcano <- function(meta_df, 
                            es_t = .4, 
                            P_t = .001,
                            label = TRUE,
                            top_label_each = 5,
                            label_size = 5,
                            point_size = 2.5,
                            point_alpha = 0.6,
                            title = NULL,
                            p_adj_method = "none",
                            custom_colors = NULL) {
  
  # Input validation
  if(!is.data.frame(meta_df)) stop("meta_df must be a data frame")
  if(!all(c("effect_size", "p_value", "name") %in% colnames(meta_df))) {
    stop("meta_df must contain columns: effect_size, p_value, and name")
  }
  
  # Handle p-value adjustment if requested
  if(p_adj_method != "none") {
    meta_df$p_value <- p.adjust(meta_df$p_value, method = p_adj_method)
  }
  
  # Default colors
  if(is.null(custom_colors)) {
    custom_colors <- c("Down" = "#44bce4", "NS" = "grey", "Up" = "#fc7474")
  }
  
  # Group the points based on thresholds
  meta_df$group <- dplyr::case_when(
    meta_df$effect_size > es_t & meta_df$p_value < P_t ~ "Up",
    meta_df$effect_size < -es_t & meta_df$p_value < P_t ~ "Down",
    TRUE ~ "NS"
  )  
  
  # Count significant findings
  sig_counts <- table(meta_df$group)
  sig_text <- paste0(
    "Up: ", sum(meta_df$group == "Up"), ", ",
    "Down: ", sum(meta_df$group == "Down"), ", ",
    "Total: ", nrow(meta_df)
  )
  
  # Basic volcano plot
  p <- ggplot(data = meta_df, 
              aes(x = effect_size, 
                  y = -log10(p_value))) +
    geom_point(size = point_size, alpha = point_alpha, 
               aes(color = group)) +
    theme_bw() + 
    theme(
      legend.position = "none",
      title = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 15, colour = "black"), 
      axis.text = element_text(size = 15, color = "black"), 
      legend.title = element_text(size = 15, colour = "black"),
      legend.text = element_text(size = 15),
      text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black")
    ) + 
    ylab("-log10(Pvalue)") + 
    xlab("Effect Size") +
    scale_color_manual(values = custom_colors) + 
    geom_vline(xintercept = c(-es_t, es_t), lty = 4, col = "black", lwd = 0.5) + 
    geom_hline(yintercept = -log10(P_t), lty = 4, col = "black", lwd = 0.5) +
    annotate("text", x = min(meta_df$effect_size, na.rm = TRUE) * 0.8, 
             y = max(-log10(meta_df$p_value), na.rm = TRUE) * 0.9, 
             label = sig_text, hjust = 0, size = 5)
  
  # Add title if provided
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  # Add labels if requested
  if(label) {
    meta_df2 <- meta_df[meta_df$group != "NS",]
    
    # Skip labeling if there are no significant points
    if(nrow(meta_df2) > 0) {
      # Get top points to label
      low_indices <- head(order(meta_df2$effect_size), min(top_label_each, nrow(meta_df2)))
      high_indices <- tail(order(meta_df2$effect_size), min(top_label_each, nrow(meta_df2)))
      forlabel_names <- c(meta_df2$name[low_indices], meta_df2$name[high_indices])
      forlabel_df <- meta_df2[meta_df2$name %in% forlabel_names,]
      
      p <- p + 
        geom_point(size = point_size + 0.5, shape = 1, data = forlabel_df) +
        ggrepel::geom_text_repel(
          data = forlabel_df,
          aes(label = name),
          size = label_size,
          color = "black",
          box.padding = 0.5,
          point.padding = 0.3,
          force = 5,
          max.overlaps = 20
        )
    }
  }
  
  p
}

