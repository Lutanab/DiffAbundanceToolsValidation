library(ggplot2)
library(ggrepel)

plot_differential_abundance2 <- function(plot_data, thresholds, number_of_samples, da_tool, pval_cap = 100, lfc_threshold = 1.5, pval_threshold = 0.01, show_theor_lfc = T) {
  
  # Adjust size levels based on highlight status
  plot_data$size <- ifelse(plot_data$highlight == "Truly Differentially Abundant", plot_data$taxa_prop, plot_data$taxa_prop / 10)
  plot_data$alpha <- ifelse(plot_data$highlight == "Truly Differentially Abundant", plot_data$taxa_prop, pmin(plot_data$taxa_prop, 0.2))
  
  # Cap the -log10(padj) values
  plot_data$capped_log_padj <- pmin(-log10(plot_data$padj), pval_cap)
  
  # Calculate the difference between blue and red points
  plot_data$difference <- abs(plot_data$log2FoldChange - plot_data$theor_lfc)
  
  median_delta <- median(plot_data$difference[plot_data$highlight == "Truly Differentially Abundant"], na.rm = TRUE)
  
  # Filter data for truly differentially abundant points
  true_da_data <- plot_data[plot_data$highlight == "Truly Differentially Abundant", ]
  lfc_threshold_list <- c(-lfc_threshold, lfc_threshold)
  
  # Filter data for correct and incorrect signals
  correct_signals <- plot_data[plot_data$is_correct == TRUE, ]
  incorrect_signals <- plot_data[plot_data$not_correct == TRUE, ]
  
  
  # Create the ggplot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = capped_log_padj, size = size, alpha = alpha)) +
    # Fill the quadrants with semi-transparent rectangles
    geom_rect(aes(xmin = -Inf, xmax = -lfc_threshold, ymin = -log10(pval_threshold), ymax = Inf), fill = "lightblue") +
    geom_rect(aes(xmin = lfc_threshold, xmax = Inf, ymin = -log10(pval_threshold), ymax = Inf), fill = "lightblue") +
    # geom_rect(aes(xmin = -lfc_threshold, xmax = lfc_threshold, ymin = -Inf, ymax = -log10(pval_threshold)), fill = "lightcoral", alpha = 0.3) +
    
    geom_point(aes(color = highlight)) +
    # Add green outline for selected points
    geom_point(data = correct_signals, aes(x = log2FoldChange, y = capped_log_padj), color = "green", size = 3, shape = 2, stroke = 1.5, alpha = 1) +
    # Add red crosses on incorrect signals
    # geom_point(data = incorrect_signals, aes(x = log2FoldChange, y = capped_log_padj), color = "orange", size = 5, shape = 4, stroke = 1.5, alpha=1) +
    # geom_point(data = true_da_data, aes(x = theor_lfc, y = capped_log_padj), color = "blue", size = 3, shape = 21, fill = "blue") +
    # geom_text(data = true_da_data, aes(x = (log2FoldChange + theor_lfc) / 2, y = capped_log_padj, 
    #                                    label = round(difference, 2)), color = "blue", vjust = -1.5, size = 3, check_overlap = TRUE) + 
    # geom_segment(data = true_da_data, aes(x = log2FoldChange, y = capped_log_padj, xend = theor_lfc, yend = capped_log_padj), color = "blue", size = 0.5) +
    geom_point(data = incorrect_signals, aes(x = log2FoldChange, y = capped_log_padj), color = "orange", size = 5, shape = 4, stroke = 1.5, alpha=1) +
    
    scale_color_manual(values = c("Truly Differentially Abundant" = "red", "Not Differentially Abundant" = "black")) +
    scale_alpha_continuous(range = c(0.3, 1), guide = "none") +  # Control opacity based on calculated alpha
    scale_size_continuous(name = "Taxa Prevalence") +  # Adjust point size range and add size legend
    scale_y_continuous(labels = scales::comma, limits = c(-0.5, pval_cap)) +  # Add space below zero on the y-axis and limit
    labs(
      title = paste0("log2FoldChange vs Adjusted P-Values, ", da_tool, " (", number_of_samples, " samples)"),
      x = "ln Fold Change",
      y = "-log10 Adjusted P-Value",
      caption = "Points size proportional to taxa proportion"
    ) +
    theme_minimal(base_size = 15) +  # Use a larger base font size for readability
    theme(
      legend.position = "right",  # Position the legend on the right
      plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold the plot title
      axis.title = element_text(face = "bold"),  # Bold the axis titles
      panel.grid.major = element_line(color = "grey80"),  # Add major grid lines
      panel.grid.minor = element_blank()  # Remove minor grid lines
    ) +
    geom_hline(yintercept = -log10(thresholds), linetype = "dashed", color = "darkgray") +  # Add significance threshold lines
    annotate("text", x = max(plot_data$log2FoldChange, na.rm = TRUE), y = -log10(thresholds), 
             label = paste0("p = ", thresholds), hjust = 1, vjust = -0.5, color = "darkgray", size = 5) +
    # Annotate the median delta
    annotate("text", x = Inf, y = Inf, label = paste("Median Delta:", round(median_delta, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "blue", fontface = "bold") +
    
    geom_vline(xintercept = lfc_threshold_list, linetype = "dashed", color = "deeppink") +  # Add vertical threshold line for log2FoldChange
    geom_vline(xintercept = lfc_threshold_list, linetype = "dashed", color = "deeppink") +   # Add vertical threshold line for -log2FoldChange
    annotate("text", x = lfc_threshold, y = max(plot_data$capped_log_padj, na.rm = TRUE),
             label = paste0("LFC = ", lfc_threshold), hjust = 1.1, vjust = 1.5, color = "deeppink", size = 5, angle = 90) +
    annotate("text", x = -lfc_threshold, y = max(plot_data$capped_log_padj, na.rm = TRUE),
             label = paste0("LFC = -", lfc_threshold), hjust = 1.1, vjust = 1.5, color = "deeppink", size = 5, angle = 90) + 
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "deeppink") + 
    annotate("text", x = max(plot_data$log2FoldChange, na.rm = TRUE), y = -log10(pval_threshold), 
             label = paste0("p = ", pval_threshold), hjust = 1, vjust = -0.5, color = "deeppink", size = 5)
  
  # Add the custom legend entries
  p <- p + 
    scale_shape_manual(
      values = c(16, 24, 15, 21), 
      name = "highlight", 
      labels = c("Not Differentially Abundant", 
                 "Truly Differentially Abundant",
                 "Blue shaded areas: non-significant differentially abundant\npoints (log2 fold change < 1.5 and > -1.5)", 
                 "Blue lines: Difference of modeled and observed LFC")
    ) +
    scale_color_manual(
      values = c("black", "red", "blue", "blue"), 
      name = "highlight", 
      labels = c("Not Differentially Abundant", 
                 "Truly Differentially Abundant",
                 "Blue shaded areas: non-significant differentially abundant\npoints (log2 fold change < 1.5 and > -1.5)", 
                 "Blue lines: Difference of modeled and observed LFC")
    ) +
    scale_size_continuous(name = "Taxa Prevalence", range = c(1, 10)) +
    guides(
      color = guide_legend(override.aes = list(size = 4)), 
      shape = guide_legend(override.aes = list(size = 4))
    ) +
    theme(legend.position = "right") +
    labs(color = "highlight", shape = "highlight") +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    )
  
  return(p)
}




diff_abundance_3d <- function(result_data){
  
  # Define RGB colors
  color_truly_diff <- "255, 0, 0"     # Red for "Truly Differentially Abundant"
  color_not_diff <- "0, 0, 0"         # Black for "Not Differentially Abundant"
  
  # Create RGBA colors based on highlight and opacity
  result_data$rgba_color <- ifelse(
    result_data$highlight == "Truly Differentially Abundant",
    paste0("rgba(", color_truly_diff, ", ", 1, ")"),
    paste0("rgba(", color_not_diff, ", ", 0.3, ")")
  )
  
  # Define the truly differentially abundant data
  true_da_data <- result_data[result_data$highlight == "Truly Differentially Abundant", ]
  true_da_data$pair_id <- seq_len(nrow(true_da_data))
  
  # Create separate dataframes for red and blue points
  red_points <- true_da_data[, c("log2FoldChange", "taxa_prop", "log10padj")]
  colnames(red_points) <- c("X", "Y", "Z")
  red_points$color <- "red"
  red_points$pair_id <- true_da_data$pair_id
  
  blue_points <- true_da_data[, c("theor_lfc", "taxa_prop", "log10padj")]
  colnames(blue_points) <- c("X", "Y", "Z")
  blue_points$color <- "blue"
  blue_points$pair_id <- true_da_data$pair_id
  
  # Combine the data
  combined_data <- rbind(red_points, blue_points)
  
  
  surfaces_opacity <- 0.03
  # Calculate finite ranges for surfaces
  finite_log2FoldChange <- result_data$log2FoldChange[is.finite(result_data$log2FoldChange)]
  finite_taxa_prop <- result_data$taxa_prop[is.finite(result_data$taxa_prop)]
  finite_log10padj <- result_data$log10padj[is.finite(result_data$log10padj)]
  
  # Create vertical surfaces at log2FoldChange = -1.5 and 1.5
  y_vals <- seq(min(finite_taxa_prop), max(finite_taxa_prop), length.out = 100)
  z_vals <- seq(min(finite_log10padj), max(finite_log10padj), length.out = 100)
  
  # Define ranges for the horizontal surfaces
  x_vals_left <- seq(min(finite_log2FoldChange), -1.5, length.out = 100)
  x_vals_right <- seq(1.5, max(finite_log2FoldChange), length.out = 100)
  y_vals_horizontal <- seq(min(finite_taxa_prop), max(finite_taxa_prop), length.out = 100)
  
  
  
  # Create the initial scatter plot
  p <- plot_ly(
    data = result_data, 
    x = ~log2FoldChange, 
    y = ~taxa_prop, 
    z = ~log10padj, 
    size = ~taxa_prop,
    
    type = 'scatter3d', 
    mode = 'markers',
    marker = list(
      sizemode = 'diameter',
      sizeref = 3,
      color = ~rgba_color,
      line = list(width = 0)
    )
    
    # colors = c("Truly Differentially Abundant" = 'red', "Not Differentially Abundant" = 'black')
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = 'log2FoldChange'),
        yaxis = list(title = 'Taxa Proportion'),
        zaxis = list(title = 'log10 padj')
      )
    ) %>% add_surface(
      x = matrix(rep(x_vals_right, each = 100), nrow = 100, ncol = 100),
      y = matrix(rep(y_vals_horizontal, times = 100), nrow = 100, ncol = 100),
      z = matrix(2, nrow = 100, ncol = 100),
      opacity = surfaces_opacity,
      showscale = FALSE,
      colorscale = list(c(0, 'blue'), c(1, 'blue'))
    ) %>%
    add_surface(
      x = matrix(-1.5, nrow = 100, ncol = 100),
      y = matrix(rep(y_vals, each = 100), nrow = 100, ncol = 100),
      z = matrix(rep(z_vals, times = 100), nrow = 100, ncol = 100),
      opacity = surfaces_opacity,
      showscale = FALSE,
      colorscale = list(c(0, 'blue'), c(1, 'blue'))
    )%>%
    add_surface(
      x = matrix(1.5, nrow = 100, ncol = 100),
      y = matrix(rep(y_vals, each = 100), nrow = 100, ncol = 100),
      z = matrix(rep(z_vals, times = 100), nrow = 100, ncol = 100),
      opacity = surfaces_opacity,
      showscale = FALSE,
      colorscale = list(c(0, 'blue'), c(1, 'blue'))
    )%>%
    add_surface(
      x = matrix(rep(x_vals_left, each = 100), nrow = 100, ncol = 100),
      y = matrix(rep(y_vals_horizontal, times = 100), nrow = 100, ncol = 100),
      z = matrix(2, nrow = 100, ncol = 100),
      opacity = 0.05,
      showscale = FALSE,
      colorscale = list(c(0, 'blue'), c(1, 'blue'))
    )%>%
    add_surface(
      x = matrix(rep(x_vals_right, each = 100), nrow = 100, ncol = 100),
      y = matrix(rep(y_vals_horizontal, times = 100), nrow = 100, ncol = 100),
      z = matrix(2, nrow = 100, ncol = 100),
      opacity = surfaces_opacity,
      showscale = FALSE,
      colorscale = list(c(0, 'blue'), c(1, 'blue'))
    )
  
  # Add blue points for truly differentially abundant rows
  p <- p %>%
    add_trace(
      data = blue_points,
      x = ~X,
      y = ~Y,
      z = ~Z,
      size = ~Y,
      type = 'scatter3d',
      mode = 'markers',
      marker = list(
        sizemode = 'diameter',
        sizeref = 3,
        color='blue',
        opacity=0.2
      ),
      inherit = FALSE,
      showlegend = FALSE
    ) %>%
    add_trace(
      data = combined_data,
      x = ~X,
      y = ~Y,
      z = ~Z,
      split = ~pair_id,
      line = list(color = "blue"),
      type = "scatter3d",
      mode = "lines",
      showlegend = FALSE,
      inherit = FALSE
    )
  
  
  return (p)
}