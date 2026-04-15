#!/usr/bin/env Rscript
# point_distribution_new.R - Standard transcript site distribution plot script (new version)
#
# Features:
#   - Analyze and plot site distribution in different transcript regions (5'UTR, CDS, 3'UTR)
#   - Support multiple output formats and resolution settings
#   - Provide professional graph design and statistical information
#   - Enhanced error handling and data validation
#   - Support custom configuration options
#
# Usage:
#   Rscript point_distribution_new.R <input_file> <output_file> [options]
#
# Input file format:
#   Tab-separated text file with the following columns:
#   V1: Region type (five_utr, cds, three_utr)
#   V2: Position
#   V3: Site count

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("--format", "-f"), type="character", default="pdf",
              help="Output file format [default: %default]", metavar="format"),
  make_option(c("--dpi", "-d"), type="integer", default=300,
              help="Output resolution [default: %default]", metavar="dpi"),
  make_option(c("--width", "-w"), type="numeric", default=84,  # 3.3 inches = 84 mm
              help="Plot width (mm) [default: %default]", metavar="width"),
  make_option(c("--height", "-H"), type="numeric", default=38,  # 1.5 inches = 38 mm
              help="Plot height (mm) [default: %default]", metavar="height"),
  make_option(c("--show-legend", "-l"), action="store_true", default=FALSE,
              help="Show legend [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=2)

input_file <- opt$args[1]
output_file <- opt$args[2]

# 处理参数，确保非空值
output_format <- if (is.null(opt$options$format)) "pdf" else opt$options$format
output_dpi <- if (is.null(opt$options$dpi)) 300 else opt$options$dpi
# 将毫米转换为英寸 (1 inch = 25.4 mm)
plot_width_mm <- if (is.null(opt$options$width)) 84 else opt$options$width
plot_height_mm <- if (is.null(opt$options$height)) 38 else opt$options$height
plot_width <- plot_width_mm / 25.4
plot_height <- plot_height_mm / 25.4
show_legend <- if (is.null(opt$options$show_legend)) FALSE else opt$options$show_legend

# Data processing function
process_data <- function(file_path) {
  # Read data
  cat("Reading input file...\n")
  dist_raw <- read.table(file_path, sep = "\t", header = FALSE)
  
  # Data validation
  if (ncol(dist_raw) < 3) {
    stop("Input file format error: at least 3 columns required")
  }
  
  # Ensure correct column names
  colnames(dist_raw) <- c("region", "position", "count")
  
  # Validate region types
  valid_regions <- c("five_utr", "cds", "three_utr")
  if (!all(dist_raw$region %in% valid_regions)) {
    invalid_regions <- unique(dist_raw$region[!dist_raw$region %in% valid_regions])
    warning(paste("Found invalid region types:", paste(invalid_regions, collapse=", ")))
    # Filter invalid data
    dist_raw <- dist_raw[dist_raw$region %in% valid_regions, ]
  }
  
  # Ensure data for each region
  if (length(unique(dist_raw$region)) < 3) {
    missing_regions <- setdiff(valid_regions, unique(dist_raw$region))
    warning(paste("Missing region data:", paste(missing_regions, collapse=", ")))
  }
  
  # Data quality check
  if (any(dist_raw$count < 0)) {
    warning("Found negative site counts, automatically filtered")
    dist_raw <- dist_raw[dist_raw$count >= 0, ]
  }
  
  if (nrow(dist_raw) == 0) {
    stop("No valid data available for analysis")
  }
  
  cat(paste("Successfully read", nrow(dist_raw), "rows of data\n"))
  return(dist_raw)
}

# Calculate statistics
calculate_statistics <- function(data) {
  cat("Calculating statistics...\n")
  
  # Calculate maximum position for each region
  x_max <- sapply(c("five_utr", "cds", "three_utr"), function(region) {
    subset_data <- data[data$region == region, ]
    if (nrow(subset_data) > 0) {
      return(ceiling(max(subset_data$position) / 5) * 5)
    } else {
      return(0)
    }
  })
  
  # Calculate maximum count for each region
  y_max <- sapply(c("five_utr", "cds", "three_utr"), function(region) {
    subset_data <- data[data$region == region, ]
    if (nrow(subset_data) > 0) {
      return(ceiling(max(subset_data$count) / 5) * 5)
    } else {
      return(0)
    }
  })
  
  # Calculate total count for each region
  sum_counts <- sapply(c("five_utr", "cds", "three_utr"), function(region) {
    subset_data <- data[data$region == region, ]
    if (nrow(subset_data) > 0) {
      return(sum(subset_data$count))
    } else {
      return(0)
    }
  })
  
  total_count <- sum(sum_counts)
  
  # Calculate proportions
  proportions <- if (total_count > 0) {
    round((sum_counts / total_count) * 100, 2)
  } else {
    c(0, 0, 0)
  }
  
  # Calculate global maximum Y value
  global_y_max <- ceiling((max(data$count, na.rm = TRUE) + 1) / 5) * 5
  
  return(list(
    x_max = x_max,
    y_max = y_max,
    sum_counts = sum_counts,
    proportions = proportions,
    total_count = total_count,
    global_y_max = global_y_max
  ))
}

# Prepare plot data
prepare_plot_data <- function(data, stats) {
  cat("Preparing plot data...\n")
  
  # Extract data for each region
  five_utr_data <- data[data$region == "five_utr", ]
  cds_data <- data[data$region == "cds", ]
  three_utr_data <- data[data$region == "three_utr", ]
  
  # Calculate text label positions
  text_positions <- data.frame(
    x = c(
      stats$x_max["five_utr"] / 2,
      stats$x_max["five_utr"] + stats$x_max["cds"] / 2,
      stats$x_max["five_utr"] + stats$x_max["cds"] + stats$x_max["three_utr"] / 2
    ),
    y = c(
      stats$y_max["five_utr"] * 0.9,
      stats$y_max["cds"] * 0.9,
      stats$y_max["three_utr"] * 0.9
    ),
    text = c(
      paste0(stats$sum_counts["five_utr"], "\n", stats$proportions["five_utr"], "%"),
      paste0(stats$sum_counts["cds"], "\n", stats$proportions["cds"], "%"),
      paste0(stats$sum_counts["three_utr"], "\n", stats$proportions["three_utr"], "%")
    ),
    region = c("five_utr", "cds", "three_utr")
  )
  
  return(list(
    five_utr = five_utr_data,
    cds = cds_data,
    three_utr = three_utr_data,
    text_data = text_positions
  ))
}

# Generate plot
generate_plot <- function(plot_data, stats, show_legend = FALSE) {
  cat("Generating plot...\n")
  
  # Define professional color scheme
  colors <- c(
    "five_utr" = "#4DAF4A",  # Green
    "cds" = "#377EB8",        # Blue
    "three_utr" = "#E41A1C"    # Red
  )
  
  # Calculate total width without gaps
  total_width_no_gaps <- sum(stats$x_max)
  
  # Calculate gap size (2% of total width with gaps)
  # Let total_width = total_width_no_gaps + 2 * gap
  # gap = 0.02 * total_width
  # Solving: total_width = total_width_no_gaps + 2 * 0.02 * total_width
  # total_width - 0.04 * total_width = total_width_no_gaps
  # 0.96 * total_width = total_width_no_gaps
  # total_width = total_width_no_gaps / 0.96
  total_width <- total_width_no_gaps / 0.96
  gap <- 0.02 * total_width
  
  # Calculate region positions with gaps
  five_utr_end <- stats$x_max["five_utr"]
  cds_start <- five_utr_end + gap
  cds_end <- cds_start + stats$x_max["cds"]
  three_utr_start <- cds_end + gap
  three_utr_end <- three_utr_start + stats$x_max["three_utr"]
  
  # Create base plot
  p <- ggplot() +
    # 5'UTR curve
    geom_line(
      data = plot_data$five_utr,
      aes(x = position, y = count),
      colour = colors["five_utr"],
      linewidth = 0.4,  # 4 pt line
      show.legend = show_legend
    ) +
    # CDS curve
    geom_line(
      data = plot_data$cds,
      aes(x = position + cds_start, y = count),
      colour = colors["cds"],
      linewidth = 0.4,  # 4 pt line
      show.legend = show_legend
    ) +
    # 3'UTR curve
    geom_line(
      data = plot_data$three_utr,
      aes(x = position + three_utr_start, y = count),
      colour = colors["three_utr"],
      linewidth = 0.4,  # 4 pt line
      show.legend = show_legend
    ) +
    # Text labels (5pt font)
    geom_text(
      data = plot_data$text_data,
      aes(
        x = ifelse(region == "five_utr", x, 
               ifelse(region == "cds", x + gap, x + 2 * gap)),
        y = y,
        label = text,
        fontface = "bold",
        vjust = 0
      ),
      colour = colors[plot_data$text_data$region],
      size = 0.833,  # 2.5 pt font (R uses points as 1/72 inch)
      show.legend = FALSE
    ) +
    # Add horizontal region lines (6 pt thickness)
    geom_segment(
      aes(x = 0, y = -1, xend = five_utr_end, yend = -1),
      colour = colors["five_utr"],
      linewidth = 0.6,  # 6 pt line
      show.legend = FALSE
    ) +
    geom_segment(
      aes(x = cds_start, y = -1, xend = cds_end, yend = -1),
      colour = colors["cds"],
      linewidth = 0.6,  # 6 pt line
      show.legend = FALSE
    ) +
    geom_segment(
      aes(x = three_utr_start, y = -1, xend = three_utr_end, yend = -1),
      colour = colors["three_utr"],
      linewidth = 0.6,  # 6 pt line
      show.legend = FALSE
    ) +
    # Add region labels at the bottom
    annotate("text", 
             x = five_utr_end / 2, 
             y = -stats$global_y_max * 0.05-1, 
             label = "5'UTR", 
             fontface = "bold", 
             size = 0.833,  # 5 pt font
             color = colors["five_utr"]
    ) +
    annotate("text", 
             x = cds_start + (stats$x_max["cds"] / 2), 
             y = -stats$global_y_max * 0.05-1, 
             label = "CDS", 
             fontface = "bold", 
             size = 0.833,  # 5 pt font
             color = colors["cds"]
    ) +
    annotate("text", 
             x = three_utr_start + (stats$x_max["three_utr"] / 2), 
             y = -stats$global_y_max * 0.05-1, 
             label = "3'UTR", 
             fontface = "bold", 
             size = 0.833,  # 5 pt font
             color = colors["three_utr"]
    ) +
    # Axis labels
    labs(
      x = NULL,
      y = "Site Count"
    ) +
    # Y-axis settings
    scale_y_continuous(
      limits = c(-stats$global_y_max * 0.15, stats$global_y_max * 1.1),
      breaks = seq(0, stats$global_y_max, by = max(1, floor(stats$global_y_max / 5))),
      expand = c(0, 0)
    ) +
    # X-axis settings
    scale_x_continuous(
      limits = c(0, total_width),
      breaks = NULL,  # Remove x-axis ticks
      expand = c(0, 0)
    ) +
    # Theme settings
    theme_bw() +
    theme(
      # Plot margins
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      # Grid lines
      panel.grid.major = element_line(color = "gray90", linewidth = 0.025),  # 0.25 pt line
      panel.grid.minor = element_blank(),
      # Panel border
      panel.border = element_rect(color = "black", linewidth = 0.0125),  # 0.125 pt line
      # Axes
      axis.line = element_line(color = "black", linewidth = 0.025),  # 0.25 pt line
      axis.ticks = element_line(color = "black", linewidth = 0.025),  # 0.25 pt line
      axis.ticks.length = unit(0.5, "mm"),  # 0.5 mm ticks
      # Axis labels
      axis.title.y = element_text(face = "bold", size = 2),  # 6 pt font
      axis.text.y = element_text(size = 1.667, color = "black"),  # 5 pt font
      axis.text.x = element_blank(),  # Remove x-axis text
      # Font - use system default font
      text = element_text(size = 1.667)  # 5 pt font
    )
  
  # If legend is needed
  if (show_legend) {
    # Create legend data
    legend_data <- data.frame(
      region = c("5'UTR", "CDS", "3'UTR"),
      color = c(colors["five_utr"], colors["cds"], colors["three_utr"])
    )
    
    # Add legend
    p <- p +
      geom_point(
        data = legend_data,
        aes(x = 0, y = 0, color = region),
        show.legend = TRUE
      ) +
      scale_color_manual(
        values = c("5'UTR" = colors["five_utr"], "CDS" = colors["cds"], "3'UTR" = colors["three_utr"])
      ) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank()
      )
  }
  
  return(p)
}

# Main function
main <- function() {
  cat("=== Site Distribution Analysis (New Version) ===\n")
  cat(paste("Input file:", input_file, "\n"))
  cat(paste("Output file:", output_file, "\n"))
  cat(paste("Output format:", output_format, "\n"))
  cat(paste("Resolution:", output_dpi, "dpi\n"))
  cat(paste("Plot size:", plot_width_mm, "x", plot_height_mm, "mm\n"))
  cat(paste("Show legend:", ifelse(show_legend, "Yes", "No"), "\n"))
  
  # Ensure show_legend is a boolean
  if (length(show_legend) == 0 || is.na(show_legend)) {
    show_legend <- FALSE
  } else if (is.character(show_legend)) {
    show_legend <- as.logical(show_legend)
    if (is.na(show_legend)) {
      show_legend <- FALSE
    }
  } else if (!is.logical(show_legend)) {
    show_legend <- FALSE
  }
  
  tryCatch({
    # Process data
    dist_raw <- process_data(input_file)
    
    # Calculate statistics
    stats <- calculate_statistics(dist_raw)
    
    # Prepare plot data
    plot_data <- prepare_plot_data(dist_raw, stats)
    
    # Generate plot
    plot <- generate_plot(plot_data, stats, show_legend)
    
    # Save plot
    cat("Saving plot...\n")
    ggsave(
      output_file,
      plot = plot,
      device = output_format,
      width = plot_width,
      height = plot_height,
      dpi = output_dpi,
      units = "in"
    )
    
    # Output statistics
    cat("\n=== Statistics ===\n")
    cat(paste("5'UTR: ", stats$sum_counts["five_utr"], " (", stats$proportions["five_utr"], "%)\n", sep=""))
    cat(paste("CDS:    ", stats$sum_counts["cds"], " (", stats$proportions["cds"], "%)\n", sep=""))
    cat(paste("3'UTR:  ", stats$sum_counts["three_utr"], " (", stats$proportions["three_utr"], "%)\n", sep=""))
    cat(paste("Total:   ", stats$total_count, "\n", sep=""))
    
    cat("\nAnalysis completed! Plot saved to:", output_file, "\n")
  }, error = function(e) {
    cat("\n[Error]: ", conditionMessage(e), "\n")
    cat("Analysis failed. Please check input file and parameter settings.\n")
    q(status = 1)
  }, warning = function(w) {
    cat("\n[Warning]: ", conditionMessage(w), "\n")
  })
}

# Run main function
main()
