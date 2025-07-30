#0. place this script under /workspace/data/Space_Hogs_shared_workspace/tree-reconciliation/specific_substrates

#1. parameters and imports
library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)
library(scales)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("plot-substrate-utils.R")
state_vec_short = c("A", "B", "C", "AC", "BC")
state_vec_long = c("marine_deep","marine_shallow", "terrestrial", "marine_deep","marine_shallow")
clock_model = "cir1"
remove_leaf_events = FALSE
type="histogram"
trial <- "master_46_res"
#TODO: update the csv file to specify what COGs to plot
#have to use the first column as index
categories_COGs = read.csv("all-cor-COGs-2024.csv", header = TRUE, sep = ",", row.names = 1)
#this is ALL the COGs, and i don't have their slopes in this file. will have to add that so we only have the high correlation COGs, i guess, or use tony's script
COG_record = read.csv(sprintf("../COG_reconciliation_summary_%s.csv", clock_model))
# this is the record of which COGs have been run
done_COGs = pull(COG_record, COG)
good_COGs = filter(COG_record, exit_status == "everything_successful") %>% pull(COG)

#2. do proportion plot
calculate_event_proportions <- function(clock_model, trial, bin_width = 10) {
  COG_record <- read.csv(sprintf("../COG_reconciliation_summary_%s.csv", clock_model))
  all_good_COGs <- filter(COG_record, exit_status == "everything_successful") %>% pull(COG)
  
  cat(sprintf("\nCalculating total events for %d COGs...\n", length(all_good_COGs)))
  
  all_events <- merge_events(all_good_COGs, location = FALSE)
  all_events <- all_events[all_events$event %in% c("spe", "hgt", "dup", "los"), ]
  
  time_bins <- seq(0, 5000, by = bin_width)
  all_events$date <- all_events$left.date
  all_events$time_bin <- cut(all_events$date, breaks = time_bins, include.lowest = TRUE)
  
  total_counts_per_bin <- all_events %>%
    group_by(time_bin) %>%
    summarise(total_count = n(), .groups = "drop")
  
  substrate_correlations <- list(
    "positive_correlation_oxygen" = "Oxygen",
    "negative_correlation_oxygen" = "Oxygen",
    "positive_correlation_iron" = "Iron", 
    "negative_correlation_iron" = "Iron",
    "positive_correlation_phosphate" = "Phosphate",
    "negative_correlation_phosphate" = "Phosphate",
    "positive_correlation_NO2NO3" = "Nitrite/Nitrate",
    "negative_correlation_NO2NO3" = "Nitrite/Nitrate",
    "Fe2+" = "Ferrous Iron",
    "Fe3+" = "Ferric Iron"
  )
  
  results <- list()
  
  for (correlation_key in names(substrate_correlations)) {
    substrate_name <- substrate_correlations[[correlation_key]]
    correlation_type <- ifelse(grepl("positive", correlation_key), "positive",
                               ifelse(grepl("negative", correlation_key), "negative", "neutral"))
    
    cat(sprintf("\nProcessing %s...\n", correlation_key))
    
    cogs <- get_COGs_of(correlation_key)
    
    if (length(cogs) > 0) {
      substrate_events <- merge_events(cogs, location = FALSE)
      substrate_events <- substrate_events[substrate_events$event %in% c("spe", "hgt", "dup", "los"), ]
      
      if (nrow(substrate_events) > 0) {
        substrate_events$date <- substrate_events$left.date
        substrate_events$time_bin <- cut(substrate_events$date, breaks = time_bins,
                                         include.lowest = TRUE)
        
        substrate_counts <- substrate_events %>%
          group_by(time_bin, event) %>%
          summarise(substrate_count = n(), .groups = "drop")
        
        proportions <- substrate_counts %>%
          left_join(total_counts_per_bin, by = "time_bin") %>%
          mutate(
            proportion = substrate_count / total_count,
            substrate = substrate_name,
            correlation_type = correlation_type,
            correlation_key = correlation_key,
            midpoint = time_bins[as.numeric(time_bin)] + bin_width/2,
            proportion_signed = ifelse(event == "los", -proportion, proportion)
          )
        
        results[[correlation_key]] <- proportions
      }
    }
  }
  
  all_proportions <- bind_rows(results)
  
  event_types <- c("spe", "hgt", "dup", "los")
  complete_grid <- expand.grid(
    correlation_key = names(substrate_correlations),
    event = event_types,
    time_bin = levels(all_events$time_bin),
    stringsAsFactors = FALSE
  )
  
  complete_grid$substrate <- unlist(substrate_correlations[complete_grid$correlation_key])
  complete_grid$correlation_type <- ifelse(grepl("positive", complete_grid$correlation_key), "positive",
                                           ifelse(grepl("negative", complete_grid$correlation_key), "negative", "neutral"))
  
  all_proportions <- complete_grid %>%
    left_join(all_proportions, by = c("correlation_key", "event", "time_bin")) %>%
    mutate(
      substrate = ifelse(is.na(substrate.y), substrate.x, substrate.y),
      correlation_type = ifelse(is.na(correlation_type.y), correlation_type.x, correlation_type.y),
      proportion = ifelse(is.na(proportion), 0, proportion),
      proportion_signed = case_when(
        event == "los" ~ -proportion,
        TRUE ~ proportion
      ),
      substrate_count = ifelse(is.na(substrate_count), 0, substrate_count),
      midpoint = time_bins[match(time_bin, unique(time_bin))] + bin_width/2
    ) %>%
    select(-substrate.x, -substrate.y, -correlation_type.x, -correlation_type.y)
  
  return(all_proportions)
}

create_nutrient_proportion_plots <- function(clock_models = c("ugam1", "cir1", "ln3"),
                                             trial = "master_46_res",
                                             bin_width = 100) {
  
  all_results <- list()
  
  for (clock_model in clock_models) {
    cat(sprintf("\n\nclock model: %s \n", clock_model))
    
    proportions <- calculate_event_proportions(clock_model, trial, bin_width)
    proportions$clock_model <- clock_model
    
    all_results[[clock_model]] <- proportions
  }
  
  all_data <- bind_rows(all_results)

#modified this too
  
  event_colors <- c(
    "marine_deep" = "#A8C4E6",
    "marine_shallow" = "#B8D8B8",
    "terrestrial" = "#E6A8A8"
  )
  
  for (cm in clock_models) {
    cm_data <- all_data %>% filter(clock_model == cm)
    
    create_symmetric_breaks <- function(data) {
      max_val <- max(abs(data$proportion_signed), na.rm = TRUE) * 100
      
      if (max_val <= 3) {
        break_interval <- 1
        max_break <- 3
      } else if (max_val <= 6) {
        break_interval <- 1
        max_break <- ceiling(max_val)
      } else if (max_val <= 10) {
        break_interval <- 2
        max_break <- ceiling(max_val / break_interval) * break_interval
      } else {
        break_interval <- 5
        max_break <- ceiling(max_val / break_interval) * break_interval
      }
      
      pos_breaks <- seq(0, max_break, by = break_interval)
      
      plot_breaks <- c(-rev(pos_breaks[-1]), pos_breaks)
      
      labels <- paste0(abs(plot_breaks), "%")
      
      return(list(breaks = plot_breaks/100, labels = labels))
    }
    
    pos_data <- cm_data %>% filter(correlation_type == "positive")
    
#i modified the below so that fill = left.origin instead of event, and modified the scale_fill_manual to match the habitats    

    if (nrow(pos_data) > 0) {
      axis_info <- create_symmetric_breaks(pos_data)
      
      p_pos <- ggplot(pos_data, aes(x = midpoint, y = proportion_signed, fill = left.origin)) +
        geom_col(position = "stack", width = bin_width * 0.8, alpha = 0.8) +
        
        scale_x_reverse(limits = c(5000, 0), breaks = seq(0, 5000, by = 500)) +
        scale_y_continuous(breaks = axis_info$breaks, labels = axis_info$labels) +
        scale_fill_manual(values = event_colors,
                          breaks = c("marine_shallow", "marine_deep", "terrestrial"),
                          labels = c("Shallow marine", "Deep marine", "Terrestrial")) +
        
        facet_wrap(~ paste("Positive Correlation:", substrate), ncol = 2, scales = "free_y") +
        
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        
        annotate('rect', xmin = 2300, xmax = 2500, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = 'gray') +
        annotate('text', x = 2400, y = Inf, label = "GOE", vjust = 2, size = 2) +
        annotate('rect', xmin = 540, xmax = 850, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = 'gray') +
        annotate('text', x = 695, y = Inf, label = "NOE", vjust = 2, size = 2) +
        
        labs(
          x = "Time (Mya)",
          y = "Proportion of All Gene Events",
          title = cm,
          fill = "Event Type"
        ) +
        
        theme_minimal() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
          strip.text = element_text(size = 10, face = "bold"),
          legend.position = "bottom"
        )
      
      ggsave(sprintf("%s/%s-plots/nutrient_proportions_positive_correlations.svg", trial, cm),
             plot = p_pos, width = 12, height = 10, device = "svg")
    }
    
    neg_data <- cm_data %>% filter(correlation_type == "negative")
    
    if (nrow(neg_data) > 0) {
      axis_info <- create_symmetric_breaks(neg_data)
      
      p_neg <- ggplot(neg_data, aes(x = midpoint, y = proportion_signed, fill = event)) +
        geom_col(position = "stack", width = bin_width * 0.8, alpha = 0.8) +
        
        scale_x_reverse(limits = c(5000, 0), breaks = seq(0, 5000, by = 500)) +
        scale_y_continuous(breaks = axis_info$breaks, labels = axis_info$labels) +
        scale_fill_manual(values = event_colors,
                          breaks = c("spe", "hgt", "dup", "los"),
                          labels = c("Speciation", "HGT", "Duplication", "Loss")) +
        
        facet_wrap(~ paste("Negative Correlation:", substrate), ncol = 2, scales = "free_y") +
        
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        
        annotate('rect', xmin = 2300, xmax = 2500, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = 'gray') +
        annotate('text', x = 2400, y = Inf, label = "GOE", vjust = 2, size = 2) +
        annotate('rect', xmin = 540, xmax = 850, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = 'gray') +
        annotate('text', x = 695, y = Inf, label = "NOE", vjust = 2, size = 2) +
        
        labs(
          x = "Time (Mya)",
          y = "Proportion of All Gene Events",
          title = cm,
          fill = "Event Type"
        ) +
        
        theme_minimal() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
          strip.text = element_text(size = 10, face = "bold"),
          legend.position = "bottom"
        )
      
      ggsave(sprintf("%s/%s-plots/nutrient_proportions_negative_correlations.svg", trial, cm),
             plot = p_neg, width = 12, height = 10, device = "svg")
    }
    
    neutral_data <- cm_data %>% filter(correlation_type == "neutral")
    
    if (nrow(neutral_data) > 0) {
      axis_info <- create_symmetric_breaks(neutral_data)
      
      p_neutral <- ggplot(neutral_data, aes(x = midpoint, y = proportion_signed, fill = event)) +
        geom_col(position = "stack", width = bin_width * 0.8, alpha = 0.8) +
        
        scale_x_reverse(limits = c(5000, 0), breaks = seq(0, 5000, by = 500)) +
        scale_y_continuous(breaks = axis_info$breaks, labels = axis_info$labels) +
        scale_fill_manual(values = event_colors,
                          breaks = c("spe", "hgt", "dup", "los"),
                          labels = c("Speciation", "HGT", "Duplication", "Loss")) +
        
        facet_wrap(~ substrate, ncol = 2, scales = "free_y") +
        
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        
        annotate('rect', xmin = 2300, xmax = 2500, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = 'gray') +
        annotate('text', x = 2400, y = Inf, label = "GOE", vjust = 2, size = 2) +
        annotate('rect', xmin = 540, xmax = 850, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = 'gray') +
        annotate('text', x = 695, y = Inf, label = "NOE", vjust = 2, size = 2) +
        
        labs(
          x = "Time (Mya)",
          y = "Proportion of All Gene Events",
          title = cm,
          fill = "Event Type"
        ) +
        
        theme_minimal() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
          strip.text = element_text(size = 10, face = "bold"),
          legend.position = "bottom"
        )
      
      ggsave(sprintf("%s/%s-plots/nutrient_proportions_iron_species.svg", trial, cm),
             plot = p_neutral, width = 10, height = 6, device = "svg")
    }
  }
  
  summary_stats <- all_data %>%
    group_by(clock_model, correlation_key, substrate, correlation_type, event) %>%
    summarise(
      mean_proportion = mean(proportion, na.rm = TRUE),
      max_proportion = max(proportion, na.rm = TRUE),
      time_of_max = midpoint[which.max(proportion)],
      .groups = "drop"
    ) %>%
    mutate(
      mean_proportion = sprintf("%.3f%%", mean_proportion * 100),
      max_proportion = sprintf("%.3f%%", max_proportion * 100)
    )
  
  write.csv(summary_stats, sprintf("%s/proportion_stats_raw.csv", trial),
            row.names = FALSE)
  
  return(list(data = all_data, stats = summary_stats))
}

remove_leaf_events <- TRUE
plot_node <- "midpoint"
create_nutrient_proportion_plots(
  clock_models = c("ugam1", "cir1", "ln3"),
  trial = "master_46_res",
  bin_width = 75)
