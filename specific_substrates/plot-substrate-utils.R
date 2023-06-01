mixed_states = rev(c("not_sure", "marine_deep", "marine_deep.marine_shallow", "marine_deep.terrestrial", "marine_shallow", "marine_shallow.terrestrial", "terrestrial"))
mixed_colors = rev(c("light grey", "blue","steel blue","yellow","light blue","orange","red"))
mixed_color_scale = setNames(mixed_colors, mixed_states)

event_states = c("dup", "spe", "hgt", "los")
event_colors = c("#06d6a0", "#ff66de","#ffd166","#118ab2")
event_color_scale = setNames(event_colors, event_states)

remove_leaf_events = FALSE

get_COGs_of = function(key) {
  COGs = str_split(categories_COGs[key,], " ")[[1]]
  bad_COGs = setdiff(COGs, good_COGs)
  print(paste("I don't have these COGs", paste(bad_COGs, collapse = " ")))
  print("These COG failed:")
  fail_COGs = intersect(bad_COGs, done_COGs)
  print(COG_record %>% filter(COG %in% fail_COGs) %>% select(COG, exit_status) )
  print(paste("I have not run to these COGs yet:", paste(setdiff(bad_COGs, done_COGs), collapse = " ")))
  good_COGs = intersect(COGs, good_COGs)
  print(paste("I have these good COGs for :", key))
  print(paste(good_COGs, collapse = " "))
  return(good_COGs)
}

read_events = function(file, location = FALSE) {
  out_cols = c("event", "left.date", "midpoint.date", "right.date", "COG")
  if (location) {
    out_cols = c(out_cols, "left.origin", "right.origin")
  }
  # print(file)
  if (file.exists(file)) {
    all = read.delim(file, na.strings = "?")
    if (remove_leaf_events) {
      all = all[all$right.date > 0, ] 
    }
    if (nrow(all) != 0) {
      COG = str_match(file, "COG\\d+")[1]
      all = all %>% filter(right.node != "UNSAMPLED") %>% mutate(COG = COG)
      return(all[out_cols]) 
    }
  }
  return()
}

merge_events = function(COGs, location = FALSE) {
  file_list = paste0(sprintf("../%s_ecceTERA_analysis/", clock_model), COGs, "_symmetric.events_event_dates.txt")
  if (location) {
    file_list = paste0(sprintf("../events_locations/T29_%s_location_merged_events/", clock_model), COGs, "-events-location.tsv")
  }
  data_list = lapply(file_list, read_events, location = location)
  merge_data = bind_rows(data_list)
  return(merge_data)
}

gen_graph = function(events, plot_title = "merged events", location = FALSE, binwidth = 50){
  top_df = filter(events, event %in% event_states[1:2])
  low_df = filter(events, event %in% event_states[3:4])
  
  if (location) {
    sel_row1 = events[, paste0(plot_node, ".origin")] %in% mixed_states[1:3]
    sel_row2 = events[, paste0(plot_node, ".origin")] %in% mixed_states[4:length(mixed_states)]
    top_df = events[sel_row1, ]
    low_df = events[sel_row2, ]
  }
  
  top = gen_g_helper(top_df, location = location, binwidth = binwidth)
  low = gen_g_helper(low_df, location = location, binwidth = binwidth)
  max_y = max(c(ggplot_build(top)$data[[1]]$y), c(ggplot_build(low)$data[[1]]$y))
  
  if (location) {
    top = top + scale_fill_manual(values = mixed_color_scale[1:3])
    low = low + scale_fill_manual(values = mixed_color_scale[4:length(mixed_states)])
  } else {
    top = top + scale_fill_manual(values = event_color_scale[1:2])
    low = low + scale_fill_manual(values = event_color_scale[3:4])
  }
  
  top = top + scale_y_continuous(limits = c(0, max_y)) +
    labs(title=plot_title) +
    theme(plot.title = element_text(size = rel(1.75), hjust = 0.5),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  
  low = low + scale_y_continuous(limits = c(max_y, 0), trans = "reverse")
  
  plot_grid(top, low, ncol = 1, align = "v", axis = "tb")
}

gen_g_helper = function(data, location = FALSE, binwidth = 50) {
  fill_with = "event"
  if (location) {
    fill_with = paste(plot_node, "origin", sep = ".")
  }
  
  g = ggplot(data, aes_string(x=paste(plot_node, "date", sep = ".")))
  if (type == "histogram") {
    g = g + geom_histogram(aes_string(y=("..count../nrow(data)"), fill=fill_with), binwidth = binwidth, boundary = 0)
  } else if (type == "density") {
    g = g + geom_density(alpha=.2, aes_string(fill=fill_with))
  } else {
    print("type can only be histogram or density, currently, you are:")
    print(type)
  }
  g + xlim(4000,0) +
    ylab("Proportion of total events") +
    theme(legend.position = "bottom",
          axis.line.y = element_line(linetype = 1)) + 
    geom_hline(yintercept = 0) +
    theme_classic()
}

get_events = function(feature, anno, location = FALSE) {
  COGs = get_COGs_of(feature)
  events = merge_events(COGs, location)
  events$substrate = anno
  return(events)
}

make_graph = function(events, out_graph, title_text, location = FALSE) {
  num_COGs = length(unique(events$COG))
  print(num_COGs)
  pt = sprintf("%s (n = %d)", title_text, num_COGs)
  gen_graph(events, plot_title = pt, location)
  ggsave(out_graph, width = 6, height = 5)
}

jitter_plot_shared = function(g, ytop) {
  g + geom_jitter(alpha = 0.3, height = 0.3) +
    scale_x_reverse(limits = c(4000, 0), breaks = seq(0, 4000, by = 500)) +
    labs(x = "Million years ago", y = "") +
    # geom_hline(yintercept = 20.5) +
    # geom_hline(yintercept = 19.5, color = "lightgray") +
    geom_jitter(alpha = 0.3, height = 0.3) +
    scale_x_reverse(limits = c(4000, 0),
                    breaks = seq(0, 4000, by = 500)) +
    annotate('rect', xmin=2500, xmax=2300, ymin=0, ymax=ytop + 1, alpha=.3, fill='gray') +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
}
