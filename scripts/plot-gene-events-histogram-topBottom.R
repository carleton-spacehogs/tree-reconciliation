library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)

COG=Sys.getenv("gene_name")
clock_model=Sys.getenv("clock_model")
sym_event_f=Sys.getenv("sym_event_date_f")

plot_title = paste(COG, "Gene Events")
out_file = sprintf("R-plots/topBottom_histogram/%s-%s-topBottom.png", COG, clock_model)

data = read.delim(sym_event_f, na.strings="?")

shared_gen_graph = function(data) {
  g = ggplot(data, aes(x=midpoint.date)) +
    geom_histogram(aes(y=after_stat(density), fill=event), binwidth = 250, boundary = 0)
  # if (type == "density") {
  #   g = ggplot(data, aes(x=midpoint.date)) +
  #    geom_density(alpha=.2, aes(fill=event))
  #}
  g + xlim(4000,0) +
    ylab("Proportion of total events") +
    theme(legend.position = "bottom",
          axis.line.y = element_line(linetype = 1)) + 
    geom_hline(yintercept = 0) +
    theme_classic()
}

gen_graph = function(events){
  top = shared_gen_graph(filter(events, event %in% c("dup", "spe"))) +
    theme(plot.title = element_text(size = rel(1.75), hjust = 0.5),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title=plot_title) + 
    scale_fill_manual(values = c("dup" = "#06d6a0", "spe" = "#ff66de"))
  bottom = shared_gen_graph(filter(events, event %in% c("hgt", "los"))) + 
    scale_y_continuous(trans = "reverse") + 
    scale_fill_manual( values = c("hgt" = "#ffd166", "los" = "#118ab2"))
  plot_grid(top, bottom, ncol = 1, align = "v", axis = "tb")
}

gen_graph(data)

ggsave(out_file, width = 6, height = 6)

print(paste("The ggplot histogram of", COG, "events has been saved to", out_file))
