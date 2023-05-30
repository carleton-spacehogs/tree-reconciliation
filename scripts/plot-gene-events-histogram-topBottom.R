library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)

COG=Sys.getenv("gene_name")
clock_model=Sys.getenv("clock_model")
sym_event_f=Sys.getenv("sym_event_date_f")

source("specific_substrates/plot-substrate-utils.R")
plot_node = "midpoint"
type="histogram"

out_file = sprintf("R-plots/topBottom_histogram/%s-%s-topBottom.png", COG, clock_model)
data = read.delim(sym_event_f, na.strings="?")

gen_graph(data, plot_title = paste(COG, "Gene Events"), location = FALSE, binwidth = 250)

ggsave(out_file, width = 6, height = 6)

print(paste("The ggplot histogram of", COG, "events has been saved to", out_file))
