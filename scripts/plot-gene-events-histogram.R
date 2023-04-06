library(ggplot2)
args <- commandArgs(TRUE)

COG=args[1]
gene_tree_method=args[2]
COG_calling_method=args[3]

sym_event_f = paste("ecceTERA_analysis/", COG, "_symmetric.events_event_dates.txt", sep = "")
plot_title = paste(COG, "Gene Events")
out_file = paste(COG_calling_method, gene_tree_method, COG, "eventsHistogram.png", sep = "-")
out_file = paste("R-plots", "histogram", out_file, sep = "/")

data <- read.delim(sym_event_f, na.strings="?")

ggplot(data, aes(x=midpoint.date)) +
	geom_histogram(aes(y=(..count../nrow(data)), fill=event), binwidth = 250, boundary = 0) + 
	xlim(4000,0) + 
	labs(x="Million years ago", 
		y="Proportion of total events",
		title=plot_title) + 
	scale_fill_manual(
		values=c("#ffd166", "#06d6a0", "#118ab2", "#ff66de"),
		name="Type of Event", 
		labels=c("Duplication", "Horizontal Gene Transfer", "Loss", "Speciation")) + 
	theme(
		legend.position = "bottom", 
		plot.title = element_text(size = rel(1.75), hjust = 0.5),
		axis.line.y = element_line(linetype = 1)) + 
	geom_hline(yintercept = 0) +
	theme_classic()

ggsave(out_file, width = 7, height = 6)

print(paste("The ggplot histogram of", COG, "events has been saved to", out_file))
