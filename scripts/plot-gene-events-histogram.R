library(ggplot2)

COG=Sys.getenv("gene_name")
gene_tree_method=Sys.getenv("gene_tree_method")
COG_calling_method=Sys.getenv("COG_calling_method")
sym_event_f=Sys.getenv("sym_event_date_f")

plot_title = paste(COG, "Gene Events")
out_file = paste(COG_calling_method, gene_tree_method, COG, "eventsHistogram.png", sep = "-")
out_file = paste("R-plots", "histogram", out_file, sep = "/")

data = read.delim(sym_event_f, na.strings="?")

ggplot(data, aes(x=midpoint.date)) +
	geom_histogram(aes(y=(..count../nrow(data)), fill=event), binwidth = 250, boundary = 0) + 
	xlim(4000,0) + 
	labs(x="Million years ago", 
		y="Proportion of total events",
		title=plot_title) + 
	scale_fill_manual(
		values = c("dup" = "#ffd166", "hgt" = "#06d6a0",
               "los" = "#118ab2", "spe" = "#ff66de")) +
	theme(
		legend.position = "bottom", 
		plot.title = element_text(size = rel(1.75), hjust = 0.5),
		axis.line.y = element_line(linetype = 1)) + 
	geom_hline(yintercept = 0) +
	theme_classic()

ggsave(out_file, width = 7, height = 6)

print(paste("The ggplot histogram of", COG, "events has been saved to", out_file))
