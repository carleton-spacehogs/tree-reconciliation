library(ggplot2)

COG=Sys.getenv("gene_name")
clock_model=Sys.getenv("clock_model")
sym_event_f=Sys.getenv("sym_event_date_f")

out_file = sprintf("R-plots/timeline/%s-%s-timeline.png", COG, clock_model)

df = read.delim(sym_event_f, na.strings="?")
df = df[order(df$midpoint.date),]
df$index = 1:length(df[,1])

dynamic_size = min(1.3, 100/nrow(df))

ggplot(df, aes(y=midpoint.date, x=index, color=event)) +
	geom_pointrange( alpha=0.5,
		aes(ymin=left.date, ymax=right.date),
		size=dynamic_size) + #dynamic size
	scale_color_manual( values = c(
		"dup" = "#ffd166",
		"hgt" = "#06d6a0",
		"los" = "#118ab2",
		"spe" = "#ff66de")) +
	coord_flip() +
	scale_y_reverse( limits = c(4000, 0),
		breaks = seq(0, 4000, by = 500)) +
	labs( y="Million years ago",
		x=paste("Event index for", COG),
		color = "Gene Event") +
	theme_classic() +
	theme(panel.grid.minor = element_blank(), axis.text=element_text(size=9))

ggsave(out_file, width = 6.5, height = 8)

print(paste("The ggplot timeline of", COG, "has been saved to", out_file))
