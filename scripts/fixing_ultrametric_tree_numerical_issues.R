library(adephylo)
library(ape)

cyano_tree = ape::read.tree(file ="ugam1_CorrCal2.71_sample_pruned.chronogram")
tip.heights<-distRoot(cyano_tree)

# learned from https://hcliedtke.github.io/R-scrapheap/be_ultrametric.html
heights.summary<-table(tip.heights)

options(digits=22) # set to maximum allowed digits
real.tree.height<-as.numeric(names(which.max(heights.summary)))
over.under<-tip.heights-real.tree.height

## extract all terminal edges for tips that do not have the final height we want:
tip.ids <- cyano_tree$edge[, 2] <= Ntip(cyano_tree)
terminal.edges <- cyano_tree$edge.length[tip.ids]

## add/subtract the extra length from the terminal branches
corrected.terminal.edges<-terminal.edges-over.under

## change the termnial edges in the phylo object
cyano_tree$edge.length[tip.ids]<-corrected.terminal.edges

ape::write.tree(cyano_tree, file="ugam1_CorrCal2.71_sample_pruned.chronogram")

