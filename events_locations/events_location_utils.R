state_vec_short = c("B", "C", "D")
state_vec_long = c("marine_deep","marine_shallow","terrestrial")
state_map2 = c("B" = "1,0,0", "C" = "0,1,0", "D" = "0,0,1")
src_ecceTERA = sprintf("../%s_ecceTERA_analysis", clock_model)

state_map = c("A" = "host_associated", "B" = "marine_deep", "C" = "marine_shallow", "D" = "terrestrial", "E" = "ungrouped")

init_env = function() {
  library(reshape2)
  library(stringr)
  library(tidyverse)
}

sort_words = function(str) {
  vec = strsplit(str,split = ' ')[[1]]
  return(paste(sort(vec), collapse = " "))
}

get_prob_matrix = function() {
  location_prob = as.data.frame(t(read.delim("node-location-probability.txt", header=FALSE)))
  colnames(location_prob) = c("node_state", "probability")
  
  # don't infer anything about the root
  first = which(location_prob$node_state == sprintf("Node0 P(%s)", state_vec_short[1]))
  
  # last 2 rows are summary statistics
  last = nrow(location_prob) - 2
  
  location_prob = location_prob[first:last, ]
  location_prob$node_state = str_replace(location_prob$node_state, "Node", "")
  location_prob = location_prob %>%
    mutate(node_state = str_replace(node_state, "Node", "")) %>%
    separate(node_state, sep = " P\\(", into = c("BT_node", "location")) %>% # BT for Bayes Traits
    mutate(location = str_replace(location, "\\)", ""),
           BT_node = as.numeric(BT_node))
  
  location_prob$loc = str_replace_all(location_prob$location, state_map)
  prob_matrix = dcast(location_prob, BT_node ~ loc, value.var = "probability")
  return(prob_matrix)
}

get_BT_coding = function() {
  bayesTraits_tree = read.delim("bayesTraits_species_tree_mapping.txt", header=FALSE,
                                col.names = c("BT_node", "leaves_count", "leaves_below"))
  bayesTraits_tree = bayesTraits_tree %>%
    mutate(BT_node = str_replace(BT_node, "Node", ""),
           leaves_below = lapply(leaves_below, sort_words))
  return(bayesTraits_tree)
}


