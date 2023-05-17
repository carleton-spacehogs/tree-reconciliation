library(coda)
library(dplyr)
# df = read.delim("C:/Users/Admin Juntao/Desktop/short_df.txt")

# shell script
# grep -nP "0.333333\t0.333333\t0.333333\t0.333333\t0.333333\t0.333333\t0.333333\t0.333333" 3-location-log.txt | tail -1
# only_matrix=3-loc_only_matrix.txt
# grep "Node0 P(B)" 3-location-log.txt -m 1 > $only_matrix
# tail -n 100000 3-location-log.txt >> $only_matrix
# cut -f13-103 $only_matrix > short_df.txt

MCMC_matrix = "3-loc_only_matrix.txt"
MCM_summary = "3-loc_confidence.csv"
num_states = 3
stationary_p = 0.05

df = read.delim(MCMC_matrix)
df = df[, -c(1:13)]
df = df[grepl("Node*P*", colnames(df))]
df = df[complete.cases(df), ]

count = 0

stats = apply(df, 2, function(x) {
  print(head(x))
  count <<- count + 1
  print(count)
  obj=heidel.diag(mcmc(x), pvalue = 0)
  as.numeric(obj)
})

df.stats = as.data.frame(t(stats))
colnames(df.stats) = c("stationarity", "start_iteration", "stationarity_pval",
                       "Halfwidth_test", "mean", "halfwidth")

out = df.stats %>% 
  mutate(Halfwidth_test = ifelse(Halfwidth_test == 1, "passed", "failed")) %>%
  select(c("stationarity_pval","Halfwidth_test","mean","halfwidth"))

write.csv(out, MCM_summary, row.names=TRUE, col.names = TRUE)


node_CI = read.csv(MCM_summary)
node_CI = node_CI %>% mutate(
  min_prob = mean - halfwidth,
  max_prob = mean + halfwidth,
  adj_stationarity_pval = p.adjust(stationarity_pval, method = "BH")
  ) %>% separate(origin, sep = ".P.", into = c("BT_node", "origin"))

rownames(node_CI) = node_CI[, 1]

if (nrow(node_CI)%%num_states != 0) { # sanity check
  print("ERROR! The number of rows is not a multiple of the number of states")
}

get_node_state(rows_for_a_node){
  stationary_rows = filter(rows_for_a_node, adj_stationarity_pval > stationary_p)
  greater_than_40 = stationary_rows$min_prob > 0.4
  if (sum(greater_than_40) == 2){
    print("case 1:")
    print(stationary_rows[greater_than_40])
    origins = stationary_rows[greater_than_40]$origin
    comb_origin = paste0(origin0, collapse = "")
    return(comb_origin)
  } else if (sum(greater_than_40 == 1)) {
    print("case 2:")
    print(stationary_rows[greater_than_40])
    return(stationary_rows[greater_than_40, "origin"])
  } else {
    return("not_sure")
  }
  
}

for (i in 1:(nrow(node_CI)/num_states)) {
  print(i)
  node_rows = (i*num_states-(num_states-1)):(i*num_states)
  print(node_CI[node_rows, ])
}




