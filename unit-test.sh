#!/bin/bash

# COG5403,30
gene=COG5403

source ./scripts/declare_file_location.sh
source ./scripts/declare_file_location.sh $gene

./reconcile.sh --stop_before_tree_building --COG $gene
echo --COG ==== done

./reconcile.sh --stop_before_tree_building --gene_sequences tmp/$gene.faa
echo --gene_sequences ==== done

./reconcile.sh --stop_before_reconciliation --alignment $alignment
echo --alignment ==== done

./reconcile.sh --gene_tree $iqtree_ufboot
echo --gene_tree $iqtree_ufboot ==== done

# # start over and do it again as a whole pipeline
rm $trimv2 $iqtree_ufboot $store_gene_tree_filename $ecceTERA_sym $chronogram_internal_nodes_f $sym_event_f $sym_event_date_f $R_plot

./reconcile.sh --COG $gene
echo full run --COG ==== done

