#/bin/python3
from Bio import Phylo
import csv
tree = Phylo.read("ecceTERA_ugam1_ChenParamsEarth_species_tree.xml", "phyloxml")

# tree.root.__dict__
# {'branch_length': None, 'id_source': None, 'name': '1730', 'width': None, '_color': None, 'node_id': None, 'events': None, 'binary_characters': None, 'date': None, 'confidences': [], 'taxonomies': [], 'sequences': [], 'distributions': [], 'references': [], 'properties': [], 'clades': [Clade(name='1729'), Clade(name='UNSAMPLED')], 'other': []}

internal_node_list = []
root_leaf_list = []

def get_all_internal_nodes(node, list):
	if node.name.isnumeric():
		internal_node_list.append(node)
		get_all_internal_nodes(node.clades[0], list)
		get_all_internal_nodes(node.clades[1], list)

get_all_internal_nodes(tree.root, internal_node_list)

def get_all_leaf(node, list):
	if node.name.isnumeric():
		get_all_leaf(node.clades[0], list)
		get_all_leaf(node.clades[1], list)
	else:
		list.append(node.name)

all_node_dict = []
for internal_node in internal_node_list:
	print(internal_node.name)
	leaf_list = []
	get_all_leaf(internal_node, leaf_list)
	leaf_list.sort()
	all_node_dict.append([internal_node.name, ' '.join(leaf_list)])

with open("leaves_of_ecceTERA_interal_nodes.csv", "w") as outfile:
	writer = csv.writer(outfile)
	writer.writerow(["ecceTERA_node","leaves_below"])
	writer.writerows(all_node_dict)
