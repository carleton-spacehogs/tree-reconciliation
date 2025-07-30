import sys
from Bio import Phylo
import rpy2.robjects as robjects
import tempfile
import os

def remove_bootstrap_and_write_nexus(input_file, output_file):
    r_script = """
    library(ape)
    tree <- read.tree("{input_file}")
    tree$node.label <- NULL 
    write.nexus(tree, file = "{output_file}")
    """.format(input_file=input_file, output_file=output_file)

    robjects.r(r_script)

def translate_tree_names(tree):
    translation_map = {}
    counter = 1
    for clade in tree.find_clades():
        if clade.is_terminal():  # Check if it's a leaf node
            original_name = clade.name
            clade.name = str(counter)
            translation_map[str(counter)] = original_name
            counter += 1
    return tree, translation_map

def convert_tree_file(input_tree_file, output_nexus_file):
    tree = Phylo.read(input_tree_file, 'newick')
    tree, translation_map = translate_tree_names(tree)

    # Create a temporary Newick file with translated names
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmpfile:
        Phylo.write(tree, tmpfile, 'newick')

    # Use R to remove bootstrap values and write Nexus
    remove_bootstrap_and_write_nexus(tmpfile.name, output_nexus_file)

    # Clean up the temporary file
    os.remove(tmpfile.name)

if __name__ == "__main__":
    input_tree_file = sys.argv[1]
    output_nexus_file = sys.argv[2]
    convert_tree_file(input_tree_file, output_nexus_file)
