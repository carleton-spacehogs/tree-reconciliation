import sys
from Bio import Phylo

def extract_genomes(tree_file, output_file):
    # Read the tree file
    tree = Phylo.read(tree_file, 'newick')

    # Extract all leaf node names
    genome_names = [clade.name for clade in tree.find_clades() if clade.is_terminal()]

    # Write to output file
    with open(output_file, 'w') as file:
        for name in genome_names:
            file.write(name + '\n')

if __name__ == "__main__":
    tree_file = sys.argv[1]
    output_file = sys.argv[2]
    extract_genomes(tree_file, output_file)
