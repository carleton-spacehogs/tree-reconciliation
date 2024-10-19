import sys

def main(blast_file, orf_file):
    try:
        # Reading the BLAST file and storing ORF names
        with open(blast_file, 'r') as infile:
            blast_orf_names = [line.split('\t')[0] for line in infile]
        
        # Output file named after the ORF file with '_BLASThits.txt' suffix
        output_file = f"{orf_file}_BLASThitstttt.txt"
        
        # Process the ORF coverage file and output matching lines
        with open(orf_file, 'r') as infile, open(output_file, 'w') as outfile:
            next(infile)  # Skip header if there is one
            for line in infile:
                cols = line.split('\t')
                orf_title = cols[4]  # Change index here if the index is cols[0]
                if orf_title in blast_orf_names:
                    outfile.write(line)
        
        print(f"Matching sequences written to {output_file}")
    
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_ORF_covg_from_BLAST_hits.py [BLAST file] [ORF coverage file]")
        sys.exit(1)
    
    blast_file = sys.argv[1]
    orf_file = sys.argv[2]
    main(blast_file, orf_file)
