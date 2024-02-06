import sys
import pandas as pd

def truncate_accession(accession):
    # Convert to string in case the accession is not in string format
    accession = str(accession)
    return accession[4:] if len(accession) > 4 else accession

def convert_predictions(value):
    # Convert prediction values '2.0', '1.0', and '0.0' to 'c', 'b', and 'a'
    if value == 2.0:
        return 'C'
    elif value == 1.0:
        return 'B'
    elif value == 0.0:
        return 'A'
    else:
        return '-'

def merge_csv(genome_csv, predictions_csv, output_csv):
    # Load the genome names
    genomes = pd.read_csv(genome_csv, header=None, names=["ncbi_genbank_assembly_accession"])
    genomes['truncated_accession'] = genomes['ncbi_genbank_assembly_accession'].apply(truncate_accession)

    # Load the predictions
    predictions = pd.read_csv(predictions_csv)
    predictions['truncated_accession'] = predictions['ncbi_genbank_assembly_accession'].apply(truncate_accession)

    # Merge the dataframes using the truncated accession numbers
    merged = pd.merge(genomes, predictions[['truncated_accession', 'Predictions']], left_on='truncated_accession', right_on='truncated_accession', how='left')

    # Drop the temporary 'truncated_accession' column
    merged = merged.drop('truncated_accession', axis=1)

    # Apply the conversion function to the 'Predictions' column
    merged['Predictions'] = merged['Predictions']   .apply(convert_predictions)

    # Save to output CSV without header
    merged.to_csv(output_csv, sep ="\t",index=False, header=False)

if __name__ == "__main__":
    genome_file = sys.argv[1]
    predictions_file = sys.argv[2]
    output_file = sys.argv[3]
    merge_csv(genome_file, predictions_file, output_file)
