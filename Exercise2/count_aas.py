import collections
import csv
import argparse

def count_amino_acids(fasta_file):
    """
    Counts the number of occurences of aminoacids in an opened fasta-file.

    :param fasta_file: textfile containing sequence or protein(s)
    :return: counter as dictionary containing the aminoacids with respective counts
    """
    amino_acid_sequence = []
    for line in fasta_file:
        if not line.startswith('>'):
            amino_acid_sequence += line.replace('\n', '')
    counter = collections.Counter(amino_acid_sequence)
    return counter


if __name__ == "__main__":
    # example command for terminal:
    # python Exercise2\count_aas.py --fasta_path Exercise2\species\human.fasta --csv_path .\Exercise2\human_aa_distribution.csv

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_path", help="path of the FASTA file with the sequences to be analyzed", type=str)
    parser.add_argument("--csv_path", help="path for the csv file in which the output will be stored", type=str)
    arguments = parser.parse_args()
    fasta_file_path = arguments.fasta_path
    csv_path = arguments.csv_path

    with open(fasta_file_path) as aminoacids_sequence:
        aa = count_amino_acids(aminoacids_sequence)
    aa_sorted = sorted(aa.items())

    print("aminoacid counts for fasta file")
    for aminoacid in aa_sorted:
        print(aminoacid[0] + ": " + str(aminoacid[1]))

    a_file = open(csv_path, "w")
    writer = csv.writer(a_file)
    writer.writerow(['aa', 'count'])
    for aa_count in aa_sorted:
        writer.writerow([aa_count[0], aa_count[1]])
