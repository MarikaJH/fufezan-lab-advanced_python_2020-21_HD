import collections
import csv

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
    organism_list = ["human", "mouse", "A_thaliana", "B_subtilis", "M_maripaludis"]

    for organism in organism_list:
        with open("species/" + organism + '.fasta') as aminoacids_sequence:
            aa = count_amino_acids(aminoacids_sequence)

        print("aminoacid counts for " + organism)
        for aminoacid in aa:
            print(aminoacid + ": " + str(aa[aminoacid]))

        a_file = open(organism + "_aa_distribution.csv", "w")
        writer = csv.writer(a_file)
        writer.writerow(['aa', 'count'])
        for key, value in aa.items():
            writer.writerow([key, value])

