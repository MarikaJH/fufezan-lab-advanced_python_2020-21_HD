import collections
import csv
import matplotlib.pyplot as plt


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


def plot_histogram(data, title):
    """
    Plots a histogram of counts and aminoacids for the data in a given dictionary with the given title
    :param data: dictionary
    :param title: title of the plot
    :return:
    """
    amino_acids = data.keys()
    amino_acids_counts = data.values()

    plt.bar(amino_acids, amino_acids_counts)
    plt.title(title)
    plt.xlabel('aminoacid')
    plt.ylabel('count')
    plt.show()


if __name__ == "__main__":
    organism_list = ["human", "mouse", "A_thaliana", "B_subtilis", "M_maripaludis"]

    for organism in organism_list:
        with open("species/" + organism + '.fasta') as aminoacids_sequence:
            aa = count_amino_acids(aminoacids_sequence)

        print("aminoacid counts for " + organism)
        for aminoacid in aa:
            print(aminoacid + ": " + str(aa[aminoacid]))

        plot_histogram(aa, 'aminoacid distribution of ' + organism)

        a_file = open(organism + "_aa_distribution.csv", "w")
        writer = csv.writer(a_file)
        writer.writerow(['aa', 'count'])
        for key, value in aa.items():
            writer.writerow([key, value])
