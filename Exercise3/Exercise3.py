import plotly
from collections import deque
import pandas as pd


def get_as_sequence_list():
    """
    Generates a list containing the 1-letter code of the aminoacid sequence of the gpcr

    Returns: aas_list as list containing the one letter code of the aminoacid sequence

    """
    amino_acid_sequence = ""
    with open("uniprot_gpcr.fasta") as gpcr_fasta:
        for line in gpcr_fasta:
            if not line.startswith('>'):
                amino_acid_sequence += line.replace('\n', '')
    aas = list(amino_acid_sequence)
    return aas


def get_hydropathy_dictionary():
    aa_properties = pd.read_csv('../data/amino_acid_properties.csv')
    aa_properties = pd.DataFrame.to_dict(aa_properties)
    aa_1_letter_code_list = aa_properties['1-letter code'].values()
    aa_hydropathy_list = aa_properties['hydropathy index (Kyte-Doolittle method)'].values()
    aa_hydropathy = dict(zip(aa_1_letter_code_list, aa_hydropathy_list))
    return aa_hydropathy


def sequence_hydropathy_list(sequence, hydropathy_dict, len_window=1):
    seq_hydropathy_list = [hydropathy_dict[aa] for aa in sequence]
    seq_hydropathy_mean = []
    window = deque([], maxlen=len_window)

    for aa in seq_hydropathy_list:
        window.append(aa)
        window_mean = sum(window) / len(window)
        seq_hydropathy_mean.append(window_mean)
    seq_hydropathy_list = seq_hydropathy_mean
    return seq_hydropathy_list


def plot_hydropathy_histogram(sequence, hydropathy, len_windows=1):
    seq_pos = []
    [seq_pos.append(pos) for pos in range(0, len(sequence))]
    graph = [plotly.graph_objs.Bar(x=seq_pos, y=hydropathy)]

    fig = plotly.graph_objs.Figure(data=graph, layout={"title": {
        "text": "Aminoacid Hydropathy For G-protein Coupled Receptor 183 With Sliding Window Of " + str(
            len_windows)}, })
    fig.update_layout(xaxis=dict(title='Sequence position'), yaxis=dict(title='Aminoacid Hydropathy'))
    fig.show()


if __name__ == "__main__":
    aa_sequence = get_as_sequence_list()
    hydropathy_dict = get_hydropathy_dictionary()
    for window_lenght in [1, 5, 10]:
        aas_hydropathy = sequence_hydropathy_list(aa_sequence, hydropathy_dict, window_lenght)
        plot_hydropathy_histogram(aa_sequence, aas_hydropathy, window_lenght)
