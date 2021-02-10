import plotly
import collections
from Exercise2 import count_aas as c_aas


def plot_histogram(data, title):
    """
    Plots a histogram of counts and aminoacids for the data in a given dictionary with the given title
    :param data: dictionary
    :param title: title of the plot
    :return:
    """

    amino_acids = list(data.keys())
    amino_acids_counts = list(data.values())

    graph = [plotly.graph_objs.Bar(x=amino_acids, y=amino_acids_counts)]
    fig = plotly.graph_objs.Figure(data=graph, layout={"title": {"text": title}})
    fig.update_layout(xaxis=dict(title='aminoacid'), yaxis=dict(title='counts'), template="plotly")
    fig.show()


if __name__ == "__main__":
    organism_list = ["human", "mouse", "A_thaliana", "B_subtilis", "M_maripaludis"]
    for organism in organism_list:
        with open("species/" + organism + '.fasta') as aminoacids_sequence:
            aa = c_aas.count_amino_acids(aminoacids_sequence)
            sorted_aa = collections.OrderedDict(sorted(aa.items()))
            plot_histogram(sorted_aa, 'aminoacid distribution of ' + organism)
