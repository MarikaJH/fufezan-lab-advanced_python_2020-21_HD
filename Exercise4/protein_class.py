import plotly.io as pio
from collections import deque
import pandas as pd
import requests
from Exercise3 import Exercise3 as Ex3
pio.renderers.default = "browser"


class Protein(object):
    def __init__(self, protein_id: str):
        """

        Args:
            protein_id: String containing the uniprot-id of a protein.
        """
        self.protein_id = protein_id

    def get_data(self):
        """
        Gets the sequence corresponding to the protein-id from uniprot and returns it as a list.

        Returns:
            List containing the protein sequence as 1-letter-codes.

        """
        url = 'https://www.uniprot.org/uniprot/' + self.protein_id + '.fasta?fil=reviewed:yes'
        r = requests.get(url)
        protein_fasta = '../data/' + self.protein_id + '.fasta'

        aminoacid_sequence = ""
        with open(protein_fasta, 'wb') as file:
            file.write(r.content)
            file.close()
        with open(protein_fasta) as file:
            for line in file:
                if not line.startswith('>'):
                    aminoacid_sequence += line.replace('\n', '')
        return list(aminoacid_sequence)

    def map(self, aa_lookup: dict, aa_property: str, len_window=1):
        """
        Gets the protein sequence using the get_data function.
        Maps the sequence with the given property from the given lookup using a sliding window of the given length.

        Args:
            aa_lookup: Nested dictionary with properties as keys and dictionaries containing the property values for
             each aminoacid.
            aa_property: String containing the name of the property to be mapped. must match one of the keys in
             aa_lookup.
            len_window: Integer giving the length of the sliding average window, default value is 1.

        Returns:
            List contaiing the property-calue for each aminoacid in the proteins' sequence.

        """
        sequence = self.get_data()
        property_dict = aa_lookup[aa_property]

        seq_property_list = [property_dict[aa] for aa in sequence]
        window = deque([], maxlen=len_window)

        if len_window > 1:
            seq_property_mean = []
            for aa in seq_property_list:
                window.append(aa)
                window_mean = sum(window) / len(window)
                seq_property_mean.append(round(window_mean, 1))
            seq_property_list = seq_property_mean
        return seq_property_list


def get_aa_properties_lookup_dict():
    """
    Gets the data for the aminoacid properties and puts them in a nested dictionary.

    Returns:
        nested dictionary with properties as keys and dictionaries containing the property values for aminoacids.

    """
    aa_properties = pd.read_csv('../data/amino_acid_properties.csv')
    aa_properties_dict = pd.DataFrame.to_dict(aa_properties)

    aa_1_letter_code_list = aa_properties_dict['1-letter code'].values()
    aa_properties_lookup = {}
    for pos, aa_property in enumerate(aa_properties.keys()):
        if pos > 2:
            property_list = aa_properties_dict[aa_property].values()
            property_dict = dict(zip(aa_1_letter_code_list, property_list))
            aa_properties_lookup[aa_property] = property_dict
    return aa_properties_lookup


if __name__ == "__main__":
    lookupdict = get_aa_properties_lookup_dict()
    protein = Protein("P32249")
    window_len = 10
    mapped_sequence = protein.map(lookupdict, "hydropathy index (Kyte-Doolittle method)", window_len)
    print(mapped_sequence)
    Ex3.plot_hydropathy_histogram(protein.get_data(), mapped_sequence, window_len)
