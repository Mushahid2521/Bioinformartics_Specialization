import numpy as np


def motif_score(list_strings):
    s = 0
    k = len(list_strings[0])
    t = len(list_strings)

    for i in range(k):
        mapping = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        for j in range(len(list_strings)):
            mapping[list_strings[j][i]] += 1

        for v in mapping.values():
            if v > 0:
                c = v / t
                s += (- c * np.log2(c))

    return s


print(motif_score([
    "TCGGGGGTTTTT",
    "CCGGTGACTTAC",
    "ACGGGGATTTTC",
    "TTGGGGACTTTT",
    "AAGGGGACTTCC",
    "TTGGGGACTTCC",
    "TCGGGGATTCAT",
    "TCGGGGATTCCT",
    "TAGGGGAACTAC",
    "TCGGGTATAACC"
]))
