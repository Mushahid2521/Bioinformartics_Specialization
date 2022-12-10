from utils import Neighbours
from utils import is_hamming, hamming_distance
from itertools import product
from collections import Counter


def motif_search_brute_force(dna, k, d):
    # k = 3
    # d = 1
    # dna = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']

    kmers = set([])
    result = set([])
    dna_concat = "".join(dna)

    for i in range(len(dna_concat) - k + 1):
        this_mer = dna_concat[i: i + k]
        neighbours = Neighbours(this_mer, d)

        for this_neighbour in neighbours:
            if this_neighbour not in kmers:
                is_found = True
                for this_dna in dna:
                    res = any([is_hamming(this_dna[j:j + k], this_neighbour, d) for j in range(len(this_dna) - k + 1)])
                    if not res:
                        is_found = False
                        break

                if is_found:
                    result.add(this_neighbour)

                kmers.add(this_neighbour)

    return result


def medianStringMotif(dna, k):
    possible_medians = product(['A', 'T', 'C', 'G'], repeat=k)
    possible_medians = ["".join(p) for p in possible_medians]

    min_distance = float('inf')
    median_string = ""
    for string in possible_medians:
        this_string_sum = 0
        for this_dna in dna:
            min_hamming = float('inf')
            for i in range(len(this_dna) - k + 1):
                dis = hamming_distance(this_dna[i:i + k], string)
                if min_hamming > dis:
                    min_hamming = dis
            this_string_sum += min_hamming
        if this_string_sum < min_distance:
            min_distance = this_string_sum
            median_string = string

    return median_string


def scoreMotif(motifs):
    score = 0
    for i in range(len(motifs[0])):
        this_list = [motifs[j][i] for j in range(len(motifs))]
        counts = dict(Counter(this_list))
        score += sum(counts.values()) - max(counts.values())

    return score


def create_profile(motifs):
    motif_profile = {'A': [], 'C': [], 'G': [], 'T': []}
    for i in range(len(motifs[0])):
        this_list = [motifs[j][i] for j in range(len(motifs))]
        counts = dict(Counter(this_list))
        for p in motif_profile.keys():
            motif_profile[p].append(counts.get(p, 0) / len(motifs))

    return motif_profile


def profile_most_probable_kmer(text, k, profile):
    max_proba = -1
    result_str = None
    for i in range(len(text) - k + 1):
        this_proba = 1
        this_mer = text[i:i + k]
        for idx, c in enumerate(this_mer):
            this_proba *= profile[c][idx]

        if this_proba > max_proba:
            max_proba = this_proba
            result_str = this_mer

    return result_str


def greedy_motif_search(dna, k, t):
    best_motifs = [dna[i][:k] for i in range(len(dna))]
    for i in range(len(dna[0]) - k + 1):
        this_mer = dna[0][i:i + k]
        motif = [this_mer]
        for j in range(1, len(dna)):
            profile_this = create_profile(motif)
            next_motif = profile_most_probable_kmer(dna[j], k, profile_this)
            motif.append(next_motif)

        if scoreMotif(motif) < scoreMotif(best_motifs):
            best_motifs = motif

    return best_motifs


with open('dataset_159_5.txt', 'r') as f:
    lines = f.readlines()
    k, t = map(int, lines[0].split())
    dna = list(map(str, lines[1].split()))
    print(" ".join(greedy_motif_search(dna, k, t)))

# if __name__ == "__main__":
# Brute Force Motif Search
# with open('dataset_156_8.txt', 'r') as f:
#     lines = f.readlines()
#     k, d = map(int, lines[0].split())
#     input_dna = lines[1].split()
#     print(" ".join(motif_search_brute_force(input_dna, k, d)))

# Consensus string search
# with open('dataset_158_9.txt', 'r') as f:
#     lines = f.readlines()
#     k = int(lines[0])
#     input_dna = []
#     for l in range(1, len(lines)):
#         input_dna.append(lines[l].strip())
#
#     print(medianStringMotif(input_dna, k))

# Profile Most Probable k mer
# with open('dataset_159_3.txt', 'r') as f:
#     lines = f.readlines()
#     lines = [line.strip() for line in lines]
#     text = str(lines[0])
#     k = int(lines[1])
#     profile = {'A': list(map(float, lines[2].split())),
#                'C': list(map(float, lines[3].split())),
#                'G': list(map(float, lines[4].split())),
#                'T': list(map(float, lines[5].split()))}
#
#     print(profile_most_probable_kmer(text, k, profile))
