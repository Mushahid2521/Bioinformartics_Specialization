import random

import numpy
from numpy.random import choice
from utils import Neighbours
from utils import is_hamming, hamming_distance
from itertools import product
from collections import Counter


def consensus_string(motifs):
    res_string = ""
    for j in range(len(motifs[0])):
        chars = [motifs[i][j] for i in range(len(motifs))]
        res_string += max(chars, key=chars.count)
    return res_string


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


def distance_between_pattern_string(pattern, dna):
    distance = 0
    for text in dna:
        this_distance = float('inf')
        for i in range(len(text) - len(pattern) + 1):
            this_mer = text[i:i + len(pattern)]
            dis = hamming_distance(pattern, this_mer)
            if dis < this_distance:
                this_distance = dis

        distance += this_distance

    return distance


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
        dict_counts = dict(Counter(this_list))
        for p in motif_profile.keys():
            motif_profile[p].append((dict_counts.get(p, 0) + 1) / (len(motifs) + 4))

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


def profile_randomly_generated_kmer(text, k, profile):
    proba_kmers = []
    for i in range(len(text) - k + 1):
        this_mer = text[i:i + k]
        proba_kmers.append(sum([profile[c][idx] for idx, c in enumerate(this_mer)]))

    proba_kmers = [p / sum(proba_kmers) for p in proba_kmers]
    idx_choice = choice([i for i in range(len(text) - k + 1)], 1, p=proba_kmers)[0]
    return text[idx_choice: idx_choice + k]


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


def randomizedMotifSearch(dna, k, t):
    ultimate_motif = None
    score = 0
    for _ in range(1000):
        best_motifs = []
        for each in dna:
            this_i = random.randint(0, len(each) - k)
            best_motifs.append(each[this_i:this_i + k])
        while True:
            profile = create_profile(best_motifs)
            motifs = [profile_most_probable_kmer(each, k, profile) for each in dna]
            if scoreMotif(motifs) < scoreMotif(best_motifs):
                best_motifs = motifs
            else:
                break

            if ultimate_motif is None or scoreMotif(ultimate_motif) > scoreMotif(best_motifs):
                ultimate_motif = best_motifs
                score = scoreMotif(ultimate_motif)
    print(score)
    return ultimate_motif


def GibbsSampler(dna, k, t, N):
    best_motifs = randomizedMotifSearch(dna, k, t)
    score = 0
    for j in range(N):
        i = random.randint(0, t - 1)
        profile = create_profile([m for p, m in enumerate(best_motifs) if p != i])
        motif_i = profile_randomly_generated_kmer(dna[i], k, profile)
        motifs = [profile_most_probable_kmer(text, k, profile) for text in dna]
        motifs[i] = motif_i

        if scoreMotif(motifs) < scoreMotif(best_motifs):
            best_motifs = motifs
            score = scoreMotif(best_motifs)

    print(score)
    return best_motifs


# MTB - Mycobacterium Tuberculosis, DosR - Dormancy Survival Regulator transcription factor binding site
with open('dataset/upstream250_MTB_genes.txt', 'r') as f:
    dna = [line.strip() for idx, line in enumerate(f.readlines()) if idx % 2 != 0]
    # # print("Median string with K= 12\n", consensus_string(medianStringMotif(dna, 12)))
    # print("Randomized Motif Search, k=12", consensus_string(randomizedMotifSearch(dna, 12, len(dna))))
    # print("Gibbs Sampler, k=12", consensus_string(GibbsSampler(dna, 12, len(dna), 1000)))
    # print("Randomized Motif Search, k=20", consensus_string(randomizedMotifSearch(dna, 20, len(dna))))
