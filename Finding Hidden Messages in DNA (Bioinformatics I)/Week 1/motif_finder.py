from utils import Neighbours
from utils import is_hamming, hamming_distance
from itertools import product


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


if __name__ == "__main__":
    # Brute Force Motif Search
    # with open('dataset_156_8.txt', 'r') as f:
    #     lines = f.readlines()
    #     k, d = map(int, lines[0].split())
    #     input_dna = lines[1].split()
    #     print(" ".join(motif_search_brute_force(input_dna, k, d)))

    # Consensus string search
    with open('dataset_158_9.txt', 'r') as f:
        lines = f.readlines()
        k = int(lines[0])
        input_dna = []
        for l in range(1, len(lines)):
            input_dna.append(lines[l].strip())

        print(medianStringMotif(input_dna, k))
