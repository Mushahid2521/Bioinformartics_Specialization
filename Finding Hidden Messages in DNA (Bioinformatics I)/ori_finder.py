chars = "ACTG"


def complementary_strand(dna):
    rev_dna = [''] * len(dna)
    for i in range(len(dna)):
        c = dna[-i - 1]
        if c == 'A':
            rev_dna[i] = 'T'
        elif c == 'T':
            rev_dna[i] = 'A'
        elif c == 'C':
            rev_dna[i] = 'G'
        elif c == 'G':
            rev_dna[i] = 'C'

    return "".join(rev_dna)


def Neighbours(pattern, d):
    if d == 0:
        return [pattern]

    suffix_neighbours = Neighbours(pattern[1:], d - 1)
    neighbours = []
    for c in chars:
        if c != pattern[0]:
            for suffix_neigh in suffix_neighbours:
                neighbours.append(c + suffix_neigh)

    if d < len(pattern):
        suffix_neighbours = Neighbours(pattern[1:], d)
        for suffix_neigh in suffix_neighbours:
            neighbours.append(pattern[0] + suffix_neigh)

    return neighbours


def frequent_k_mers_with_mismatch(genome, k, d):
    mapping = {}
    max_cnt = 0
    frequents = []

    for i in range(len(genome) - k + 1):
        current_text = genome[i:i + k]
        neighbours = []
        for dis in range(d + 1):
            neighbours.extend(Neighbours(current_text, dis))

        for neighbour in neighbours:
            cnt1 = mapping.get(neighbour, 0)
            cnt1 += 1
            mapping[neighbour] = cnt1
            if cnt1 > max_cnt:
                max_cnt = cnt1

            cnt2 = mapping.get(complementary_strand(neighbour), 0)
            cnt2 += 1
            mapping[complementary_strand(neighbour)] = cnt2
            if cnt2 > max_cnt:
                max_cnt = cnt2

    for k, v in mapping.items():
        if v == max_cnt:
            frequents.append(k)

    return frequents


def skew(genome):
    skewness = [0]
    current_min = len(genome)
    min_positions = []
    for i, c in enumerate(genome):
        if c == 'G':
            skewness.append(skewness[-1] + 1)
        elif c == 'C':
            skewness.append(skewness[-1] - 1)
        else:
            skewness.append(skewness[-1])

        if skewness[-1] < current_min:
            current_min = skewness[-1]

    for idx, val in enumerate(skewness):
        if val == current_min:
            min_positions.append(idx)

    return min_positions


def find_ori(genome_file):
    with open(genome_file, 'r') as file:
        genome_seq = "".join(gene.strip() for gene in file.readlines()[1:])
        skew_point = skew(genome_seq)
        print(skew_point)
        ori_region = genome_seq[skew_point[0] - 250: skew_point[-1] + 250]
        possible_oris = frequent_k_mers_with_mismatch(ori_region, 9, 1)
        print(possible_oris)

# find_ori('Salmonella_enterica.txt')
# find_ori('E_coli.txt')

# print(frequent_k_mers_with_mismatch("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))
