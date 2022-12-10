def frequentWords(text, k):
    kmers = {}
    highest_kmers = []
    maximum_cnt = 0
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        if kmer in kmers:
            kmers[kmer] += 1
            if kmers[kmer] > maximum_cnt:
                maximum_cnt = kmers[kmer]
        else:
            kmers[kmer] = 1
            if kmers[kmer] > maximum_cnt:
                maximum_cnt = kmers[kmer]

    for key in kmers.keys():
        if kmers[key] == maximum_cnt:
            highest_kmers.append(key)

    return highest_kmers


# with open('dataset_2_13.txt', 'r') as file:
#     lines = file.readlines()
#     print(" ".join(frequentWords(lines[0].strip(), int(lines[1].strip()))))


def position_of_kmer(dna, pattern):
    pos = []
    k = len(pattern)
    for i in range(len(dna) - len(pattern) + 1):
        if dna[i:i + k] == pattern:
            pos.append(str(i))

    return " ".join(pos)


