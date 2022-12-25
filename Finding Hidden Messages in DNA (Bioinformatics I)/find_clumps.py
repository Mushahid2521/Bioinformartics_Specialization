def frequentWords(text, k, t):
    valid_kmers = []
    kmers = {}
    kmer_pos = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]

        if kmer in kmers:
            kmers[kmer] += 1
            pos = kmer_pos.get(kmer, [])
            pos.append(i)
            kmer_pos[kmer] = pos
            if kmers[kmer] >= t:
                valid_kmers.append(kmer)

        else:
            kmers[kmer] = 1
            pos = kmer_pos.get(kmer, [])
            pos.append(i)
            kmer_pos[kmer] = pos
            if kmers[kmer] >= t:
                valid_kmers.append(kmer)

    return valid_kmers, kmer_pos


def findClumps(dna, K, L, t):
    patterns = set([])
    valid_kmers, kmer_pos = frequentWords(dna, K, t)
    print(valid_kmers)
    for kmer in valid_kmers:
        print(kmer_pos[kmer])

    return patterns


with open('dataset/E_coli.txt', 'r') as file:
    lines = file.readlines()
    dna = lines[0].strip()
    K, L, t = map(int, lines[1].strip().split(" "))
    print(" ".join(findClumps("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", 5, 50, 4)))

