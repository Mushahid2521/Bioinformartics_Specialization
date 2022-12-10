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


with open('dataset_3_2.txt', 'r') as file:
    dna = file.readlines()[0].strip()
    print(complementary_strand(dna))
