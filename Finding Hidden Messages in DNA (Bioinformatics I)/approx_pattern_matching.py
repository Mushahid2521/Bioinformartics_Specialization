def approx_pattern_matching(pattern, genome, d):
    k = len(pattern)
    pos = []

    def hamming(pt1, pt2, dis):
        for c1, c2 in zip(pt1, pt2):
            if c1 != c2:
                dis -= 1

        return True if dis >= 0 else False

    for i in range(len(genome) - k + 1):
        current_seq = genome[i:i + k]
        if hamming(pattern, current_seq, d):
            pos.append(i)

    return len(pos)



with open('dataset/dataset_9_6.txt', 'r') as file:
    lines = [s.strip() for s in file.readlines()]
    print(approx_pattern_matching(lines[0], lines[1], int(lines[2])))
