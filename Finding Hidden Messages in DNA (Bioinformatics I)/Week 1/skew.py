def hamming_distance(input1, input2):
    dis = 0
    for i in range(len(input1)):
        if input2[i] != input1[i]:
            dis += 1

    return dis


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


#
#
# with open('dataset_7_10.txt', 'r') as file:
#     genome = file.readlines()[0]
#     print(" ".join(map(str, skew(genome))))

with open('dataset_9_3.txt', 'r') as file:
    lines = file.readlines()
    print(hamming_distance(lines[0].strip(), lines[1].strip()))
