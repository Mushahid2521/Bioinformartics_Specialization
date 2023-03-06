from collections import defaultdict


def composition_problem(k, text):
    return " ".join([text[i:i + k] for i in range(len(text) - k + 1)])


def path_to_genome(path):
    genome = path[0]
    for i in range(1, len(path)):
        genome += path[i][-1]
    return genome


def represent_in_graph(patterns):
    adjacent_dict = defaultdict(set)
    for item in patterns:
        adjacent_dict[item[:-1]].add(item)

    with open('output.txt', 'w') as f:
        for item in patterns:
            suffixes = adjacent_dict[item[1:]]
            if suffixes:
                f.write(f"{item}: {' '.join(suffixes)}\n")


# represent_in_graph(['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT', 'GGCAC'])

with open('./datasets/dataset_198_10.txt', 'r') as f:
    lines = f.readlines()
    patterns = list(map(str, lines[0].split()))
    represent_in_graph(patterns)
