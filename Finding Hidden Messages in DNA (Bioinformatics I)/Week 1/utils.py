chars = "ACTG"


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


def is_hamming(pt1, pt2, dis):
    for c1, c2 in zip(pt1, pt2):
        if c1 != c2:
            dis -= 1

    return True if dis >= 0 else False


def hamming_distance(pt1, pt2):
    dis = 0
    for c1, c2 in zip(pt1, pt2):
        if c1 != c2:
            dis += 1

    return dis
