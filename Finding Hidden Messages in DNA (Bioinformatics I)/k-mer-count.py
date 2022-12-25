def count(text, pattern):
    cnt = 0
    k = len(pattern)
    for i in range(len(text) - k + 1):
        if text[i: i + k] == pattern:
            cnt += 1

    return cnt


with open('dataset/dataset_2_6.txt', 'r') as file:
    lines = file.readlines()
    print(lines[0])
    print(count(lines[0].strip(), lines[1].strip()))
