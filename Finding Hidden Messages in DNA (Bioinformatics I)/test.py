def valid(p):
    cnt = 0
    for c in p:
        if c == '(':
            cnt += 1
        elif c == ')':
            if cnt < 1:
                return False
            cnt -= 1

    return True if cnt == 0 else False


def generate(A, n, ans):
    if len(A) == 2 * n:
        if valid(A):
            ans.append("".join(A))

    else:
        A.append('(')
        generate(A, n, ans)
        A.pop()
        A.append(')')
        generate(A, n, ans)
        A.pop()


ans = []
generate([], 3, ans)
print(ans)
