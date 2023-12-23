def method_gauss(A, b):
    size1 = len(A)
    for i in range(size1):
        index_tmp = i
        for j in range(i + 1, size1):
            if abs(A[j][i]) > abs(A[index_tmp][i]):
                index_tmp = j

        A[i], A[index_tmp] = A[index_tmp], A[i]
        b[i], b[index_tmp] = b[index_tmp], b[i]

        for j in range(i + 1, size1):
            ratio = A[j][i] / A[i][i]
            for k in range(i, size1):
                A[j][k] -= ratio * A[i][k]
            b[j] -= ratio * b[i]
    x = [0] * size1
    for i in range(size1 - 1, -1, -1):
        x[i] = b[i] / A[i][i]
        for j in range(i + 1, size1):
            x[i] -= A[i][j] * x[j] / A[i][i]
    return x


def method_yacobi(A, b, eps):
    size = len(b)
    c = [0 for _ in range(size)]
    n = 0

    while True:
        xn1 = c.copy()
        for i in range(size):
            sum = 0
            for j in range(size):
                if j != i:
                    sum += A[i][j] * xn1[j]
            c[i] = (b[i] - sum) / A[i][i]
        n += 1

        if max([abs(c[i] - xn1[i]) for i in range(size)]) <= eps:
            break

    return c, n

def method_gauss_zeydel(A, b, eps):
    t = len(b)
    c = [0 for _ in range(t)]
    n = 0

    while True:
        xn = c.copy()
        for i in range(t):
            sum = 0
            for j in range(t):
                if j != i:
                    sum += A[i][j] * c[j] # тут уже обновленный на прошлой итерации c
            c[i] = (b[i] - sum) / A[i][i]
        n += 1

        if max([abs(c[i] - xn[i]) for i in range(t)]) <= eps:
            break

    return c, n

eps = 0.00005
A = [[-0.19, -0.31, -0.51], [-0.11, 0.89, 0.22], [1.61, -0.09, 0.2]]
b = [-2.34, 2.33, 2.03]
print("Метод Гаусса:", tuple(method_gauss(A, b)), sep='\n')
x, n = method_yacobi(A, b, eps=eps)
print('Якоби:', x, 'n:', n)

x, n = method_gauss_zeydel(A, b, eps=eps)
print('Гаусса-Зейдель:', x, 'n:', n)