print("\n ---------------1-------------------- \n")
def neville(x, y, target):
    n = len(x)
    Q = [[0] * n for i in range(n)]
    for i in range(n):
        Q[i][0] = y[i]
    for j in range(1, n):
        for i in range(n - j):
            Q[i][j] = ((target - x[i+j]) * Q[i][j-1] + (x[i] - target) * Q[i+1][j-1]) / (x[i] - x[i+j])
    return Q[0][-1]

x = [3.6, 3.8, 3.9]
y = [1.675, 1.436, 1.318]
target = 3.7
result = neville(x, y, target)

print("The 2nd degree interpolating value for f(3.7) is: ", result)
print("\n --------------2--------------------- \n")

def newton_forward(x, y):
    n = len(x)
    diff_table = [[0] * (n-i) for i in range(n)]
    for i in range(n):
        diff_table[i][0] = y[i]
    for j in range(1, n):
        for i in range(n-j):
            diff_table[i][j] = diff_table[i+1][j-1] - diff_table[i][j-1]
    return diff_table

def compute_coefficients(x, diff_table, degree):
    coefficient = diff_table[0][0]
    term = 1
    for i in range(1, degree+1):
        term = term * (x - i + 1) / i
        coefficient += term * diff_table[0][i]
    return coefficient

x = [7.2, 7.4, 7.5, 7.6]
y = [23.5492, 25.3913, 26.8224, 27.4589]
diff_table = newton_forward(x, y)

for degree in range(1, 4):
    print(f"Degree {degree} polynomial approximation:")
    print("P(x) = ", end="")
    for i in range(degree+1):
        coefficient = compute_coefficients(x[0], diff_table, i)
        if i > 0:
            print(" + ", end="")
        print(f"{coefficient:.4f}", end="")
        for j in range(i):
            print(f"(x - {x[j]:.1f})", end="")
    print("\n")

print("\n ---------------3-------------------- \n")

def newton_forward(x, y):
    n = len(x)
    diff_table = [[0] * (n-i) for i in range(n)]
    for i in range(n):
        diff_table[i][0] = y[i]
    for j in range(1, n):
        for i in range(n-j):
            diff_table[i][j] = diff_table[i+1][j-1] - diff_table[i][j-1]
    return diff_table

def compute_coefficients(x, diff_table, degree):
    coefficient = diff_table[0][0]
    term = 1
    for i in range(1, degree+1):
        term = term * (x - i + 1) / i
        coefficient += term * diff_table[0][i]
    return coefficient

x = [7.2, 7.4, 7.5, 7.6]
y = [23.5492, 25.3913, 26.8224, 27.4589]
diff_table = newton_forward(x, y)

print("Degree 1 polynomial approximation:")
degree = 1
P = compute_coefficients(x[0], diff_table, degree)
print(f"P({7.3}) = {P:.4f}")

print("\nDegree 2 polynomial approximation:")
degree = 2
P = compute_coefficients(x[0], diff_table, degree)
P += compute_coefficients(x[1], diff_table, degree-1) * (7.3 - x[0])
print(f"P({7.3}) = {P:.4f}")

print("\nDegree 3 polynomial approximation:")
degree = 3
P = compute_coefficients(x[0], diff_table, degree)
P += compute_coefficients(x[1], diff_table, degree-1) * (7.3 - x[0])
P += compute_coefficients(x[2], diff_table, degree-2) * (7.3 - x[0]) * (7.3 - x[1])
print(f"P({7.3}) = {P:.4f}")

print("\n -----------------4------------------ \n")

import numpy as np

x = np.array([3.6, 3.8, 3.9])
y = np.array([1.675, 1.436, 1.318])
y_prime = np.array([-1.195, -1.188, -1.182])

# Construct the divided difference table
n = len(x)
table = np.zeros((2*n, 2*n))
table[::2, 0] = x
table[1::2, 0] = x
table[::2, 1] = y
table[1::2, 1] = y_prime

for j in range(2, 2*n):
    for i in range(2*n-j):
        if table[i, 0] == table[i+j-1, 0]:
            table[i, j] = table[i+1, j-1] / np.math.factorial(j-1)
        else:
            table[i, j] = (table[i+1, j-1] - table[i, j-1]) / (table[i+j-1, 0] - table[i, 0])

# Print the Hermite polynomial approximation matrix
print("Hermite polynomial approximation matrix:")
print(table)

print("\n ----------------5------------------- \n")


from scipy import interpolate

x = np.array([2, 5, 8, 10])
y = np.array([3, 5, 7, 9])

# Use cubic spline interpolation to find the coefficients of the cubic polynomials
spline = interpolate.CubicSpline(x, y, bc_type='natural')
a = spline.c
b = spline.derivative(nu=1)
d = np.zeros_like(a)
d[:-1] = 3 * (a[1:] - a[:-1]) / (x[1:] - x[:-1]) - 2 * b[:-1] - b[1:]
d[-1] = 0

# Construct the matrix A and vector b
n = len(x) - 1
A = np.zeros((4*n, 4*n))
b = np.zeros(4*n)

for i in range(n):
    # Add the equations for the ith cubic polynomial to the matrix A and vector b
    A[4*i, 4*i:4*i+4] = np.array([x[i]**3, x[i]**2, x[i], 1])
    A[4*i+1, 4*i:4*i+4] = np.array([x[i+1]**3, x[i+1]**2, x[i+1], 1])
    A[4*i+2, 4*i:4*i+4] = np.array([3*x[i+1]**2, 2*x[i+1], 1, 0])
    A[4*i+3, 4*i:4*i+4] = np.array([3*x[i+1]**2, 2*x[i+1], 1, 0])
    b[4*i] = a[i]
    b[4*i+1] = a[i+1]
    b[4*i+2] = b[i+1]
    b[4*i+3] = b[i+1]

# Solve the system of linear equations to find the coefficients of the cubic polynomials
x = np.linalg.solve(A, b)

# Print the matrix A, vector b, and vector x
print("Matrix A:")
print(A)
print("Vector b:")
print(b)
print("Vector x:")
print(x)

