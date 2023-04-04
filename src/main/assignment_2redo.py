import numpy as np
np.set_printoptions(precision=7, suppress=True, linewidth=100)

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

print(result)
print("\n --------------2--------------------- \n")


def divided_difference_table(x_points, y_points):
    size = len(x_points)
    matrix = np.zeros((size, size))
    list = []

    for i in range(size):
        matrix[i][0] = y_points[i]

    for i in range(1, size):
        for j in range(1, i + 1):
            matrix[i][j] = (matrix[i][j - 1] - matrix[i - 1][j - 1]) / (x_points[i] - x_points[i - j])

            if i == j:
                list.append(matrix[i][j])

    print(list)
    return matrix


def get_approximate_result(matrix, x_points, value, start):
    reoccuring_x_span = 1
    reoccuring_px_result = start

    for index in range(1, len(matrix)):
        polynomial_coefficient = matrix[index][index]

        reoccuring_x_span *= (value - x_points[index - 1])

        mlt_operation = polynomial_coefficient * reoccuring_x_span

        reoccuring_px_result += mlt_operation

    print(reoccuring_px_result)

x_points = [7.2, 7.4, 7.5, 7.6]
y_points = [23.5492, 25.3913, 26.8224, 27.4589]
divided_table = divided_difference_table(x_points, y_points)


print("\n ---------------3-------------------- \n")

get_approximate_result(divided_table, x_points, 7.3, y_points[0])
print()

print("\n ---------------4-------------------- \n")


def apply_div_diff(matrix):
    size = len(matrix)
    for i in range(2, size):
        for j in range(2, i + 2):

            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue

            # something get left and diag left
            left = matrix[i][j - 1]
            diag_left = matrix[i - 1][j - 1]
            numerator = left - diag_left

            denominator = matrix[i][0] - matrix[i - j + 1][0]

            operation = numerator / denominator
            matrix[i][j] = operation
    return matrix


def hermite_interpolation(x_points, y_points, slopes):
    # main difference with hermite's method , using instances with x

    num_of_points = len(x_points)
    matrix = np.zeros((num_of_points * 2, num_of_points * 2))

    # populate x values

    index = 0
    for x in range(0, len(matrix), 2):
        matrix[x][0] = x_points[index]
        matrix[x + 1][0] = x_points[index]
        index += 1

    # prepopulate y values
    index = 0
    for x in range(0, len(matrix), 2):
        matrix[x][1] = y_points[index]
        matrix[x + 1][1] = y_points[index]
        index += 1

    # prepopulate with derivatives (every other row)
    index = 0
    for x in range(1, len(matrix), 2):
        matrix[x][2] = slopes[index]
        index += 1

    apply_div_diff(matrix)
    print(matrix)

x_points = [3.6, 3.8, 3.9]
y_points = [1.675, 1.436, 1.318]
slopes = [-1.195, -1.188, -1.182]
hermite_interpolation(x_points, y_points, slopes)

print("\n ---------------5-------------------- \n")

def cubic_spline_interpolation(x_points, y_points):
    size = len(x_points)
    matrix = np.zeros((size, size))
    matrix[0][0] = matrix[size - 1][size - 1] = 1

    for i in range(1, size - 1):
        index = i - 1
        for j in range(index, size, 2):
            matrix[i][j] = x_points[index + 1] - x_points[index]
            index += 1

    for i in range(1, size - 1):
        matrix[i][i] = 2 * ((x_points[i + 1] - x_points[i]) + (x_points[i] - x_points[i - 1]))

    print(np.matrix(matrix), "\n")

    spline_condition = np.zeros((size))

    for i in range (1, size - 1):
        first_term = (3 / (x_points[i + 1] - x_points[i])) * (y_points[i + 1] - y_points[i])
        second_term = (3 / (x_points[i] - x_points[i - 1])) * (y_points[i] - y_points[i - 1])
        spline_condition[i] = first_term - second_term

    print(np.array(spline_condition), "\n")
    print(np.array(np.linalg.solve(matrix, spline_condition)))

x_points, y_points = [2, 5, 8, 10], [3, 5, 7, 9]
cubic_spline_interpolation(x_points, y_points)
