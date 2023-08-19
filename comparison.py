import numpy as np

def levenshtein_distance(original, solution):
    """Returns Levenstheina distance between two sequences"""
    distances = np.zeros((len(original) + 1, len(solution) + 1))

    for o in range(len(original) + 1):
        distances[o][0] = o
    for s in range(len(solution) + 1):
        distances[0][s] = s

    a = 0
    b = 0
    c = 0

    for o in range(1, len(original) + 1):
        for s in range(1, len(solution) + 1):
            if original[o - 1] == solution[s - 1]:
                distances[o][s] = distances[o - 1][s - 1]
            else:
                a = distances[o][s - 1]
                b = distances[o - 1][s]
                c = distances[o - 1][s - 1]

                if a <= b and a <= c:
                    distances[o][s] = a + 1
                elif b <= a and b <= c:
                    distances[o][s] = b + 1
                else:
                    distances[o][s] = c + 1

    return distances[len(original)][len(solution)]

def similarity(distance, original, solution):
    """Returns % of similarity between original and result sequence"""
    return round((1.0 - distance / max(len(original), len(solution))) * 100, 2)

def pairwise_alignment(seq1, seq2):
    """Aligns two sequences - Needlman's algorithm"""
    match = 1
    mismatch = -1
    gap = -2

    n = len(seq1)
    m = len(seq2)
    scoring_matrix = [[0 for j in range(m)] for i in range(n)]

    for i in range(n):
        for j in range(m):
            if seq1[i] == seq2[j]:
                scoring_matrix[i][j] = match
            else:
                scoring_matrix[i][j] = mismatch

    matrix = [[0 for j in range(m+1)] for i in range(n+1)]

    for i in range(n+1):
        matrix[i][0] = i * gap
    for j in range(m+1):
        matrix[0][j] = j * gap

    for i in range(1, n+1):
        for j in range(1, m+1):
            matrix[i][j] = max(matrix[i-1][j-1] + scoring_matrix[i-1][j-1],
                               matrix[i-1][j] + gap,
                               matrix[i][j-1] + gap)

    aligned1 = ""
    aligned2 = ""
    while(n>0 and m>0):

        if(n>0 and m>0 and matrix[n][m] == matrix[n-1][m-1] + scoring_matrix[n-1][m-1]):
            aligned1 = seq1[n-1] + aligned1
            aligned2 = seq2[m-1] + aligned2
            n -= 1
            m -= 1

        elif(n>0 and matrix[n][m] == matrix[n-1][m] + gap):
            aligned1 = seq1[n-1] + aligned1
            aligned2 = " " + aligned2
            n -= 1

        else:
            aligned1 = " " + aligned1
            aligned2 = seq2[m-1] + aligned2
            m -= 1

    return aligned1, aligned2
