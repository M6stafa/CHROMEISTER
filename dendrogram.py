import plotly.figure_factory as ff
import numpy as np
import scipy.spatial.distance as ssd


def neighbor_joining(D):
    # -*- coding: utf-8 -*-

    # Implement the Neighbor Joining Algorithm
    # http://rosalind.info/problems/ba7e


    # imports
    import math
    from copy import deepcopy


    # Helper Functions
    def neighbor_joining(D, m, removed = []):
        n = len(D) - len(removed)
        D_len = len(D)

        if n == 2:
            # T ← tree consisting of a single edge of length D1,2
            T = {}
            i, j = [i for i in range(D_len) if i not in removed]
            T[i] = {j: D[i][j]}
            T[j] = {i: D[i][j]}
            return T

        # D* ← neighbor-joining matrix constructed from the distance matrix D
        # find elements i and j such that D*i,j is a minimum non-diagonal element of D*
        total_distance = [0 for i in range(D_len)]
        for i in range(D_len):
            new_sum = 0
            for j in range(D_len):
                if j not in removed:
                    new_sum += D[i][j]
            total_distance[i] = new_sum
        Dstar = [[0 for i in range(D_len)] for j in range(D_len)]
        min_value = math.inf
        indexes = None
        for i in range(D_len):
            for j in range(i + 1, D_len):
                if i in removed or j in removed:
                    continue

                Dstar[i][j] = Dstar[j][i] = (n - 2) * D[i][j] - total_distance[i] - total_distance[j]
                if Dstar[i][j] < min_value:
                    min_value = Dstar[i][j]
                    indexes = (i, j)

        # Δ ← (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2)
        i, j = indexes
        delta = (total_distance[i] - total_distance[j]) / (n - 2)

        # limbLengthi ← (1/2)(Di,j + Δ)
        # limbLengthj ← (1/2)(Di,j - Δ)
        limb_len_i = (D[i][j] + delta) / 2
        limb_len_j = (D[i][j] - delta) / 2

        # add a new row/column m to D so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k
        D.append([0 for k in range(len(D))])
        for k in range(len(D)):
            D[k].append(0)
        for k in range(len(D) - 1):
            if k in removed:
                continue
            D[k][m] = D[m][k] = (D[k][i] + D[k][j] - D[i][j]) / 2

        # D ← D with rows i and j removed
        # D ← D with columns i and j removed
        removed.append(i)
        removed.append(j)

        # recursive call
        T = neighbor_joining(D, m + 1, removed)

        # add two new limbs (connecting node m with leaves i and j) to the tree T
        # assign length limbLengthi to Limb(i)
        # assign length limbLengthj to Limb(j)
        T[m][i] = limb_len_i
        T[m][j] = limb_len_j
        T[i] = { m: limb_len_i }
        T[j] = { m: limb_len_j }

        return T


    # Main
    # n = int(input())
    # D = []
    # for i in range(n):
    #     D.append(list(map(int, input().split(' '))))

    tree = neighbor_joining(D, len(D))

    # for v in sorted(tree.keys()):
    #     for w in sorted(tree[v].keys()):
    #         print(f'{v}->{w}:{tree[v][w]:.3f}')

    return tree


names = [
    '[1] A/turkey/Ontario/FAV110-4/2009(H1N1)',
    '[2] A/mallard/Nova_Scotia/00088/2010(H1N1)',
    '[3] A/thick-billed_murre/Canada/1871/2011(H1N1)',
    '[4] A/duck/Guangxi/030D/2009(H1N1)',
    '[5] A/mallard/France/691/2002(H1N1)',
    '[6] A/duck/Hokkaido/w73/2007(H1N1)',
    '[7] A/pintail/Miyagi/1472/2008(H1N1)',
    '[8] A/mallard/Korea/KNU_YP09/2009(H1N1)',
    '[9] A/mallard/Maryland/352/2002(H1N1)',
    '[10] A/mallard/MD/26/2003(H1N1)',
    '[11] A/dunlin/Alaska/44421-660/2008(H1N1)',
    '[12] A/mallard/Minnesota/Sg-00620/2008(H1N1)',
    '[13] A/turkey/Virginia/4135/2014(H1N1)',
    '[14] A/chicken/Eastern_China/XH222/2008(H5N1)',
    '[15] A/duck/Eastern_China/JS017/2009(H5N1)',
    '[16] A/chicken/Yunnan/chuxiong01/2005(H5N1)',
    '[17] A/chicken/Germany/R3234/2007(H5N1)',
    '[18] A/domestic_duck/Germany/R1772/07(H5N1)',
    '[19] A/wild_bird/Hong_Kong/07035-1/2011(H5N1)',
    '[20] A/Chicken/Hong_Kong/822.1/01_(H5N1)',
    '[21] A/chicken/Miyazaki/10/2011(H5N1)',
    '[22] A/chicken/Korea/es/2003(H5N1)',
    '[23] A/mandarin_duck/Korea/K10-483/2010(H5N1)',
    '[24] A/turkey/VA/505477-18/2007(H5N1)',
    '[25] A/American_black_duck/NB/2538/2007(H7N3)',
    '[26] A/American_black_duck/New_Brunswick/02490/2007(H7N3)',
    '[27] A/American_green-winged_teal/California/44242-906/2007(H7N3)',
    '[28] A/avian/Delaware_Bay/226/2006(H7N3)',
    '[29] A/chicken/British_Columbia/GSC_human_B/04(H7N3)',
    '[30] A/chicken/Rizhao/713/2013(H7N9)',
    '[31] A/chicken/Jiangsu/1021/2013(H7N9)',
    '[32] A/duck/Jiangxi/3096/2009(H7N9)',
    '[33] A/wild_duck/Korea/SH19-47/2010(H7N9)',
    '[34] A/turkey/Minnesota/1/1988(H7N9)',
    '[35] A/mallard/Minnesota/AI09-3770/2009(H7N9)',
    '[36] A/mallard/Postdam/178-4/83(H2N2)',
    '[37] A/duck/Hong_Kong/319/1978(H2N2)',
    '[38] A/emperor_goose/Alaska/44297-260/2007(H2N2)',
]

colors = [
    'blue',
    'blue',
    'blue',
    'blue',
    'blue',
    'blue',
    'blue',
    'blue',
    'blue',
    'blue',
    'blue',
    'blue',
    'blue',
    'red',
    'red',
    'red',
    'red',
    'red',
    'red',
    'red',
    'red',
    'red',
    'red',
    'red',
    'green',
    'green',
    'green',
    'green',
    'green',
    'black',
    'black',
    'black',
    'black',
    'black',
    'black',
    'magenta',
    'magenta',
    'magenta',
]

X = {0: {1: 0.781, 2: 0.719, 3: 0.801, 4: 0.71, 5: 0.753, 6: 0.698, 7: 0.73, 8: 0.741, 9: 0.705, 10: 0.805, 11: 0.747, 12: 0.108, 13: 0.804, 14: 0.816, 15: 0.852, 16: 0.806, 17: 0.864, 18: 0.828, 19: 0.747, 20: 0.828, 21: 0.804, 22: 0.828, 23: 0.722, 24: 0.985, 25: 0.968, 26: 0.982, 27: 0.968, 28: 0.985, 29: 0.982, 30: 0.982, 31: 0.968, 32: 0.996, 33: 0.986, 34: 0.996, 35: 0.976, 36: 0.986, 37: 0.986}, 1: {0: 0.781, 2: 0.132, 3: 0.702, 4: 0.748, 5: 0.786, 6: 0.727, 7: 0.733, 8: 0.421, 9: 0.408, 10: 0.118, 11: 0.108, 12: 0.862, 13: 0.78, 14: 0.781, 15: 0.781, 16: 0.828, 17: 0.792, 18: 0.84, 19: 0.766, 20: 0.805, 21: 0.757, 22: 0.816, 23: 0.487, 24: 0.974, 25: 0.996, 26: 0.996, 27: 0.996, 28: 0.985, 29: 0.946, 30: 0.946, 31: 0.926, 32: 0.968, 33: 0.986, 34: 0.982, 35: 0.962, 36: 0.952, 37: 0.986}, 2: {0: 0.719, 1: 0.132, 3: 0.685, 4: 0.741, 5: 0.775, 6: 0.751, 7: 0.727, 8: 0.407, 9: 0.39, 10: 0.104, 11: 0.08, 12: 0.816, 13: 0.793, 14: 0.793, 15: 0.806, 16: 0.839, 17: 0.806, 18: 0.84, 19: 0.782, 20: 0.818, 21: 0.772, 22: 0.829, 23: 0.541, 24: 0.963, 25: 0.996, 26: 0.996, 27: 0.996, 28: 0.974, 29: 0.954, 30: 0.954, 31: 0.94, 32: 0.982, 33: 0.986, 34: 0.996, 35: 0.964, 36: 0.954, 37: 0.966}, 3: {0: 0.801, 1: 0.702, 2: 0.685, 4: 0.648, 5: 0.348, 6: 0.406, 7: 0.663, 8: 0.805, 9: 0.835, 10: 0.851, 11: 0.813, 12: 0.849, 13: 0.432, 14: 0.449, 15: 0.411, 16: 0.51, 17: 0.452, 18: 0.586, 19: 0.417, 20: 0.502, 21: 0.37, 22: 0.503, 23: 0.756, 24: 0.988, 25: 0.988, 26: 0.988, 27: 0.98, 28: 0.988, 29: 0.982, 30: 0.988, 31: 0.972, 32: 0.996, 33: 0.975, 34: 0.996, 35: 0.974, 36: 0.981, 37: 0.968}, 4: {0: 0.71, 1: 0.748, 2: 0.741, 3: 0.648, 5: 0.394, 6: 0.446, 7: 0.499, 8: 0.719, 9: 0.716, 10: 0.731, 11: 0.742, 12: 0.697, 13: 0.554, 14: 0.566, 15: 0.554, 16: 0.544, 17: 0.566, 18: 0.648, 19: 0.513, 20: 0.602, 21: 0.506, 22: 0.603, 23: 0.743, 24: 0.996, 25: 0.967, 26: 0.967, 27: 0.982, 28: 0.996, 29: 0.982, 30: 0.982, 31: 0.982, 32: 0.996, 33: 0.966, 34: 0.996, 35: 0.976, 36: 0.986, 37: 0.966}, 5: {0: 0.753, 1: 0.786, 2: 0.775, 3: 0.348, 4: 0.394, 6: 0.267, 7: 0.478, 8: 0.708, 9: 0.733, 10: 0.755, 11: 0.709, 12: 0.768, 13: 0.604, 14: 0.616, 15: 0.638, 16: 0.621, 17: 0.627, 18: 0.698, 19: 0.527, 20: 0.686, 21: 0.554, 22: 0.686, 23: 0.794, 24: 0.985, 25: 0.996, 26: 0.996, 27: 0.982, 28: 0.985, 29: 0.996, 30: 0.996, 31: 0.968, 32: 0.996, 33: 0.976, 34: 0.996, 35: 0.952, 36: 0.962, 37: 0.976}, 6: {0: 0.698, 1: 0.727, 2: 0.751, 3: 0.406, 4: 0.446, 5: 0.267, 7: 0.466, 8: 0.672, 9: 0.722, 10: 0.741, 11: 0.726, 12: 0.711, 13: 0.591, 14: 0.592, 15: 0.58, 16: 0.625, 17: 0.591, 18: 0.662, 19: 0.532, 20: 0.639, 21: 0.508, 22: 0.64, 23: 0.723, 24: 0.985, 25: 0.996, 26: 0.996, 27: 0.982, 28: 0.985, 29: 0.982, 30: 0.982, 31: 0.954, 32: 0.996, 33: 0.966, 34: 0.996, 35: 0.966, 36: 0.976, 37: 0.976}, 7: {0: 0.73, 1: 0.733, 2: 0.727, 3: 0.663, 4: 0.499, 5: 0.478, 6: 0.466, 8: 0.733, 9: 0.758, 10: 0.738, 11: 0.767, 12: 0.731, 13: 0.602, 14: 0.602, 15: 0.553, 16: 0.618, 17: 0.614, 18: 0.696, 19: 0.553, 20: 0.663, 21: 0.529, 22: 0.674, 23: 0.778, 24: 0.974, 25: 0.967, 26: 0.967, 27: 0.954, 28: 0.985, 29: 0.982, 30: 0.982, 31: 0.954, 32: 0.996, 33: 0.965, 34: 0.996, 35: 0.976, 36: 0.986, 37: 0.966}, 8: {0: 0.741, 1: 0.421, 2: 0.407, 3: 0.805, 4: 0.719, 5: 0.708, 6: 0.672, 7: 0.733, 9: 0.221, 10: 0.409, 11: 0.416, 12: 0.819, 13: 0.794, 14: 0.794, 15: 0.806, 16: 0.803, 17: 0.782, 18: 0.828, 19: 0.779, 20: 0.839, 21: 0.793, 22: 0.839, 23: 0.316, 24: 0.974, 25: 0.982, 26: 0.996, 27: 0.996, 28: 0.974, 29: 0.968, 30: 0.968, 31: 0.94, 32: 0.996, 33: 0.986, 34: 0.982, 35: 0.974, 36: 0.974, 37: 0.986}, 9: {0: 0.705, 1: 0.408, 2: 0.39, 3: 0.835, 4: 0.716, 5: 0.733, 6: 0.722, 7: 0.758, 8: 0.221, 10: 0.396, 11: 0.392, 12: 0.786, 13: 0.759, 14: 0.759, 15: 0.772, 16: 0.769, 17: 0.771, 18: 0.84, 19: 0.745, 20: 0.838, 21: 0.76, 22: 0.849, 23: 0.389, 24: 0.974, 25: 0.982, 26: 0.996, 27: 0.996, 28: 0.974, 29: 0.968, 30: 0.968, 31: 0.926, 32: 0.996, 33: 0.986, 34: 0.996, 35: 0.974, 36: 0.974, 37: 0.986}, 10: {0: 0.805, 1: 0.118, 2: 0.104, 3: 0.851, 4: 0.731, 5: 0.755, 6: 0.741, 7: 0.738, 8: 0.409, 9: 0.396, 11: 0.069, 12: 0.859, 13: 0.792, 14: 0.793, 15: 0.805, 16: 0.814, 17: 0.804, 18: 0.852, 19: 0.781, 20: 0.817, 21: 0.769, 22: 0.828, 23: 0.503, 24: 0.963, 25: 0.996, 26: 0.996, 27: 0.996, 28: 0.974, 29: 0.954, 30: 0.954, 31: 0.94, 32: 0.968, 33: 0.986, 34: 0.982, 35: 0.976, 36: 0.966, 37: 0.966}, 11: {0: 0.747, 1: 0.108, 2: 0.08, 3: 0.813, 4: 0.742, 5: 0.709, 6: 0.726, 7: 0.767, 8: 0.416, 9: 0.392, 10: 0.069, 12: 0.839, 13: 0.793, 14: 0.793, 15: 0.806, 16: 0.839, 17: 0.784, 18: 0.84, 19: 0.794, 20: 0.796, 21: 0.772, 22: 0.807, 23: 0.551, 24: 0.952, 25: 0.996, 26: 0.996, 27: 0.996, 28: 0.963, 29: 0.954, 30: 0.954, 31: 0.94, 32: 0.968, 33: 0.986, 34: 0.996, 35: 0.964, 36: 0.954, 37: 0.966}, 12: {0: 0.108, 1: 0.862, 2: 0.816, 3: 0.849, 4: 0.697, 5: 0.768, 6: 0.711, 7: 0.731, 8: 0.819, 9: 0.786, 10: 0.859, 11: 0.839, 13: 0.804, 14: 0.816, 15: 0.84, 16: 0.806, 17: 0.853, 18: 0.828, 19: 0.748, 20: 0.804, 21: 0.78, 22: 0.804, 23: 0.707, 24: 0.985, 25: 0.968, 26: 0.982, 27: 0.968, 28: 0.985, 29: 0.982, 30: 0.982, 31: 0.968, 32: 0.996, 33: 0.986, 34: 0.996, 35: 0.976, 36: 0.986, 37: 0.986}, 13: {0: 0.804, 1: 0.78, 2: 0.793, 3: 0.432, 4: 0.554, 5: 0.604, 6: 0.591, 7: 0.602, 8: 0.794, 9: 0.759, 10: 0.792, 11: 0.793, 12: 0.804, 14: 0.099, 15: 0.267, 16: 0.36, 17: 0.342, 18: 0.387, 19: 0.38, 20: 0.365, 21: 0.286, 22: 0.365, 23: 0.781, 24: 0.986, 25: 0.988, 26: 0.988, 27: 0.98, 28: 0.986, 29: 0.974, 30: 0.98, 31: 0.98, 32: 0.996, 33: 0.965, 34: 0.988, 35: 0.986, 36: 0.986, 37: 0.974}, 14: {0: 0.816, 1: 0.781, 2: 0.793, 3: 0.449, 4: 0.566, 5: 0.616, 6: 0.592, 7: 0.602, 8: 0.794, 9: 0.759, 10: 0.793, 11: 0.793, 12: 0.816, 13: 0.099, 15: 0.28, 16: 0.347, 17: 0.328, 18: 0.386, 19: 0.38, 20: 0.38, 21: 0.275, 22: 0.38, 23: 0.792, 24: 0.986, 25: 0.988, 26: 0.988, 27: 0.98, 28: 0.986, 29: 0.974, 30: 0.98, 31: 0.98, 32: 0.996, 33: 0.965, 34: 0.988, 35: 0.986, 36: 0.986, 37: 0.964}, 15: {0: 0.852, 1: 0.781, 2: 0.806, 3: 0.411, 4: 0.554, 5: 0.638, 6: 0.58, 7: 0.553, 8: 0.806, 9: 0.772, 10: 0.805, 11: 0.806, 12: 0.84, 13: 0.267, 14: 0.28, 16: 0.395, 17: 0.341, 18: 0.455, 19: 0.396, 20: 0.452, 21: 0.274, 22: 0.455, 23: 0.782, 24: 0.976, 25: 0.988, 26: 0.988, 27: 0.988, 28: 0.976, 29: 0.985, 30: 0.988, 31: 0.988, 32: 0.996, 33: 0.965, 34: 0.988, 35: 0.976, 36: 0.986, 37: 0.976}, 16: {0: 0.806, 1: 0.828, 2: 0.839, 3: 0.51, 4: 0.544, 5: 0.621, 6: 0.625, 7: 0.618, 8: 0.803, 9: 0.769, 10: 0.814, 11: 0.839, 12: 0.806, 13: 0.36, 14: 0.347, 15: 0.395, 17: 0.01, 18: 0.475, 19: 0.349, 20: 0.443, 21: 0.186, 22: 0.429, 23: 0.756, 24: 0.972, 25: 0.98, 26: 0.98, 27: 0.98, 28: 0.98, 29: 0.974, 30: 0.98, 31: 0.98, 32: 0.996, 33: 0.975, 34: 0.988, 35: 0.966, 36: 0.976, 37: 0.966}, 17: {0: 0.864, 1: 0.792, 2: 0.806, 3: 0.452, 4: 0.566, 5: 0.627, 6: 0.591, 7: 0.614, 8: 0.782, 9: 0.771, 10: 0.804, 11: 0.784, 12: 0.853, 13: 0.342, 14: 0.328, 15: 0.341, 16: 0.01, 18: 0.475, 19: 0.383, 20: 0.442, 21: 0.185, 22: 0.428, 23: 0.782, 24: 0.966, 25: 0.98, 26: 0.98, 27: 0.98, 28: 0.976, 29: 0.974, 30: 0.98, 31: 0.98, 32: 0.996, 33: 0.975, 34: 0.988, 35: 0.976, 36: 0.986, 37: 0.964}, 18: {0: 0.828, 1: 0.84, 2: 0.84, 3: 0.586, 4: 0.648, 5: 0.698, 6: 0.662, 7: 0.696, 8: 0.828, 9: 0.84, 10: 0.852, 11: 0.84, 12: 0.828, 13: 0.387, 14: 0.386, 15: 0.455, 16: 0.475, 17: 0.475, 19: 0.489, 20: 0.44, 21: 0.427, 22: 0.44, 23: 0.865, 24: 0.976, 25: 0.98, 26: 0.98, 27: 0.98, 28: 0.966, 29: 0.985, 30: 0.988, 31: 0.98, 32: 0.996, 33: 0.986, 34: 0.988, 35: 0.976, 36: 0.986, 37: 0.984}, 19: {0: 0.747, 1: 0.766, 2: 0.782, 3: 0.417, 4: 0.513, 5: 0.527, 6: 0.532, 7: 0.553, 8: 0.779, 9: 0.745, 10: 0.781, 11: 0.794, 12: 0.748, 13: 0.38, 14: 0.38, 15: 0.396, 16: 0.349, 17: 0.383, 18: 0.489, 20: 0.467, 21: 0.259, 22: 0.465, 23: 0.747, 24: 0.988, 25: 0.98, 26: 0.98, 27: 0.98, 28: 0.988, 29: 0.985, 30: 0.988, 31: 0.988, 32: 0.996, 33: 0.965, 34: 0.988, 35: 0.976, 36: 0.986, 37: 0.963}, 20: {0: 0.828, 1: 0.805, 2: 0.818, 3: 0.502, 4: 0.602, 5: 0.686, 6: 0.639, 7: 0.663, 8: 0.839, 9: 0.838, 10: 0.817, 11: 0.796, 12: 0.804, 13: 0.365, 14: 0.38, 15: 0.452, 16: 0.443, 17: 0.442, 18: 0.44, 19: 0.467, 21: 0.375, 22: 0.037, 23: 0.793, 24: 0.976, 25: 0.988, 26: 0.988, 27: 0.98, 28: 0.966, 29: 0.963, 30: 0.972, 31: 0.972, 32: 0.988, 33: 0.976, 34: 0.996, 35: 0.986, 36: 0.996, 37: 0.974}, 21: {0: 0.804, 1: 0.757, 2: 0.772, 3: 0.37, 4: 0.506, 5: 0.554, 6: 0.508, 7: 0.529, 8: 0.793, 9: 0.76, 10: 0.769, 11: 0.772, 12: 0.78, 13: 0.286, 14: 0.275, 15: 0.274, 16: 0.186, 17: 0.185, 18: 0.427, 19: 0.259, 20: 0.375, 22: 0.367, 23: 0.746, 24: 0.986, 25: 0.98, 26: 0.98, 27: 0.988, 28: 0.986, 29: 0.974, 30: 0.98, 31: 0.972, 32: 0.996, 33: 0.965, 34: 0.98, 35: 0.976, 36: 0.986, 37: 0.964}, 22: {0: 0.828, 1: 0.816, 2: 0.829, 3: 0.503, 4: 0.603, 5: 0.686, 6: 0.64, 7: 0.674, 8: 0.839, 9: 0.849, 10: 0.828, 11: 0.807, 12: 0.804, 13: 0.365, 14: 0.38, 15: 0.455, 16: 0.429, 17: 0.428, 18: 0.44, 19: 0.465, 20: 0.037, 21: 0.367, 23: 0.793, 24: 0.976, 25: 0.988, 26: 0.988, 27: 0.98, 28: 0.966, 29: 0.963, 30: 0.972, 31: 0.972, 32: 0.988, 33: 0.976, 34: 0.996, 35: 0.986, 36: 0.996, 37: 0.974}, 23: {0: 0.722, 1: 0.487, 2: 0.541, 3: 0.756, 4: 0.743, 5: 0.794, 6: 0.723, 7: 0.778, 8: 0.316, 9: 0.389, 10: 0.503, 11: 0.551, 12: 0.707, 13: 0.781, 14: 0.792, 15: 0.782, 16: 0.756, 17: 0.782, 18: 0.865, 19: 0.747, 20: 0.793, 21: 0.746, 22: 0.793, 24: 0.98, 25: 0.988, 26: 0.996, 27: 0.996, 28: 0.98, 29: 0.974, 30: 0.98, 31: 0.964, 32: 0.996, 33: 0.986, 34: 0.996, 35: 0.962, 36: 0.972, 37: 0.976}, 24: {0: 0.985, 1: 0.974, 2: 0.963, 3: 0.988, 4: 0.996, 5: 0.985, 6: 0.985, 7: 0.974, 8: 0.974, 9: 0.974, 10: 0.963, 11: 0.952, 12: 0.985, 13: 0.986, 14: 0.986, 15: 0.976, 16: 0.972, 17: 0.966, 18: 0.976, 19: 0.988, 20: 0.976, 21: 0.986, 22: 0.976, 23: 0.98, 25: 0.033, 26: 0.156, 27: 0.324, 28: 0.272, 29: 0.974, 30: 0.968, 31: 0.968, 32: 0.968, 33: 0.972, 34: 0.954, 35: 0.948, 36: 0.972, 37: 0.96}, 25: {0: 0.968, 1: 0.996, 2: 0.996, 3: 0.988, 4: 0.967, 5: 0.996, 6: 0.996, 7: 0.967, 8: 0.982, 9: 0.982, 10: 0.996, 11: 0.996, 12: 0.968, 13: 0.988, 14: 0.988, 15: 0.988, 16: 0.98, 17: 0.98, 18: 0.98, 19: 0.98, 20: 0.988, 21: 0.98, 22: 0.988, 23: 0.988, 24: 0.033, 26: 0.146, 27: 0.328, 28: 0.325, 29: 0.968, 30: 0.968, 31: 0.968, 32: 0.968, 33: 0.975, 34: 0.954, 35: 0.956, 36: 0.976, 37: 0.966}, 26: {0: 0.982, 1: 0.996, 2: 0.996, 3: 0.988, 4: 0.967, 5: 0.996, 6: 0.996, 7: 0.967, 8: 0.996, 9: 0.996, 10: 0.996, 11: 0.996, 12: 0.982, 13: 0.988, 14: 0.988, 15: 0.988, 16: 0.98, 17: 0.98, 18: 0.98, 19: 0.98, 20: 0.988, 21: 0.98, 22: 0.988, 23: 0.996, 24: 0.156, 25: 0.146, 27: 0.257, 28: 0.289, 29: 0.954, 30: 0.954, 31: 0.982, 32: 0.954, 33: 0.975, 34: 0.954, 35: 0.956, 36: 0.976, 37: 0.966}, 27: {0: 0.968, 1: 0.996, 2: 0.996, 3: 0.98, 4: 0.982, 5: 0.982, 6: 0.982, 7: 0.954, 8: 0.996, 9: 0.996, 10: 0.996, 11: 0.996, 12: 0.968, 13: 0.98, 14: 0.98, 15: 0.988, 16: 0.98, 17: 0.98, 18: 0.98, 19: 0.98, 20: 0.98, 21: 0.988, 22: 0.98, 23: 0.996, 24: 0.324, 25: 0.328, 26: 0.257, 28: 0.211, 29: 0.968, 30: 0.968, 31: 0.968, 32: 0.968, 33: 0.973, 34: 0.982, 35: 0.956, 36: 0.976, 37: 0.996}, 28: {0: 0.985, 1: 0.985, 2: 0.974, 3: 0.988, 4: 0.996, 5: 0.985, 6: 0.985, 7: 0.985, 8: 0.974, 9: 0.974, 10: 0.974, 11: 0.963, 12: 0.985, 13: 0.986, 14: 0.986, 15: 0.976, 16: 0.98, 17: 0.976, 18: 0.966, 19: 0.988, 20: 0.966, 21: 0.986, 22: 0.966, 23: 0.98, 24: 0.272, 25: 0.325, 26: 0.289, 27: 0.211, 29: 0.962, 30: 0.956, 31: 0.968, 32: 0.968, 33: 0.972, 34: 0.982, 35: 0.948, 36: 0.972, 37: 0.984}, 29: {0: 0.982, 1: 0.946, 2: 0.954, 3: 0.982, 4: 0.982, 5: 0.996, 6: 0.982, 7: 0.982, 8: 0.968, 9: 0.968, 10: 0.954, 11: 0.954, 12: 0.982, 13: 0.974, 14: 0.974, 15: 0.985, 16: 0.974, 17: 0.974, 18: 0.985, 19: 0.985, 20: 0.963, 21: 0.974, 22: 0.963, 23: 0.974, 24: 0.974, 25: 0.968, 26: 0.954, 27: 0.968, 28: 0.962, 30: 0.213, 31: 0.481, 32: 0.252, 33: 0.8, 34: 0.803, 35: 0.975, 36: 0.982, 37: 0.975}, 30: {0: 0.982, 1: 0.946, 2: 0.954, 3: 0.988, 4: 0.982, 5: 0.996, 6: 0.982, 7: 0.982, 8: 0.968, 9: 0.968, 10: 0.954, 11: 0.954, 12: 0.982, 13: 0.98, 14: 0.98, 15: 0.988, 16: 0.98, 17: 0.98, 18: 0.988, 19: 0.988, 20: 0.972, 21: 0.98, 22: 0.972, 23: 0.98, 24: 0.968, 25: 0.968, 26: 0.954, 27: 0.968, 28: 0.956, 29: 0.213, 31: 0.329, 32: 0.098, 33: 0.702, 34: 0.671, 35: 0.966, 36: 0.976, 37: 0.966}, 31: {0: 0.968, 1: 0.926, 2: 0.94, 3: 0.972, 4: 0.982, 5: 0.968, 6: 0.954, 7: 0.954, 8: 0.94, 9: 0.926, 10: 0.94, 11: 0.94, 12: 0.968, 13: 0.98, 14: 0.98, 15: 0.988, 16: 0.98, 17: 0.98, 18: 0.98, 19: 0.988, 20: 0.972, 21: 0.972, 22: 0.972, 23: 0.964, 24: 0.968, 25: 0.968, 26: 0.982, 27: 0.968, 28: 0.968, 29: 0.481, 30: 0.329, 32: 0.216, 33: 0.701, 34: 0.607, 35: 0.976, 36: 0.976, 37: 0.976}, 32: {0: 0.996, 1: 0.968, 2: 0.982, 3: 0.996, 4: 0.996, 5: 0.996, 6: 0.996, 7: 0.996, 8: 0.996, 9: 0.996, 10: 0.968, 11: 0.968, 12: 0.996, 13: 0.996, 14: 0.996, 15: 0.996, 16: 0.996, 17: 0.996, 18: 0.996, 19: 0.996, 20: 0.988, 21: 0.996, 22: 0.988, 23: 0.996, 24: 0.968, 25: 0.968, 26: 0.954, 27: 0.968, 28: 0.968, 29: 0.252, 30: 0.098, 31: 0.216, 33: 0.672, 34: 0.64, 35: 0.966, 36: 0.976, 37: 0.966}, 33: {0: 0.986, 1: 0.986, 2: 0.986, 3: 0.975, 4: 0.966, 5: 0.976, 6: 0.966, 7: 0.965, 8: 0.986, 9: 0.986, 10: 0.986, 11: 0.986, 12: 0.986, 13: 0.965, 14: 0.965, 15: 0.965, 16: 0.975, 17: 0.975, 18: 0.986, 19: 0.965, 20: 0.976, 21: 0.965, 22: 0.976, 23: 0.986, 24: 0.972, 25: 0.975, 26: 0.975, 27: 0.973, 28: 0.972, 29: 0.8, 30: 0.702, 31: 0.701, 32: 0.672, 34: 0.37, 35: 0.983, 36: 0.97, 37: 0.983}, 34: {0: 0.996, 1: 0.982, 2: 0.996, 3: 0.996, 4: 0.996, 5: 0.996, 6: 0.996, 7: 0.996, 8: 0.982, 9: 0.996, 10: 0.982, 11: 0.996, 12: 0.996, 13: 0.988, 14: 0.988, 15: 0.988, 16: 0.988, 17: 0.988, 18: 0.988, 19: 0.988, 20: 0.996, 21: 0.98, 22: 0.996, 23: 0.996, 24: 0.954, 25: 0.954, 26: 0.954, 27: 0.982, 28: 0.982, 29: 0.803, 30: 0.671, 31: 0.607, 32: 0.64, 33: 0.37, 35: 0.986, 36: 0.976, 37: 0.976}, 35: {0: 0.976, 1: 0.962, 2: 0.964, 3: 0.974, 4: 0.976, 5: 0.952, 6: 0.966, 7: 0.976, 8: 0.974, 9: 0.974, 10: 0.976, 11: 0.964, 12: 0.976, 13: 0.986, 14: 0.986, 15: 0.976, 16: 0.966, 17: 0.976, 18: 0.976, 19: 0.976, 20: 0.986, 21: 0.976, 22: 0.986, 23: 0.962, 24: 0.948, 25: 0.956, 26: 0.956, 27: 0.956, 28: 0.948, 29: 0.975, 30: 0.966, 31: 0.976, 32: 0.966, 33: 0.983, 34: 0.986, 36: 0.546, 37: 0.71}, 36: {0: 0.986, 1: 0.952, 2: 0.954, 3: 0.981, 4: 0.986, 5: 0.962, 6: 0.976, 7: 0.986, 8: 0.974, 9: 0.974, 10: 0.966, 11: 0.954, 12: 0.986, 13: 0.986, 14: 0.986, 15: 0.986, 16: 0.976, 17: 0.986, 18: 0.986, 19: 0.986, 20: 0.996, 21: 0.986, 22: 0.996, 23: 0.972, 24: 0.972, 25: 0.976, 26: 0.976, 27: 0.976, 28: 0.972, 29: 0.982, 30: 0.976, 31: 0.976, 32: 0.976, 33: 0.97, 34: 0.976, 35: 0.546, 37: 0.749}, 37: {0: 0.986, 1: 0.986, 2: 0.966, 3: 0.968, 4: 0.966, 5: 0.976, 6: 0.976, 7: 0.966, 8: 0.986, 9: 0.986, 10: 0.966, 11: 0.966, 12: 0.986, 13: 0.974, 14: 0.964, 15: 0.976, 16: 0.966, 17: 0.964, 18: 0.984, 19: 0.963, 20: 0.974, 21: 0.964, 22: 0.974, 23: 0.976, 24: 0.96, 25: 0.966, 26: 0.966, 27: 0.996, 28: 0.984, 29: 0.975, 30: 0.966, 31: 0.976, 32: 0.966, 33: 0.983, 34: 0.976, 35: 0.71, 36: 0.749}}


dist = [[0 for _ in range(len(names))] for _ in range(len(names))]
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        dist[i][j] = dist[j][i] = X[i][j]


T = neighbor_joining(dist)
# print(T)
# print(type(T))
# print(max(T.keys()))
# for i in range(max(T.keys()) + 1):
#     if i not in T:
#         print(i)

def DFS(T, start, end):
    seen = {start: True}
    stack = [(start, 0)]
    while len(stack) > 0:
        node, d = stack.pop()
        for child in T[node]:
            if child in seen:
                continue

            if child == end:
                return d + T[node][child]

            seen[child] = True
            stack.append((child, d + T[node][child]))

    return None


dist_np = np.zeros((len(names), len(names)))
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        dist_np[i][j] = dist_np[j][i] = DFS(T, i, j)


def colorize(text, color):
    return fr'$\color{{{color}}}{{{text}}}$'


ticktext = [colorize(names[i], colors[i]) for i in range(len(names))]
# print(ticktext)

fig = ff.create_dendrogram(dist_np, distfun= lambda x: ssd.squareform(x), orientation='left', labels=ticktext)
fig.update_layout(
    width = 2000,
    height = 800,
    margin = dict(l=300),
)
fig.show()
