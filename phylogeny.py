# In the name of God
# Bio Algorithm Project (Phase 3)
# Mohammad Mehdi Heydari (98209094)
# Mostafa Najafi (98209218)


# imports
from os import path
from multiprocessing import Pool
import json
import math
import pickle

import numpy as np
from tqdm import tqdm
import plotly.figure_factory as ff
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import linkage

from CHROMEISTER import CHROMEISTER


# Functions
def get_score(args):
    i, j, query_path, db_path = args
    chro = CHROMEISTER(query_path, db_path, kmer_len=16, kmer_key_len=8, z=4)
    chro.run(omit_lsgrs=True, output_dir=None, verbose=0)

    return i, j, chro.score


def create_dist_matrix(genomes_dir, meta, verbose=1):
    # run processes
    with open(meta, 'r') as f:
        meta = json.load(f)

    dist_matrix = [[math.inf for _ in range(len(meta))] for _ in range(len(meta))]
    for i in range(len(meta)):
        dist_matrix[i][i] = 0

    input_args = []

    for i in range(len(meta)):
        path_i = path.join(genomes_dir, meta[i]['file_name'])

        for j in range(i + 1, len(meta)):
            path_j = path.join(genomes_dir, meta[j]['file_name'])
            input_args.append((i, j, path_i, path_j))
            input_args.append((i, j, path_j, path_i))

    if verbose: progress_bar = tqdm(total=len(input_args))

    with Pool() as pool:
        for i, j, score in pool.imap_unordered(get_score, input_args):
            dist_matrix[i][j] = dist_matrix[j][i] = min(score, dist_matrix[i][j])

            if verbose: progress_bar.update(1)

    if verbose: progress_bar.close()

    return dist_matrix, meta


def create_meg_content(dist_matrix, meta):
    output = ''
    output += f'#mega\n!Format DataType=Distance DataFormat=LowerLeft NTaxa={len(meta)};\n\n'

    for i in range(len(meta)):
        output += f'[{i + 1}] #{meta[i]["genome_name"].replace(" ", "")}\n'
    output += '\n'

    output += f'[{" ".join([str(i + 1) for i in range(len(meta))])}]\n'
    for i in range(len(meta)):
        output += f'[{i + 1}] {" ".join([str(dist_matrix[i][j]) for j in range(0, i)])}\n'

    return output


def neighbor_joining(dist_matrix_np):
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

            linkage = [[i, j, D[i][j]/2, D[i][j]/2]]

            return T, linkage

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
        T, linkage = neighbor_joining(D, m + 1, removed)

        # add two new limbs (connecting node m with leaves i and j) to the tree T
        # assign length limbLengthi to Limb(i)
        # assign length limbLengthj to Limb(j)
        T[m][i] = limb_len_i
        T[m][j] = limb_len_j
        T[i] = { m: limb_len_i }
        T[j] = { m: limb_len_j }

        # update linkage
        linkage.append([i, j, limb_len_i, limb_len_j])

        return T, linkage


    # main
    dist_matrix = dist_matrix_np.tolist()
    n = len(dist_matrix)

    tree, linkage = neighbor_joining(dist_matrix, n)

    node_height = [0 for i in range(n)]
    linkage = linkage[::-1]
    num_orginals = 0

    for i in range(len(linkage)):
        node1, node2, limb_len1, limb_len2 = linkage[i]
        mid_node = n + i

        if 0 <= node1 < n: num_orginals += 1
        if 0 <= node2 < n: num_orginals += 1

        dist = max(node_height[node1] + limb_len1, node_height[node2] + limb_len2)
        node_height.append(dist) # mid_node height

        linkage[i] = [node1, node2, dist, num_orginals]

    return tree, linkage


def draw_plot(dist_matrix, meta, algortihm, save_path, title, show=False):
    def colorize(text, color):
        return fr'$\color{{{color}}}{{\verb|⬤ {text}|}}$'


    dist_matrix = np.array(dist_matrix)
    unique_colors = set([m['color'] for m in meta])

    if algortihm == 'UPGMA':
        algo_name = 'UPGMA'
        Z = linkage(ssd.squareform(dist_matrix), method='average', optimal_ordering=False)
    else:  # NJ (neighbor_joining)
        algo_name = 'Neighbor Joining'
        Z = neighbor_joining(dist_matrix)[1]


    fig = ff.create_dendrogram(
        dist_matrix,
        # distfun = distfun,
        linkagefun = lambda x: Z,
        orientation = 'right',
        labels = [colorize(m['genome_name'], m['color']) for m in meta],
        color_threshold = math.inf,  # Z[-len(unique_colors) + 1][2],
        colorscale = ['#212121'] * 8,
    )

    fig.update_layout(
        title = f'{title} ({algo_name})',
        title_x = 0.5,
        width = 2000,
        height = len(meta) * 35,
        paper_bgcolor = 'rgba(255, 255, 255, 1)',
        plot_bgcolor = 'rgba(245, 245, 245, 1)',
    )
    fig.update_yaxes(
        range = [0, len(meta) * 10],
        side='right',
        ticks = '',
        tickfont = dict(size=15),
    )
    fig.update_xaxes(
        range = [-1.01, 0],
        tickmode = 'array',
        tickvals = np.arange(-1, 0.1, 0.1),
        ticktext = list(map(lambda x: f'{x:.1f}', np.arange(1, -0.1, -0.1))),
    )

    if show:
        fig.show()

    if save_path is not None:
        fig.write_image(save_path)


if __name__ == '__main__':
    # imports
    import argparse


    # Functions
    def dir_path_type(dir_path):
        if path.exists(dir_path) and path.isdir(dir_path):
            return dir_path
        raise argparse.ArgumentTypeError('The path does not exist or isn\'t directory!')


    def file_path_type(file_path):
        if path.exists(file_path) and path.isfile(file_path):
            return file_path
        raise argparse.ArgumentTypeError('The path does not exist or isn\'t file!')


    # Parse args
    parser = argparse.ArgumentParser(description='Create phylogeny distance matrix based on CHROMEISTER scores')
    parser.add_argument('mode', type=str, choices=['compute', 'result'])
    parser.add_argument('--meta', required=True, type=file_path_type, help='/path/to/meta.json')
    parser.add_argument('--dist-mat', type=str, required=True, help='/path/to/dist_matrix.pickle')
    parser.add_argument('--genomes-dir', type=dir_path_type, help='/path/to/genomes/directory for compute mode')
    parser.add_argument('--plot', type=str, default=None, help='/path/to/plot.png')
    parser.add_argument('--plot-title', type=str, default='', help='plot title')
    parser.add_argument('--no-show', action='store_true', help='Don\'t show the plot')
    parser.add_argument('--meg', type=str, default=None, help='/path/to/output.meg')
    parser.add_argument('--algorithm', type=str, choices=['UPGMA', 'NJ'], default='UPGMA', help='phlogenetic tree creation algortihm (NJ = neighbor joining)')
    parser.add_argument('--verbose', '-v', type=int, choices=[0, 1], default=1, help='Print some info during run!')

    args = parser.parse_args()

    if args.verbose:
        # Show args
        print('Arguments:')
        for arg in vars(args):
            print(f'{arg} = {getattr(args, arg)}')

    # Compute
    if args.mode == 'compute':
        dist_matrix, meta = create_dist_matrix(args.genomes_dir, args.meta, args.verbose)

        # Save distance matrix
        with open(args.dist_mat, 'wb') as f:
            pickle.dump(dist_matrix, f)

    # Results
    if args.mode == 'result':
        with open(args.dist_mat, 'rb') as f:
            dist_matrix = pickle.load(f)

        with open(args.meta, 'r') as f:
            meta = json.load(f)

    draw_plot(dist_matrix, meta, args.algorithm, args.plot, args.plot_title, show=not args.no_show)

    if args.meg is not None:
        meg_content = create_meg_content(dist_matrix, meta)

        with open(args.meg, 'w') as f:
            f.write(meg_content)

