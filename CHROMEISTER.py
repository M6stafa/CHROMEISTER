# In the name of God
# Bio Algorithm Project (Phase 3)
# Mohammad Mehdi Heydari (98209094)
# Mostafa Najafi (98209218)


# imports
import os
import re
import pathlib
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from tqdm.auto import tqdm


# Constants
VALID_NUCLEOTIDES = ['A', 'C', 'G', 'T']
COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}


# Classes
class Hit:
    def __init__(self, kmer_hash, kmer_pos):
        self.hit_count = 0
        self.hash = kmer_hash
        self.pos_in_x = kmer_pos
        self.pos_in_y = None
        self.repetition = False


class CHROMEISTER:
    def __init__(self, query_path, db_path, kmer_len=32, kmer_key_len=12, hash='z-based',
                 z=4, hit_matrix_dim=1000, query_name=None, db_name=None):

        self.query_path = query_path
        self.db_path = db_path
        self.kmer_len = kmer_len
        self.kmer_key_len = kmer_key_len
        self.hash = hash
        self.z = z
        self.hit_matrix_dim = hit_matrix_dim

        # Check and parse hash
        if hash != 'z-based' and re.match(f'^[\*1]{{{kmer_len}}}$', hash) is None:
            raise Exception('Invalid hash. Hash must be "z-based" or [*1]{k-mer length}')
        self._hash_idx = self._find_hash_idx()

        self.query_name = query_name or pathlib.Path(query_path).stem
        self.db_name = db_name or pathlib.Path(db_path).stem

        self.query_headers, self.query_len = self._get_fasta_info(self.query_path)
        self.db_headers, self.db_len = self._get_fasta_info(self.db_path)

        self.hit_table = None
        self.hit_matrix = None
        self.dot_plot = None
        self.score = None
        self.LSGRs = None

        self._pow4 = {
            'A': [0 for _ in range(kmer_len + 1)],
            'C': [1 * (4 ** i) for i in range(kmer_len + 1)],
            'G': [2 * (4 ** i) for i in range(kmer_len + 1)],
            'T': [3 * (4 ** i) for i in range(kmer_len + 1)],
        }
        self._progress_update = 10000


    def run(self, filter_mode='one', diag_len=4, neighbour_dist=1, kernel_width=2, dist_th=1.5,
            omit_lsgrs=False, sampling_value=5, diag_separation=10, hsp_th=3,
            save_output='mdpHLhs', output_dir='./outputs', verbose=0):

        if verbose: print('Creating Hit table...')
        self.create_hit_table(verbose=verbose)
        if verbose: print('Creating Hit matrix...')
        self.create_hit_matrix(filter_mode=filter_mode, verbose=verbose)
        if verbose: print('Creating dot-plot...')
        self.create_dot_plot(diag_len=diag_len, neighbour_dist=neighbour_dist, kernel_width=kernel_width)
        if verbose: print('Computing score...')
        self.compute_score(dist_th=dist_th, verbose=verbose)
        if not omit_lsgrs:
            if verbose: print('Finding LSGRs...')
            self.find_lsgrs(sampling_value=sampling_value, diag_separation=diag_separation, hsp_th=hsp_th, verbose=verbose)
        if output_dir is not None:
            if verbose: print('Saving Outputs...')
            self.write_outputs(output_dir, save_output=save_output, verbose=verbose)


    def create_hit_table(self, verbose=0):
        self.hit_table = {}

        # Read Database(db)
        curr_kmer = ''
        curr_pos = 0
        hit_counter = 0
        if verbose:
            print('Loading Database...', flush=True)
            progress_bar = tqdm(total=self.db_len)

        for c in self._read_fasta(self.db_path):
            if c == '>':
                curr_kmer = ''
                continue

            # New character
            curr_pos += 1
            if verbose and ((curr_pos + 1) % self._progress_update == 0):
                progress_bar.update(self._progress_update)

            if c in VALID_NUCLEOTIDES:
                curr_kmer += c
            else:  # It can be anything (including N, Y, X ...)
                curr_kmer = ''
                continue

            # If kmer become bigger than kmer_len, truncate it
            if len(curr_kmer) > self.kmer_len:
                curr_kmer = curr_kmer[-self.kmer_len:]

            # If reach the kmer size, store it in hit table
            if len(curr_kmer) == self.kmer_len:
                kmer_key = curr_kmer[:self.kmer_key_len]

                if kmer_key in self.hit_table:
                    if self.hit_table[kmer_key].repetition == False:
                        hit_counter -= 1
                        self.hit_table[kmer_key].repetition = True
                else:
                    hit_counter += 1
                    self.hit_table[kmer_key] = Hit(self._kmer_hash(curr_kmer), curr_pos)

                # Non overlapping kmers
                curr_kmer = ''

        if verbose:
            progress_bar.update(self.db_len - progress_bar.n)
            progress_bar.close()
            print(f'\nFound {hit_counter} unique and non-overlapping kmer in database\n', flush=True)

        # Read query
        curr_kmer = ''
        curr_kmer_rc = ''  # Reverse Complement
        curr_pos = 0
        hit_counter = 0
        if verbose:
            print('Loading Query...', flush=True)
            progress_bar = tqdm(total=self.query_len)

        for c in self._read_fasta(self.query_path):
            if c == '>':
                curr_kmer = ''
                curr_kmer_rc = ''
                continue

            # New character
            curr_pos += 1
            if verbose and ((curr_pos + 1) % self._progress_update == 0):
                progress_bar.update(self._progress_update)

            if c in VALID_NUCLEOTIDES:
                curr_kmer += c
                curr_kmer_rc = COMPLEMENT[c] + curr_kmer_rc
            else:  # It can be anything (including N, Y, X ...)
                curr_kmer = ''
                curr_kmer_rc = ''
                continue

            # If kmer become bigger than kmer_len, truncate it
            if len(curr_kmer) > self.kmer_len:
                curr_kmer = curr_kmer[-self.kmer_len:]
                curr_kmer_rc = curr_kmer_rc[:self.kmer_len]

            # If reach the kmer size, store it in hit table
            if len(curr_kmer) == self.kmer_len:
                for kmer in [curr_kmer, curr_kmer_rc]:
                    kmer_key = kmer[:self.kmer_key_len]

                    if kmer_key in self.hit_table:
                        hit = self.hit_table[kmer_key]

                        if (hit.repetition == False) and (hit.hash == self._kmer_hash(kmer)):
                            hit_counter += 1
                            hit.hit_count += 1
                            hit.pos_in_y = curr_pos

        if verbose:
            progress_bar.update(self.query_len - progress_bar.n)
            progress_bar.close()
            print(f'\nFound {hit_counter} hits with database\n', flush=True)


    def create_hit_matrix(self, filter_mode='one', verbose=0):
        self.filter_mode = filter_mode
        self.hit_matrix = np.zeros((self.hit_matrix_dim, self.hit_matrix_dim), dtype=np.int)

        pixel_size_query = self.hit_matrix_dim / self.query_len
        ratio_query = self.query_len / self.hit_matrix_dim

        pixel_size_db = self.hit_matrix_dim / self.db_len
        ratio_db = self.db_len / self.hit_matrix_dim

        i_r_fix = max(1.0, self.kmer_len * pixel_size_query)
        j_r_fix = max(1.0, self.kmer_len * pixel_size_db)

        if filter_mode == 'one':
            self.hit_count_value = 1
        else:  # filter_mode == 'min'
            # find minimum
            min_hit_count = np.inf
            for hit in self.hit_table.values():
                if hit.pos_in_y is not None and hit.hit_count < min_hit_count:
                    min_hit_count = hit.hit_count

            # find minumum after subtraction
            sub_min_hit_count = np.inf
            for hit in self.hit_table.values():
                if hit.pos_in_y is not None and 0 < (hit.hit_count - min_hit_count) < sub_min_hit_count:
                    sub_min_hit_count = hit.hit_count - min_hit_count
            if sub_min_hit_count == np.inf:
                sub_min_hit_count = 0

            self.hit_count_value = min_hit_count + sub_min_hit_count

        for hit in self.hit_table.values():
            if hit.hit_count == self.hit_count_value:
                # We plot it
                # Convert scale to hit_matrix
                redir_db = min(self.hit_matrix_dim - 1, int(hit.pos_in_x / ratio_db))
                redir_query = min(self.hit_matrix_dim - 1, int(hit.pos_in_y / ratio_query))
                i_r, j_r = i_r_fix, j_r_fix

                while (int(i_r) >= 1) and (int(j_r) >= 1):
                    if (redir_query - int(i_r) > 0) and (redir_db - int(j_r) > 0):
                        self.hit_matrix[redir_query - int(i_r)][redir_db - int(j_r)] += 1
                    else:
                        self.hit_matrix[redir_query][redir_db] += 1
                        break

                    i_r -= min(1.0, pixel_size_query)
                    j_r -= min(1.0, pixel_size_db)

        if verbose:
            print(f'Found {np.sum(self.hit_matrix)} unique hits for hit counts = {self.hit_count_value} and hash = {self.hash}')


    def create_dot_plot(self, diag_len=4, neighbour_dist=1, kernel_width=2):
        score_density = np.zeros_like(self.hit_matrix)

        ################ Start New Idea ################
        # Author Idea:
        # Only keeps max of rows and columns, make lines
        # thin and hard to detect
        # Our Idea:
        # keeps the max of rows and columns and their
        # close neighbours if have value > 0
        ################# End New Idea #################
        # Only keep max of rows
        for i in range(self.hit_matrix_dim):
            cmax_pos = np.argmax(self.hit_matrix[i, :])
            if self.hit_matrix[i, cmax_pos] > 0:
                # score_density[i, :] = 0  # not needed, we start from np.zeros
                for k in range(-neighbour_dist, neighbour_dist + 1):
                    idx = max(0, min(self.hit_matrix_dim - 1, i + k))
                    if self.hit_matrix[idx, cmax_pos] > 0:
                        score_density[idx, cmax_pos] = 1

        # Only keep max of columns
        for i in range(self.hit_matrix_dim):
            rmax_pos = np.argmax(self.hit_matrix[:, i])
            if self.hit_matrix[rmax_pos, i] > 0:
                score_density[:, i] = 0
                for k in range(-neighbour_dist, neighbour_dist + 1):
                    idx = max(0, min(self.hit_matrix_dim - 1, i + k))
                    if self.hit_matrix[rmax_pos, idx] > 0:
                        score_density[rmax_pos, idx] = 1

        # Diagonal expandation
        self.dot_plot = np.copy(score_density)
        diag_len_2 = diag_len // 2

        # Main diagonal expandation
        for i in range(diag_len + 1, self.hit_matrix_dim - (diag_len + 1)):
            for j in range(diag_len + 1, self.hit_matrix_dim - (diag_len + 1)):
                value = 0
                for k in range(-diag_len_2, diag_len_2 + 1):
                    if score_density[i + k, j + k] > 0:
                        value += 1

                if value > diag_len:
                    for k in range(1, diag_len + 2):
                        self.dot_plot[i + k, j + k] = 1
                        self.dot_plot[i - k, j - k] = 1

        # Anti diagonal expandation
        for i in range(diag_len + 1, self.hit_matrix_dim - (diag_len + 1)):
            for j in range(diag_len + 1, self.hit_matrix_dim - (diag_len + 1)):
                value = 0
                for k in range(-diag_len_2, diag_len_2 + 1):
                    if score_density[i - k, j + k] > 0:
                        value += 1

                if value > diag_len:
                    for k in range(1, diag_len + 2):
                        self.dot_plot[i - k, j + k] = 1
                        self.dot_plot[i + k, j - k] = 1

        # Kernel to remove single points
        ################ Start New Idea ################
        # Author Idea:
        # Only remove points with no neighbours
        # Our Idea:
        # Bigger kernel and remove points with low
        # number of neighbours
        ################# End New Idea #################
        kernel_th = kernel_width + 1
        for i in range(self.hit_matrix_dim):
            for j in range(self.hit_matrix_dim):
                min_i = max(0, i - kernel_width)
                max_i = min(self.hit_matrix_dim, i + kernel_width + 1)
                min_j = max(0, j - kernel_width)
                max_j = min(self.hit_matrix_dim, j + kernel_width + 1)

                value = np.sum(self.dot_plot[min_i:max_i, min_j:max_j])

                if value < kernel_th:
                    self.dot_plot[i, j] = 0


    def compute_score(self, dist_th=1.5, verbose=0):
        def dvec(i):
            a1 = np.argmax(self.dot_plot[:, i - 1])
            a2 = np.argmax(self.dot_plot[:, i])

            return abs(a2 - a1)

        self.score = 0

        dvec1 = dvec(1)
        dvec2 = dvec(2)
        dvec3 = dvec(3)

        for i in range(4, self.hit_matrix_dim):
            distance = np.mean([dvec1, dvec2, dvec3])

            if distance > dist_th or distance == 0:
                self.score += self.hit_matrix_dim

            dvec1 = dvec2
            dvec2 = dvec3
            dvec3 = dvec(i)

        self.score /= (self.hit_matrix_dim ** 2)

        if verbose:
            print(f'Score = {self.score:.03f}')


    def find_lsgrs(self, sampling_value=5, diag_separation=10, hsp_th=3, verbose=0):
        self.submat = self._downsample(self.dot_plot, sampling_value)
        ################ Start New Idea ################
        # Author Idea:
        # Calculate the main diagnol and anti diagonal HSPs together
        # Can't find close HSPs sometimes
        # Our Idea:
        main_diag_hsps = self._growing_regions(self.submat, 'right', th=hsp_th)
        antidiag_hsps = self._growing_regions(self.submat, 'left', th=hsp_th)
        hsps = np.vstack((main_diag_hsps, antidiag_hsps))
        ################# End New Idea #################
        events, event_types = self._detect_events(hsps, sampling_value, diag_separation)

        self.LSGRs = pd.DataFrame(
            data = [list(events[i]) + [event_types[i]] for i in range(len(events))],
            columns = ['x1', 'y1', 'x2', 'y2', 'len', 'event'],
        )

        if verbose:
            print('LSGRs:')
            print(self.LSGRs)


    def write_outputs(self, output_dir, save_output='mdpHLhs', verbose=0):
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # Save hit matrix (m)
        if 'm' in save_output:
            with open(os.path.join(output_dir, 'hit_matrix.mat'), 'w') as f:
                f.write(f'{self.db_len}\n')
                f.write(f'{self.query_len}\n')
                for i in range(self.hit_matrix_dim):
                    f.write(f'{" ".join(map(str, self.hit_matrix[i]))}\n')

        # Save dot plot (d)
        if 'd' in save_output:
            with open(os.path.join(output_dir, 'hits-XY.hits'), 'w') as f:
                f.write('X Y\n')
                for x in range(self.hit_matrix_dim):
                    for y in range(self.hit_matrix_dim):
                        if self.dot_plot[y][x] > 0:
                            f.write(f'{x} {y}\n')

        # Save plots (p and H)
        save_png = 'p' in save_output
        save_html = 'H' in save_output

        if save_png or save_html:
            self.draw_dot_plot(draw_lsgrs=False, show=False, output_path=os.path.join(output_dir, 'dot-plot'), save_png=save_png, save_html=save_html)
            if self.LSGRs is not None:
                self.draw_dot_plot(draw_lsgrs=True, show=False, output_path=os.path.join(output_dir, 'dot-plot with LSGRs'), save_png=save_png, save_html=save_html)

        # Save LSGRs (L)
        if 'L' in save_output and self.LSGRs is not None:
            with open(os.path.join(output_dir, 'events.txt'), 'w') as f:
                f.write(f'{self.db_len},{self.query_len}\n')
                f.write(self.LSGRs.to_csv(index=False, line_terminator='\n'))
                f.write('0,0,0,0,0,Null event\n')

        # Save headers (h)
        if 'h' in save_output:
            with open(os.path.join(output_dir, 'fastas_headers.csv'), 'w') as f:
                f.write(f'{os.path.split(self.query_path)[1]} (Y-Axis)\n')
                f.write(self.query_headers.to_csv(index=False, line_terminator='\n'))
                f.write('#\n')
                f.write(f'{os.path.split(self.db_path)[1]} (X-Axis)\n')
                f.write(self.db_headers.to_csv(index=False, line_terminator='\n'))
                f.write('#\n')

        # Save score (s)
        if 's' in save_output:
            with open(os.path.join(output_dir, 'score.txt'), 'w') as f:
                f.write(f'query: {self.query_name}, {self.query_len}\n')
                f.write(f'database: {self.db_name}, {self.db_len}\n')
                f.write(f'{self.score}\n')

        if verbose:
            print(f'Outputs saved in {os.path.abspath(output_dir)}')


    def draw_dot_plot(self, draw_lsgrs=False, show=False, output_path=None, save_png=True, save_html=True):
        fig = px.imshow(self.dot_plot, color_continuous_scale='gray')

        if draw_lsgrs:
            for index, lsgr in self.LSGRs.iterrows():
                fig.add_shape(
                    type="line",
                    xref="x",
                    yref="y",
                    x0=lsgr['x1'],
                    y0=lsgr['y1'],
                    x1=lsgr['x2'],
                    y1=lsgr['y2'],
                    line=dict(
                        color="red",
                        width=3,
                    ),
                )

        fig.update_xaxes(title=dict(text=self.db_name, standoff=0), side='top', tickvals=[10], ticktext=['►'])
        fig.update_yaxes(title=self.query_name, tickvals=[10], ticktext=['▼'])
        fig.update_layout(
            title = f'hash = {self.hash}, filter mode = {self.filter_mode}, hit counts = {self.hit_count_value}, Score = {self.score:.03f}',
            width = max(1000, self.hit_matrix_dim) + 200,
            height = max(1000, self.hit_matrix_dim) + 200 + 100,
            coloraxis_showscale = False,
            margin = dict(t=100),
        )

        if show:
            fig.show()
        else:
            if save_png:
                fig.write_image(f'{output_path}.png')
            if save_html:
                fig.write_html(f'{output_path}.html')


    # Downsample dot plot
    def _downsample(self, mat, downscale):
        l = len(mat)
        size = int(np.ceil(l / downscale))
        m = np.zeros((size, size))

        for i in range(l):
            for j in range(l):
                m_up = max(0, i - 1)
                m_down = min(l, i + 2)
                m_left = max(0, j - 1)
                m_right = min(l, j + 2)

                if np.sum(mat[m_up:m_down, m_left:m_right] > 0):
                    m[i // downscale, j // downscale] = 1

        return m


    # Find HSPs
    def _growing_regions(self, mat, direction, reward=6, penalty=-15, side_penalty=-3, max_hsps=500, th=3, wsize=7):
        if wsize % 2 == 0:
            raise Exception('wsize must be odd!')

        if direction not in ['left', 'right']:
            raise Exception('direction must be \'left\' or \'right\'!')

        mat = np.copy(mat)
        l = len(mat)
        HSPs = np.zeros((max_hsps, 5), dtype=np.int)
        num_hsps = 0  # number of founded HSPs
        lH = round(wsize / 2) - 1
        rH = round(wsize / 2) + 1

        i = 0
        while i < (l - th):
            value = np.max(mat[i, :]) * reward
            if value == 0:
                i += 1
                continue

            pos = np.argmax(mat[i, :])
            # these two hold ending frag
            end_frag = pos
            j = i

            j_valid = i  # has value in window

            count_penalties = 1

            while (value > 0) and (j < (l - 1)):
                # Reset position used
                for x in range(max(0, end_frag - 2), min(l - 1, end_frag + 2) + 1):
                    mat[max(0, j - 1), x] = 0
                    mat[max(0, j), x] = 0

                # Go for next
                j += 1
                # Check next
                if direction == 'left':
                    m_left = max(0, end_frag - lH)
                    m_right = min(l, end_frag + 1)
                else:
                    m_left = max(0, end_frag)
                    m_right = min(l, end_frag + lH + 1)

                window = mat[j, m_left:m_right]

                v = np.max(window)
                selected = np.argmax(window)

                # Make it rather go diagonally
                chose_diagonal = False
                if len(window) >= 2:
                    if direction == 'left' and v == window[-2]:
                        selected = len(window) - 2
                        chose_diagonal = True
                    if direction == 'right' and v == window[1]:
                        selected = 1
                        chose_diagonal = True

                if v == 0:  # If no similarity is found
                    value += count_penalties * penalty
                    count_penalties += 1
                else:  # Similarity is found
                    end_frag = m_left + selected  # To make the indexing

                    ################ Start New Idea ################
                    # Author Idea:
                    # Use the last row (most of the times it hasn't
                    # value != 0 and point to irrelevant position.
                    # Our Idea:
                    # Save the last row that has value != 0
                    j_valid = j
                    ################# End New Idea #################

                    if not chose_diagonal:
                        value += count_penalties * side_penalty
                        count_penalties += 1
                    else:
                        value += reward
                        count_penalties = 1

            if j_valid - i > th:
                HSPs[num_hsps, 0] = pos
                HSPs[num_hsps, 1] = i
                HSPs[num_hsps, 2] = end_frag
                HSPs[num_hsps, 3] = j_valid
                HSPs[num_hsps, 4] = abs(i - j_valid)

                num_hsps += 1

            if num_hsps >= max_hsps:
                break

        return HSPs[:num_hsps]


    # detect_events
    def _detect_events(self, HSPs, sampling, diag_separation):
        events = np.zeros_like(HSPs)
        event_types = ['' for _ in range(len(HSPs))]

        for i in range(len(HSPs)):
            is_inverted = False
            is_diagonal = True

            if HSPs[i, 0] > HSPs[i, 2]:
                is_inverted = True

            ################ Start New Idea ################
            # Author Idea:
            # if (abs(HSPs[i, 0] - HSPs[i, 1]) > diag_separation) and (abs(HSPs[i, 2] - HSPs[i, 3]) > diag_separation):
            #     is_diagonal = False
            # Our Idea:
            x_mid = HSPs[i, 0] + ((HSPs[i, 2] - HSPs[i, 0]) / 2)
            y_mid = HSPs[i, 1] + ((HSPs[i, 3] - HSPs[i, 1]) / 2)
            if abs(x_mid - y_mid) > diag_separation:
                is_diagonal = False
            ################# End New Idea #################

            events[i] = HSPs[i] * sampling

            if is_diagonal and (not is_inverted): event_types[i] = 'synteny block'
            if is_diagonal and is_inverted: event_types[i] = 'inversion'
            if (not is_diagonal) and (not is_inverted): event_types[i] = 'transposition'
            if (not is_diagonal) and is_inverted: event_types[i] = 'inverted transposition'

        return events, event_types


    # Read file generator
    def _read_fasta(self, file_path):
        with open(file_path, 'r') as f:
            while True:
                c = f.read(1)

                if not c: break  # End of file

                # Skip info lines
                if c == '>':
                    f.readline()
                    yield c  # signal new block
                elif c != '\n':
                    yield c.upper()  # new nucleotide


    def _get_fasta_info(self, file_path):
        num_nucleotides = 0
        headers = [[0, '', 0, 0]]
        index = 0

        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    index += 1
                    headers[index - 1][2] = num_nucleotides - headers[index - 1][2]
                    headers[index - 1][3] = num_nucleotides
                    headers.append([index, line[1:], num_nucleotides, num_nucleotides])
                else:
                    num_nucleotides += len(line)

        headers[index][2] = num_nucleotides - headers[index][2]
        headers[index][3] = num_nucleotides

        headers = pd.DataFrame(data=headers, columns=['ID', 'header_name', 'length (bp)', 'accumulated_length'])

        return headers.iloc[1:], num_nucleotides


    def _find_hash_idx(self):
        idx = []

        if self.hash == 'z-based':
            self.hash = ''
            for i in range(self.kmer_len):
                if i % self.z == 0:
                    self.hash += '1'
                else:
                    self.hash += '*'

        # hash = [*1]{kmer_len}
        for i, c in enumerate(self.hash):
            if c == '1':
                idx.append(i)

        return idx


    # Calculate kmer hash
    def _kmer_hash(self, kmer):
        value = 0

        for i in self._hash_idx:
            value += self._pow4[kmer[i]][self.kmer_len - (i + 1)]

        return value


# Run from console
if __name__ == '__main__':
    # imports
    import argparse
    from os import path


    # helper functions
    def file_path_type(file_path):
        if path.exists(file_path) and path.isfile(file_path):
            return file_path
        raise argparse.ArgumentTypeError('The path does not exist or isn\'t file!')


    # Parse args
    parser = argparse.ArgumentParser(description='CHROMEISTER')
    parser.add_argument('--db', required=True, metavar='DatabasePath', type=file_path_type, help='/path/to/db.fasta')
    parser.add_argument('--query', required=True, metavar='QueryPath', type=file_path_type, help='/path/to/query.fasta')
    parser.add_argument('--db-name', type=str, default=None, help='Database name used in outputs instead of filename')
    parser.add_argument('--query-name', type=str, default=None, help='Query name used in outputs instead of filename')
    parser.add_argument('--kmer-len', '-kmer', type=int, default=32, help='k-mer length')
    parser.add_argument('--kmer-key-len', '-kmer-key', type=int, default=12, help='k-mer key length')
    parser.add_argument('--hash', type=str, default='z-based',
                                  help='How use nucleotides to calculate k-mer hash, can be "z-based" or string of [*1]{k-mer length}')
    parser.add_argument('--z', '-z', type=int, default=4, help='The "z" in hash function')
    parser.add_argument('--dimension', '-dim', type=int, default=1000, help='size of the hit matrix')
    parser.add_argument('--out-dir', '-out', type=str, default='./outputs', help='/path/to/output/directory, can be "None" for not saving anything')
    parser.add_argument('--filter-mode', type=str, choices=['one', 'min'], default='one', help='How to filter hit counts for creating hit matrix')
    parser.add_argument('--diag-len', type=int, default=4, help='diagonal length for expandation')
    parser.add_argument('--neighbour-dist', type=int, default=1, help='how far neighbours keep from max in rows and columns')
    parser.add_argument('--kernel-width', type=int, default=2, help='kernel width for removing points = 2 * value + 1')
    parser.add_argument('--dist-th', type=float, default=1.5, help='distance threshold used in compute score')
    parser.add_argument('--omit-lsgrs', action='store_true', default=False, help='Don\'t calculate and save LSGRs')
    parser.add_argument('--sampling-value', '-sampling', type=int, default=3, help='downsample factor for find HSPs better')
    parser.add_argument('--diag-separation', type=int, default=10, help='dot-plot downsample factor to find HSPs better')
    parser.add_argument('--hsp-th', type=int, default=3, help='min size of HSP')
    parser.add_argument('--save-output', type=str, default='mdpHLhs', help='Which output need to save: m: hit_matrix.mat, d: hits-XY.hits, p: plots.png, H: plots.html, L: events.txt, h: fastas_headers.csv, s: score.txt')
    parser.add_argument('--verbose', '-v', type=int, choices=[0, 1], default=1, help='Print some info during run!')

    args = parser.parse_args()

    if args.out_dir == 'None':
        args.out_dir = None

    if args.out_dir is not None:
        # Write args in output
        if not os.path.exists(args.out_dir):
            os.mkdir(args.out_dir)

        with open(os.path.join(args.out_dir, 'args.txt'), 'w') as f:
            for arg in vars(args):
                f.write(f'{arg} = {getattr(args, arg)}\n')

    if args.verbose:
        # Show args
        print('Arguments:')
        for arg in vars(args):
            print(f'{arg} = {getattr(args, arg)}')

    # Run CHROMEISTER
    if args.verbose: print('\nRunning CHROMEISTER...')

    chromeister = CHROMEISTER(
        query_path = args.query,
        db_path = args.db,
        kmer_len = args.kmer_len,
        kmer_key_len = args.kmer_key_len,
        hash = args.hash,
        z = args.z,
        hit_matrix_dim = args.dimension,
        query_name = args.query_name,
        db_name = args.db_name,
    )

    chromeister.run(
        filter_mode = args.filter_mode,
        diag_len = args.diag_len,
        neighbour_dist = args.neighbour_dist,
        kernel_width = args.kernel_width,
        dist_th = args.dist_th,
        omit_lsgrs = args.omit_lsgrs,
        sampling_value = args.sampling_value,
        diag_separation = args.diag_separation,
        hsp_th = args.hsp_th,
        output_dir = args.out_dir,
        save_output = args.save_output,
        verbose = args.verbose,
    )
