#!/usr/bin/python3
from typing import List

from which_pyqt import PYQT_VER

import math

INF = math.inf

K = 7

D = 3

NO_ALIGNMENT = "No Alignment Possible"

LEFT = "left"

DIAGONAL = "diagonal"

DOWN = "down"

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import random

MAXINDELS = 5

MATCH = -3
INDEL = 5
SUB = 1

NUM_CHARACTERS = 100


def get_num_chars(alignment_1, alignment_2):
    alignment_1 = alignment_1[::-1]
    alignment_1 = alignment_1[::-1]
    alignment_2 = alignment_2[::-1]

    alignment1 = alignment_1[:NUM_CHARACTERS]
    alignment2 = alignment_2[:NUM_CHARACTERS]
    return alignment1, alignment2


def left_banded(alignment_1, alignment_2, sequence_2, seq_to_pos, x):
    x = x - 1
    seq_to_pos = seq_to_pos - 1
    alignment_1 = alignment_1 + '-'
    alignment_2 = alignment_2 + sequence_2[seq_to_pos]
    return alignment_1, alignment_2, seq_to_pos, x


def down_banded(alignment_1, alignment_2, sequence_1, x, y):
    y = y - 1
    x = x + 1
    alignment_1 = alignment_1 + sequence_1[x]
    alignment_2 = alignment_2 + '-'
    return alignment_1, alignment_2, x, y


def diagonal_banded(alignment_1, alignment_2, sequence_1, sequence_2, seq_to_pos, y):
    y = y - 1
    seq_to_pos = seq_to_pos - 1
    alignment_1 = alignment_1 + sequence_1[y]
    alignment_2 = alignment_2 + sequence_2[seq_to_pos]
    return alignment_1, alignment_2, seq_to_pos, y


def down_unrestricted(alignment_1, alignment_2, sequence_1, y):
    y = y - 1
    alignment_1 = alignment_1 + sequence_1[y]
    alignment_2 = alignment_2 + "-"
    return alignment_1, alignment_2, y


def diagonal_unrestricted(alignment_1, alignment_2, sequence_1, sequence_2, x, y):
    x = x - 1
    y = y - 1
    alignment_1 = alignment_1 + sequence_1[y]
    alignment_2 = alignment_2 + sequence_2[x]
    return alignment_1, alignment_2, x, y


def left_unrestricted(alignment_1, alignment_2, sequence_2, x):
    x = x - 1
    alignment_1 = alignment_1 + "-"
    alignment_2 = alignment_2 + sequence_2[x]
    return alignment_1, alignment_2, x


class GeneSequencing:
    matrix: List[List[int]]
    banded_matrix: List[List[int]]

    def __init__(self):
        self.banded = None
        self.MaxCharactersToAlign = None
        pass

    # Unrestricted:
        # Time: O(unrestricted_map) + O(unrestricted_align) = O(n) + O(n^2) = O(n + n^2) = O(n^2)
        # Space: O(unrestricted_map) + O(unrestricted_align) = O(1) + O(n^2) = O(1 + n^2) = O(n^2)
    # Banded:
        # Time: O(banded_map) + O(banded_align) = O(n) + O(kn) = O(n + kn) = O(kn)
        # Space: O(banded_map) + O(banded_align) = O(1) + O(kn) = O(1 + kn) = O(kn)
    def align(self, sequence_1, sequence_2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        sequence_1, sequence_2 = sequence_1[:align_length], sequence_2[:align_length]

        if not banded:  # Unrestricted
            align1, align2 = self.generate_unrestricted_map(sequence_1, sequence_2)  # O(n) time, O(1) space
            score = self.unrestricted_alignment(sequence_1, sequence_2)  # O(n^2) for time and space
        else:  # Banded
            pos = self.banded_alignment(sequence_1, sequence_2, K)  # O(kn) for time and space
            if pos == INF:
                return {'align_cost': pos, 'seqi_first100': NO_ALIGNMENT, 'seqj_first100': NO_ALIGNMENT}
            score = (self.banded_matrix[len(sequence_1)][pos])[0]
            align1, align2 = self.generate_banded_map(sequence_1, sequence_2, pos)  # O(n) time, O(1) space

        alignment1, alignment2 = get_num_chars(align1, align2)  # const

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}

    # Time: O(2n) = O(n)
    # Space: constant (see unrestricted func for more detail)
    def generate_banded_map(self, sequence_1, sequence_2, current_position):
        alignment_1 = alignment_2 = ""
        seq_to_pos = len(sequence_2) - 1
        x_val, y_val = current_position, len(sequence_1) - 1

        if abs(len(sequence_1) - len(sequence_2)) >= D:
            return NO_ALIGNMENT, NO_ALIGNMENT

        while True:  # n + n = 2n
            cur_val = self.banded_matrix[y_val][x_val]  # everything in here is const time and space
            direction = cur_val[1]

            if direction == DOWN:
                alignment_1, alignment_2, x_val, y_val = down_banded(alignment_1, alignment_2, sequence_1, x_val, y_val)
            elif direction == LEFT:
                alignment_1, alignment_2, seq_to_pos, x_val = left_banded(alignment_1, alignment_2, sequence_2, seq_to_pos,
                                                                          x_val)
            else:
                alignment_1, alignment_2, seq_to_pos, y_val = diagonal_banded(alignment_1, alignment_2, sequence_1, sequence_2,
                                                                              seq_to_pos, y_val)
            if y_val == 0 and x_val == D:
                break
        return alignment_1, alignment_2

    # Time: O(n(100 + c + k)) = O(n(100 k)) = O(kn(100)) = O(kn)
    # Space: O(kn)
    def banded_alignment(self, sequence_1, sequence_2, k):
        global direction
        self.banded_matrix = [[0 for j in range(0, k)] for i in range(0, len(sequence_1))]
        sequence_2, sequence_1 = sequence_2, ' ' + sequence_1

        for i in range(0, 4):
            self.banded_matrix[i][D - i] = (i * 5, DOWN)
            self.banded_matrix[0][D + i] = (i * 5, LEFT)

        seq_to_beginning = 0
        limit_max = False
        total_chars = 4
        starting_index = D
        position = 0
        r0 = range(1, len(sequence_1))
        r1 = range(0, total_chars)
        r2 = range(0, 6)

        for i in r0:  # n
            if seq_to_beginning + total_chars + 1 >= len(sequence_2):
                temp_seq = sequence_2[seq_to_beginning:len(sequence_2) + 1]

            else:
                temp_seq = sequence_2[seq_to_beginning:total_chars + seq_to_beginning]

            for j in r1:  # 100
                match = False  # everything in here is constant

                if sequence_1[i] == temp_seq[j]:
                    match = True
                dist_to_beat = INF

                if j - 1 >= starting_index:
                    left_val = self.banded_matrix[i][j - 1]
                    dist_to_beat = left_val[0] + 5
                    direction = LEFT

                index_j = starting_index + j
                if index_j + 1 < K:
                    top_val = self.banded_matrix[i - 1][index_j + 1]
                    top = top_val[0] + 5

                    if top >= dist_to_beat:
                        continue
                    dist_to_beat = top
                    direction = DOWN

                diag_val = self.banded_matrix[i - 1][index_j]

                if match:
                    top_left = diag_val[0] - D

                else:
                    top_left = diag_val[0] + 1

                if top_left >= dist_to_beat:
                    pass
                else:

                    dist_to_beat = top_left
                    direction = DIAGONAL

                self.banded_matrix[i][index_j] = [dist_to_beat, direction]

            if limit_max:
                seq_to_beginning = seq_to_beginning + 1

            if total_chars < k and not limit_max:
                total_chars = total_chars + 1

                if total_chars == k:
                    limit_max = True

            if i >= len(sequence_2) - D:
                total_chars = total_chars - 1

            if starting_index != 0:
                starting_index = starting_index - 1

        sequence_1_len = len(sequence_1) - 1
        lowest_val = (self.banded_matrix[sequence_1_len][0])[0]

        chars_total = total_chars + 1
        if chars_total > K:
            return INF

        for i in r1:  # 100
            temp = (self.banded_matrix[sequence_1_len][i])[0]

            if temp <= lowest_val:
                lowest_val = temp
                position = i

        for i in r2:  # k
            if self.banded_matrix[sequence_1_len][i + 1] == 0:
                break

            position = i

        return position + 1

    # Time: O(2n) = O(n)
    # Space: constant
    def generate_unrestricted_map(self, seq_1, seq_2):
        alignment_1 = alignment_2 = ""
        x_value, y_value = len(seq_2), len(seq_1)

        while True:  # Worst case scenario is traverse all the way up and over matrix: n + n
            direction = (self.matrix[y_value][x_value])[1]

            if direction == LEFT:
                alignment_1, alignment_2, x_value = left_unrestricted(alignment_1, alignment_2, seq_2, x_value) # const time

            elif direction == DIAGONAL:
                alignment_1, alignment_2, x_value, y_value = diagonal_unrestricted(alignment_1, alignment_2, seq_1,
                                                                                   seq_2, x_value, y_value)  # const time
            else:
                alignment_1, alignment_2, y_value = down_unrestricted(alignment_1, alignment_2, seq_1, y_value)  # const time

            if y_value == 0 and x_value == 0:
                break

        return alignment_1, alignment_2

    # Time: O(n + m +n^2 + c) = O(n + n +n^2) = O(2n +n^2) = O(n^2)
    # Space: O(n^2)
    def unrestricted_alignment(self, seq_1, seq_2):
        seq_2, seq_1 = ' ' + seq_2, ' ' + seq_1

        self.matrix = [[0 for j in range(0, len(seq_2))] for i in range(0, len(seq_1))]

        for i in range(0, len(seq_2)):  # n
            i_ = i * 5
            self.matrix[0][i] = (i_, LEFT)  # const
            # print(i_)
        for i in range(0, len(seq_1)): # m = n
            i_ = i * 5
            self.matrix[i][0] = (i_, DOWN)  # const
            # print(i_)
        for i in range(1, len(seq_1)):
            for j in range(1, len(seq_2)):  # n * m = n * n = n^2
                match = False
                if seq_2[j] == seq_1[i]:
                    match = True

                direction = LEFT
                left_val = self.matrix[i][j - 1]
                best_distance = left_val[0] + 5

                top_val = self.matrix[i - 1][j]
                top = top_val[0] + 5  # TOP
                if top <= best_distance:
                    best_distance = top
                    direction = DOWN

                top_left_val = self.matrix[i - 1][j - 1]
                if match:
                    top_left = top_left_val[0] - D
                else:
                    top_left = top_left_val[0] + 1
                if top_left > best_distance:
                    pass
                else:
                    best_distance = top_left
                    direction = DIAGONAL
                self.matrix[i][j] = (best_distance, direction)
        return self.matrix[len(seq_1) - 1][len(seq_2) - 1][0]
