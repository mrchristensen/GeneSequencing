#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

# Other Constants
NUM_CHARACTERS = 100


def get_char_alignment(alignment_1, alignment_2):
	alignment_1 = alignment_1[::-1]
	alignment_1 = alignment_1[::-1]
	alignment_2 = alignment_2[::-1]

	alignment1 = alignment_1[:NUM_CHARACTERS]
	alignment2 = alignment_2[:NUM_CHARACTERS]
	return alignment1, alignment2


class GeneSequencing:

	def __init__(self):
		self.banded = None
		self.MaxCharactersToAlign = None
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		seq1, seq2 = seq1[:align_length], seq2[:align_length]
		if not banded:  # Unrestricted
			score = self.unrestricted_alignment(seq1, seq2)
			align1, align2 = self.generate_unrestricted_map(seq1, seq2)
		else:  # Banded
			pos = self.banded_alignment(seq1, seq2, 7)
			if pos == math.inf:
				score = math.inf
				align1 = align2 = "No Possible Alignment"
				return {'align_cost': score, 'seqi_first100': align1, 'seqj_first100': align2}
			score = (self.bandedMatrix[len(seq1)][pos])[0]

			align1, align2 = self.generate_banded_map(seq1, seq2, pos)

		alignment1, alignment2 = get_char_alignment(align1, align2)

		return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}

	def generate_banded_map(self, seq1, seq2, pos):
		alignment_1 = alignment_2 = ""
		seq_to_pos = len(seq2) - 1
		x_val, y_val = pos, len(seq1) - 1

		if abs(len(seq1) - len(seq2)) >= 3:
			return "No Possible Alignment", "No Possible Alignment"
		while True:
			curVal = self.bandedMatrix[y_val][x_val]
			direction = curVal[1]

			if direction == "down":
				y_val = y_val - 1
				x_val = x_val + 1
				alignment_1 += seq1[x_val]
				alignment_2 += '-'
			elif direction == "left":
				x_val = x_val - 1
				seq_to_pos -= 1
				alignment_1 += '-'
				alignment_2 += seq2[seq_to_pos]
			else:
				y_val = y_val - 1
				seq_to_pos -= 1
				alignment_1 += seq1[y_val]
				alignment_2 += seq2[seq_to_pos]
			if y_val == 0 and x_val == 3:
				break
		return alignment_1, alignment_2

	def banded_alignment(self, seq1, seq2, k):
		seq1 = ' ' + seq1
		seq2 = seq2
		self.bandedMatrix = [[0 for j in range(0, k)] for i in range(0, len(seq1))]
		for i in range(0,4):
			self.bandedMatrix[0][3 + i] = (i * 5, "left")
		for i in range(0,4):
			self.bandedMatrix[i][3 - i] = (i * 5, "down")

		startingIndex = 3
		totalChars = 4
		seq2Beginning = 0
		max = False
		#need to keep updating sequence 2 as going along matrix
		for i in range(1, len(seq1)):
			stringCheck = seq2Beginning + totalChars + 1 >= len(seq2)
			if stringCheck:
				seq2Copy = seq2[seq2Beginning:len(seq2) + 1]
			else:
				seq2Copy = seq2[seq2Beginning:totalChars + seq2Beginning]
			for j in range(0, totalChars):
				isMatch = False
				if seq1[i] == seq2Copy[j]:
					isMatch = True
				bestDistance = math.inf
				if j - 1 >= startingIndex:
					leftVal = self.bandedMatrix[i][j - 1]
					bestDistance = leftVal[0] + 5
					direction = "left"
				#TOP
				if startingIndex + j + 1 < 7:
					topVal = self.bandedMatrix[i - 1][startingIndex + j + 1] # TOP
					top = topVal[0] + 5
					if top < bestDistance:
						bestDistance = top
						direction = "down"
				#DIAG
				diagVal = self.bandedMatrix[i - 1][startingIndex + j]
				if isMatch:
					topLeft = diagVal[0] - 3
				else:
					topLeft = diagVal[0] + 1  # Going diag, so substitution(1)
				if topLeft < bestDistance:
					bestDistance = topLeft
					direction = "diag"
				self.bandedMatrix[i][startingIndex + j] = (bestDistance, direction)
			if max:
				seq2Beginning += 1
			if totalChars < k and not max:
				totalChars += 1
				if totalChars == k:
					max = True
			if i >= len(seq2) - 3:
				totalChars -= 1
			if startingIndex != 0:
				startingIndex -= 1
		lowestVal = self.bandedMatrix[len(seq1) - 1][0]
		lowestVal = lowestVal[0]
		position = 0
		if (totalChars + 1 > 7):
			return math.inf
		for i in range(1, totalChars + 1):
			temp = self.bandedMatrix[len(seq1) - 1][i]
			if temp[0] <= lowestVal:
				lowestVal = temp[0]
				position = i
		for i in range(0, 6):
			if self.bandedMatrix[len(seq1) - 1][i + 1] == 0:
				break
			position = i
		return position + 1

	def generate_unrestricted_map(self, seq1, seq2):
		align = ""
		align2 = ""
		yValue = len(seq1)
		xValue = len(seq2)
		while True:
			curVal = self.matrix[yValue][xValue]
			direction = curVal[1]
			if direction == "left":
				xValue -= 1
				align += "-"
				align2 += seq2[xValue]
			elif direction == "diag":
				xValue -= 1
				yValue -= 1
				align += seq1[yValue]
				align2 += seq2[xValue]
			else:
				yValue -= 1
				align += seq1[yValue]
				align2 += "-"
			if yValue == 0 and xValue == 0:
				break
		return align, align2

	def unrestricted_alignment(self, seq1, seq2):
		seq1 = ' ' + seq1
		seq2 = ' ' + seq2
		self.matrix = [[0 for j in range(0, len(seq2))] for i in range(0, len(seq1))]

		#Set up base case
		for i in range(0, len(seq2)):
			self.matrix[0][i] = (i * 5, "left")
		for i in range(0, len(seq1)):
			self.matrix[i][0] = (i * 5, "down")
		for i in range(1, len(seq1)):
			for j in range(1, len(seq2)):
				isMatch = False
				if seq2[j] == seq1[i]:
					isMatch = True

				direction = "left"
				leftVal = self.matrix[i][j - 1]
				bestDistance = leftVal[0] + 5  # Going right, so insertion(5)

				topVal = self.matrix[i - 1][j]
				top = topVal[0] + 5  #TOP
				if top <= bestDistance:
					bestDistance = top
					direction = "down"

				topLeftVal = self.matrix[i - 1][j - 1]
				if isMatch:
					topLeft = topLeftVal[0] - 3
				else:
					topLeft = topLeftVal[0] + 1  # Going diag, so substitution(1)
				if topLeft <= bestDistance:
					bestDistance = topLeft
					direction = "diag"
				self.matrix[i][j] = (bestDistance, direction)
		return self.matrix[len(seq1) - 1][len(seq2) - 1][0]

