import numpy as np
import BW
import sys


#Takes in two sequences and their length
#Returns the pairwise alignment
def pairwise_align(seq1, seq2, l1, l2):
	match, mismatch, gap = 5, -1, -4
	A = np.empty(shape=(l1+1, l2+1), dtype=int)
	T = np.zeros(shape=(l1+1, l2+1), dtype=int)
	#INITIALIZE
	A[0][0] = gap
	for j in range(1,l2+1):
		A[0][j] = A[0][j-1] + gap
	for i in range(1, l1+1):
		A[i][0] = A[i-1][0] + gap
	#RECURRENCE
	for i in range(1, l1+1):
		for j in range(1, l2+1):
			if seq1[i-1] == seq2[j-1]:
				score = match
			else:
				score = mismatch
			vals = np.array([A[i-1][j]+gap, A[i][j-1]+gap, A[i-1][j-1]+score])
			A[i][j] = max(vals)
			T[i][j] = np.argmax(vals) #0, vertical, 1 horizontal, 2 diagonal
	#TRACEBACK
	aligned = np.array([[4 for _ in range(max(l1, l2))] for _ in range(2)])
	i, j = l1-1, l2-1
	for k in range(max(l1, l2)-1, -1, -1):
		if (T[i][j]) == 0:
			aligned[0][k] = seq1[i]
			aligned[1][k] = 4
			i -= 1
		elif (T[i][j]) == 1:
			aligned[0][k] = 4
			aligned[1][k] = seq2[j]
			j -= 1
		else:
			aligned[0][k] = seq1[i]
			aligned[1][k] = seq2[j]
			i -= 1
			j -= 1
	return aligned

#Hamming distance function
def distance(seq1, seq2):
	val = 0
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			val += 1
	return val

#Pairwise distance matrix of alignments
def get_pairs(seqs, seqLengths):
	distances = np.zeros(shape=(len(seqs), len(seqs)), dtype=int)
	for i in range(len(seqs)):
		for j in range(len(seqs)):
			if (i != j and i > j):
				aligned = pairwise_align(seqs[i], seqs[j], seqLengths[i], seqLengths[j])
				distanceVal = distance(aligned[0], aligned[1])
				distances[i][j] = distanceVal
				distances[j][i] = distanceVal
	return distances

#Perform alignment of two aligned sequence sets
def merge_align(align1, align2, l1, l2):
	match, mismatch, gap = 5, -1, -4
	A = np.empty(shape=(l1+1, l2+1), dtype=int)
	T = np.zeros(shape=(l1+1, l2+1), dtype=int)
	#INITIALIZE
	A[0][0] = gap
	for j in range(1,l2+1):
		A[0][j] = A[0][j-1] + gap
	for i in range(1, l1+1):
		A[i][0] = A[i-1][0] + gap
	#RECURRENCE
	for i in range(1, l1+1):
		for j in range(1, l2+1):
			#Check if columns completely match
			if (set(align1[:,i-1]) == set(align2[:,j-1])):
				if (4 in align1[:,i-1] or 4 in align2[:,j-1]):
					score = match + gap
				else:
					score = match
			else:
				score = mismatch
			vals = np.array([A[i-1][j]+gap, A[i][j-1]+gap, A[i-1][j-1]+score])
			A[i][j] = max(vals)
			T[i][j] = np.argmax(vals) #0, vertical, 1 horizontal, 2 diagonal
	#TRACEBACK
	aligned = np.array([[4 for _ in range(max(l1, l2))] for _ in range(len(align1)+len(align2))])
	i, j = l1-1, l2-1
	for k in range(max(l1, l2)-1, -1, -1):
		if (T[i][j] == 0):
			for l in range(len(align1)):
				aligned[l][k] = align1[l][i]
			for m in range(len(align1)-1, len(aligned)):
				aligned[m][k] = 4
			i -= 1
		elif (T[i][j] == 1):
			for l in range(len(align1)):
				aligned[l][k] = 4
			for m in range(len(align2)):
				aligned[len(align1)+m][k] = align2[m][j]
			j -= 1
		else:
			for l in range(len(align1)):
				aligned[l][k] = align1[l][i]
			for m in range(len(align2)):
				aligned[len(align1)+m][k] = align2[m][j]
			i -= 1
			j -= 1
	return aligned

def result(seqs, seqLengths, distances):
	if (len(seqs) == 1):
		return seqs.pop(0)
	else:
		vals = list()
		lengths = list()
		for i in range(int(len(seqs)/2)):
			vals.append(merge_align(seqs[i], seqs[i+1], seqLengths[i], seqLengths[i+1]))
			lengths.append(len(vals[-1][0]))
		return result(vals, lengths, distances)

if __name__ == "__main__":
	#0-A, 1-C, 2-G, 3-T, 4-GAP
	seqs, seqLengths = BW.parse(sys.argv[1])
	distances = get_pairs(seqs, seqLengths)
	multi_align = result([pairwise_align(seqs[i], seqs[i+1], seqLengths[i], seqLengths[i+1]) for i in range(int(len(seqs)/2))], [y for y in seqLengths], distances)
	print(multi_align)
	np.savetxt(sys.argv[1][:-4]+".align", multi_align)
