import numpy as np
cimport numpy as np
import BW
import sys


#Takes in two sequences and their length, l1 >= l2
#Returns the pairwise semi-global alignment
cpdef pairwise_align(seq1, seq2, l1, l2):
	cdef int match = 5
	cdef int mismatch = -1
	cdef int gap = -4
	cdef long [:,:] A = np.empty(shape=(l1+1, l2+1), dtype=int)
	cdef long [:,:] T = np.zeros(shape=(l1+1, l2+1), dtype=int)
	cdef int i
	cdef int j
	cdef int k
	cdef int score
	cdef long [:] vals
	cdef long [:,:] aligned
	#INITIALIZE
	A[0][0] = 0
	for j in range(1,l2+1):
		A[0][j] = A[0][j-1] + gap
	for i in range(1, l1+1):
		A[i][0] = 0
	#RECURRENCE
	for i in range(1, l1+1):
		for j in range(1, l2+1):
			if seq1[i-1] == seq2[j-1]:
				score = match
			else:
				score = mismatch
			vals = np.array([A[i-1][j]+gap, A[i][j-1]+gap, A[i-1][j-1]+score])
			A[i][j] = np.amax(vals)
			T[i][j] = np.argmax(vals) #0, vertical, 1 horizontal, 2 diagonal
	print(np.array([x for x in A]))
	print(np.array([x for x in T]))
	#TRACEBACK
	aligned = np.array([[4 for _ in range(l1)] for _ in range(2)])
	for i in range(l1):
		aligned[0][i] = seq1[i]
	i, j = np.argmax(A[:,l2]), l2
	while (i > 0):
		if (T[i][j]) == 0:
			aligned[1][i-1] = 4
			i -= 1
		elif (T[i][j]) == 1:
			aligned[0][i-1] = 4
			j -= 1
		else:
			aligned[1][i-1] = seq2[j-1]
			i -= 1
			j -= 1
	print(np.array([x for x in aligned]))
	return aligned

#Hamming distance function
cpdef distance(seq1, seq2):
	cdef int val = 0
	cdef int i
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			val += 1
	return val

#Pairwise distance matrix of alignments
cpdef get_pairs(seqs, seqLengths):
	cdef long [:,:] distances = np.zeros(shape=(len(seqs), len(seqs)), dtype=int)
	cdef int i
	cdef int j
	for i in range(len(seqs)):
		for j in range(len(seqs)):
			if (i != j and i > j):
				aligned = pairwise_align(seqs[i], seqs[j], seqLengths[i], seqLengths[j])
				distanceVal = distance(aligned[0], aligned[1])
				distances[i][j] = distanceVal
				distances[j][i] = distanceVal
	return distances

#Perform alignment of two aligned sequence sets
cpdef merge_align(align1, align2, l1, l2):
	cdef int match = 5
	cdef int mismatch = -1
	cdef int gap = -2
	cdef long [:,:] A = np.empty(shape=(l1+1, l2+1), dtype=int)
	cdef long [:,:] T = np.zeros(shape=(l1+1, l2+1), dtype=int)
	cdef int i
	cdef int j
	cdef int k
	cdef int l
	cdef int m
	cdef int score
	cdef long [:] vals
	cdef long [:,:] aligned

	#INITIALIZE
	A[0][0] = 0
	for j in range(1,l2+1):
		if l1 >= l2:
			A[0][j] = A[0][j-1] + gap
		else:
			A[0][j] = 0
	for i in range(1, l1+1):
		if l2 > l1:
			A[i][0] = A[i-1][0] + gap
		else:
			A[i][0] = 0
	#RECURRENCE
	for i in range(1, l1+1):
		for j in range(1, l2+1):
			#Check if columns completely match
			if (set(align1[:,i-1]) == set(align2[:,j-1])):
				if (4 in align1[:,i-1] or 4 in align2[:,j-1]):
					score = match + gap
				else:
					score = match + match
			elif (len(set.intersection(set(align1[:,i-1]), set(align2[:,j-1]))) > 1):
				score = match
			else:
				score = mismatch
			vals = np.array([A[i-1][j]+gap, A[i][j-1]+gap, A[i-1][j-1]+score])
			A[i][j] = max(vals)
			T[i][j] = np.argmax(vals) #0, vertical, 1 horizontal, 2 diagonal
	#TRACEBACK
	aligned = np.array([[4 for _ in range(max(l1, l2))] for _ in range(len(align1)+len(align2))])
	if l1 >= l2:
		i, j = np.argmax(A[:,l2]), l2
	else:
		i, j = l1, np.argmax(A[l1])
	for k in range(max(l1, l2)-1, -1, -1):
		if (T[i][j] == 0):
			for l in range(len(align1)):
				aligned[l][k] = align1[l][i-1]
			for m in range(len(align1)-1, len(aligned)):
				aligned[m][k] = 4
			i -= 1
		elif (T[i][j] == 1):
			for l in range(len(align1)):
				aligned[l][k] = 4
			for m in range(len(align2)):
				aligned[len(align1)+m][k] = align2[m][j-1]
			j -= 1
		else:
			for l in range(len(align1)):
				aligned[l][k] = align1[l][i-1]
			for m in range(len(align2)):
				aligned[len(align1)+m][k] = align2[m][j-1]
			i -= 1
			j -= 1
	print(np.array([x for x in aligned]))
	return aligned

cpdef result(seqs, seqLengths, distances):
	cdef list vals
	cdef list lengths
	cdef int i
	if (len(seqs) == 1):
		return seqs.pop(0)
	else:
		vals = list()
		lengths = list()
		for i in range(int(len(seqs)/2)):
			vals.append(merge_align(seqs[i], seqs[i+1], seqLengths[i], seqLengths[i+1]))
			lengths.append(len(vals[-1][0]))
		return result(vals, lengths, distances)

def main(fasta):
	#0-A, 1-C, 2-G, 3-T, 4-GAP
	seqs, seqLengths = BW.parse(fasta)
	distances = get_pairs(seqs, seqLengths)
	seqList = list()
	for i in range(int(len(seqs)/2)):
		if seqLengths[i] >= seqLengths[i+1]:
			seqList.append(pairwise_align(seqs[i], seqs[i+1], seqLengths[i], seqLengths[i+1]))
		else:
			seqList.append(pairwise_align(seqs[i+1], seqs[i], seqLengths[i+1], seqLengths[i]))
	multi_align = result(seqList, [y for y in seqLengths], distances)
	return multi_align

if __name__ == "__main__":
	main(sys.argv[1])

