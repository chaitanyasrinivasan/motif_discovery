import numpy as np
import cygibbs
cimport numpy as np
import sys

'''
This program performs a polytime progressive alignment by iteratively merging
pairwise local alignments of sequences.
'''

#Takes in two sequences and their lengths, l1 >= l2
#Returns the pairwise local alignment
cpdef pairwise_align(seq1, seq2, l1, l2):
	cdef int match = 5
	cdef int mismatch = -1
	cdef int gap = -4 
	cdef long [:,:] A = np.empty(shape=(l1+1, l2+1), dtype=int)
	cdef long [:,:] T = np.zeros(shape=(l1+1, l2+1), dtype=int)
	cdef int i
	cdef int j
	cdef int score
	cdef long [:] vals
	cdef long [:,:] aligned
	cdef str str1 = ""
	cdef str str2 = ""
	cdef int starti
	cdef int startj
	#INITIALIZE ALIGNMENT MATRIX
	A[0][0] = 0
	for j in range(1,l2+1):
		A[0][j] = 0
	for i in range(1, l1+1):
		A[i][0] = 0
	#GET OPTIMAL ALIGNMENT SCORE
	for i in range(1, l1+1):
		for j in range(1, l2+1):
			if seq1[i-1] == seq2[j-1]:
				score = match
			else:
				score = mismatch
			vals = np.array([A[i-1][j]+gap, A[i][j-1]+gap, A[i-1][j-1]+score, 0])
			A[i][j] = np.amax(vals)
			T[i][j] = np.argmax(vals) #0, vertical, 1 horizontal, 2 diagonal
	#Get coordinates of max value in matrix
	i, j = np.unravel_index(np.argmax(A), np.shape(A))
	#readjust indices to read from seq1, seq2
	i -= 1
	j -= 1
	if (i >= j):
		offset = i-j
		aligned = np.array([[4 for _ in range(max(l1, offset+l2))] for _ in range(2)])
		for k in range(l1):
			aligned[0][k] = seq1[k]
		for k in range(l2):
			aligned[1][k+offset] = seq2[k]
	else:
		offset = j-i
		aligned = np.array([[4 for _ in range(max(l1+offset, l2))] for _ in range(2)])
		for k in range(l1):
			aligned[0][k+offset] = seq1[k]
		for k in range(l2):
			aligned[1][k] = seq2[k]
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
	cdef dict dists = dict()
	cdef int minVal
	cdef int minArg
	cdef int oddSeq = -1
	#COMPUTE DISTANCES
	for i in range(len(seqs)):
		for j in range(len(seqs)):
			if (i > j):
				if (seqLengths[i] >= seqLengths[j]):
					aligned = pairwise_align(seqs[i], seqs[j], seqLengths[i], seqLengths[j])
				else:
					aligned = pairwise_align(seqs[j], seqs[i], seqLengths[j], seqLengths[i])
				distanceVal = distance(aligned[0], aligned[1])
				distances[i][j] = distanceVal
				distances[j][i] = distanceVal
	#PRIMS ALGORITHM
	for i in range(len(seqs)):
		if len(dists) == (len(seqs) -1):
			oddSeq = i
		if (i not in dists):
			minVal = len(seqs[0])+1 #unobtainable value
			for j in range(len(seqs)):
				if (i != j and distances[i][j] < minVal and j not in dists):
					minVal = distances[i][j]
					minArg = j
			dists[i] = minArg
			dists[minArg] = i
	return dists, oddSeq

#Perform alignment of two aligned sequence sets
cpdef merge_align(align1, align2, l1, l2):
	cdef int match = 5
	cdef int mismatch = -1
	cdef int gap = -4
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
	cdef set common = set()
	#INITIALIZE
	A[0][0] = 0
	for j in range(1,l2+1):
		A[0][j] = A[0][j-1] + gap
	for i in range(1, l1+1):
		A[i][0] = 0
	#RECURRENCE
	for i in range(1, l1+1):
		for j in range(1, l2+1):
			#Check if columns completely match
			if (set(align1[:,i-1]) == set(align2[:,j-1]) and set(align1[:,i-1]) != set([4])):
				score = match
			else:
				score = mismatch
			vals = np.array([A[i-1][j]+gap, A[i][j-1]+gap, A[i-1][j-1]+score])
			A[i][j] = max(vals)
			T[i][j] = np.argmax(vals) #0, vertical, 1 horizontal, 2 diagonal
	#TRACEBACK
	i, j = np.argmax(A[:,l2]), l2
	aligned = np.array([[4 for _ in range(l1+(l2-i+1))] for _ in range(len(align1)+len(align2))])
	for l in range(len(align1)):
		for m in range(i-1, l1):
			aligned[l][m+(l2-i+1)] = align1[l][m]
	while (i > 0):
		if (T[i][j] == 0):
			for m in range(len(align2)):
				aligned[m+len(align1)][i-1] = align2[m][j-1]
			i -= 1
		elif (T[i][j] == 1):
			for m in range(len(align1)):
				aligned[m][i-1] = align1[m][i-1]
			j -= 1
		else:
			for m in range(len(align1)):
				aligned[m][i-1] = align1[m][i-1]
			for m in range(len(align2)):
				aligned[m+len(align1)][i-1] = align2[m][j-1]
			i -= 1
			j -= 1
	return aligned

cpdef build_roots(seqs, dists, seqLengths, oddSeq):
	noDoubles = set()
	seqList = list()
	lengthList = list()
	for key in dists:
		if (key not in noDoubles):
			if seqLengths[key] >= seqLengths[dists[key]]:
				seqList.append(pairwise_align(seqs[key], seqs[dists[key]], seqLengths[key], seqLengths[dists[key]]))
			else:
				seqList.append(pairwise_align(seqs[dists[key]], seqs[key], seqLengths[dists[key]], seqLengths[key]))
			lengthList.append(len(seqList[-1][0]))
			noDoubles.add(key)
			noDoubles.add(dists[key])
	if (oddSeq != -1):
		seqList.append(seqs[oddSeq])
	return seqList, lengthList

cpdef merge_leaves(seqList, lengthList):
	cdef list vals
	cdef list lengths
	cdef int i
	if (len(seqList) == 1):
		return seqList.pop(0)
	else:
		vals = list()
		lengths = list()
		for i in range(int(len(seqList)/2)):
			if (lengthList[i] >= lengthList[i+1]):
				vals.append(merge_align(seqList[2*i], seqList[2*i+1], lengthList[2*i], lengthList[2*i+1]))
			else:
				vals.append(merge_align(seqList[2*i+1], seqList[2*i], lengthList[2*i+1], lengthList[2*i]))
			lengths.append(len(vals[-1][0]))
		if (len(seqList) % 2 == 1): #add merge remaining sequence
			vals.append(seqList[-1])
			lengths.append(len(vals[-1][0]))
		return merge_leaves(vals, lengths)

def main(fasta):
	#0-A, 1-C, 2-G, 3-T, 4-GAP
	seqs, freq, seqLengths = cygibbs.parse(np.loadtxt(fasta, dtype="str"))
	dists, oddSeq = get_pairs(seqs, seqLengths)
	seqList, lengthList = build_roots(seqs, dists, seqLengths, oddSeq)
	multi_align = merge_leaves(seqList, lengthList)
	print(np.array([x for x in multi_align]))
	#return multi_align

if __name__ == "__main__":
	main(sys.argv[1])
