import numpy as np
cimport numpy as np
import BW
import sys


#Takes in two sequences and their length, l1 >= l2
#Returns the pairwise local alignment
cpdef pairwise_align(seq1, seq2, l1, l2):
	print(np.array([x for x in seq1]))
	print(np.array([x for x in seq2]))
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
	#INITIALIZE
	A[0][0] = 0
	for j in range(1,l2+1):
		A[0][j] = 0
	for i in range(1, l1+1):
		A[i][0] = 0
	#RECURRENCE
	for i in range(1, l1+1):
		for j in range(1, l2+1):
			if seq1[i-1] == seq2[j-1]:
				score = match
			else:
				score = mismatch
			vals = np.array([A[i-1][j]+gap, A[i][j-1]+gap, A[i-1][j-1]+score, 0])
			A[i][j] = np.amax(vals)
			T[i][j] = np.argmax(vals) #0, vertical, 1 horizontal, 2 diagonal
	#TRACEBACK LOCAL ALIGNMENT
	print(np.array([x for x in A]))		
	print(np.array([x for x in T]))
	#Get coordinates of max value in matrix
	starti, startj = np.unravel_index(np.argmax(A), np.shape(A))
	i, j = starti, startj
	print(i, j)
	while (A[i][j] > 0):
		if (T[i][j]) == 0:
			str2 += str(4)
			i -= 1
		elif (T[i][j]) == 1:
			str1 += str(4)
			j -= 1
		else:
			str1 += str(seq1[i-1])
			str2 += str(seq2[j-1])
			i -= 1
			j -= 1
	aligned = np.array([[4 for _ in range(i+len(str1)+l1-starti)] for _ in range(2)])
	print(seq1[:l1-starti])
	print(str1[::-1])
	print(seq1[i:])
	print(seq2[:startj])
	print(str2[::-1])
	print(seq2[j:])
	print(np.array([x for x in aligned]))
	sys.exit()
	#return aligned

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
		if (i in dists):
			break
		minVal = len(seqs[0])+1
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
	print(np.array([x for x in A]))
	print(np.array([x for x in T]))
	#TRACEBACK

	i, j = np.argmax(A[:,l2]), l2
	aligned = np.array([[4 for _ in range(l1+(l2-i+1))] for _ in range(len(align1)+len(align2))])
	for l in range(len(align1)):
		for m in range(i-1, l1):
			aligned[l][m+(l2-i+1)] = align1[l][m]
	while (i > 0):
		print(i, j)
		print(np.array([x for x in aligned]))
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
	print(np.array([x for x in align1]))
	print(np.array([x for x in align2]))
	print(np.array([x for x in aligned]))
	sys.exit()
	return aligned

cpdef result(seqList, lengthList):
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
				vals.append(merge_align(seqList[i], seqList[i+1], lengthList[i], lengthList[i+1]))
			else:
				vals.append(merge_align(seqList[i+1], seqList[i], lengthList[i+1], lengthList[i]))
			lengths.append(len(vals[-1][0]))
		return result(vals, lengths)

def main(fasta):
	#0-A, 1-C, 2-G, 3-T, 4-GAP
	seqs, seqLengths = BW.parse(fasta)
	dists, oddSeq = get_pairs(seqs, seqLengths)
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
	multi_align = result(seqList, lengthList)
	print(np.array([x for x in multi_align]))
	#return multi_align

if __name__ == "__main__":
	main(sys.argv[1])
