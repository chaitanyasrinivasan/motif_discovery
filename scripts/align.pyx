import numpy as np
import cygibbs
cimport numpy as np
import sys

'''
Chaitanya Srinivasan

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
	#TRACEBACK
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

#Hamming distance function, returns distance
cpdef distance(seq1, seq2):
	cdef int val = 0
	cdef int i
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			val += 1
	return val

# Returns MST wrt distance of pairwise alignments
cpdef get_pairs(seqs, seqLengths):
	cdef long [:,:] distances = np.zeros(shape=(len(seqs), len(seqs)), dtype=int)
	cdef int i
	cdef int j
	cdef dict dists = dict()
	cdef int minVal
	cdef int minArg
	cdef int oddSeq = -1
	#COMPUTE PAIRWISE DISTANCES
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
	#PRIMS MST ALGORITHM 
	for i in range(len(seqs)):
		if len(dists) == (len(seqs) -1):
			oddSeq = i 
		# can only have unique pairs
		if (i not in dists):
			minVal = len(seqs[0])+1 #unobtainable value
			for j in range(len(seqs)):
				if (i != j and distances[i][j] < minVal and j not in dists):
					minVal = distances[i][j]
					minArg = j
			dists[i] = minArg
			dists[minArg] = i
	return dists, oddSeq

#Perform local alignment of two aligned sequence sets, l1 >= l2
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

	#INITIALIZE
	A[0][0] = 0
	for j in range(1,l2+1):
		A[0][j] = 0
	for i in range(1, l1+1):
		A[i][0] = 0
	#RECURRENCE
	for i in range(1, l1+1):
		for j in range(1, l2+1):
			#Check if columns completely match
			if (set(align1[:,i-1]) == set(align2[:,j-1])):
				score = match
			else:
				score = mismatch
			vals = np.array([A[i-1][j]+gap, A[i][j-1]+gap, A[i-1][j-1]+score, 0])
			A[i][j] = max(vals)
			T[i][j] = np.argmax(vals) #0, vertical, 1 horizontal, 2 diagonal
	#TRACEBACK
	i, j = np.unravel_index(np.argmax(A), np.shape(A))
	i -= 1
	j -= 1
	if (i >= j):
		offset = i-j
		aligned = np.array([[4 for _ in range(max(l1, offset+l2))] for _ in range(len(align1) + len(align2))])
		for l in range(len(align1)):
			for k in range(l1):
				aligned[l][k] = align1[l][k]
		for l in range(len(align2)):
			for k in range(l2):
				aligned[l+len(align1)][k+offset] = align2[l][k]
	else:
		offset = j-i
		aligned = np.array([[4 for _ in range(max(l1+offset, l2))] for _ in range(len(align1) + len(align2))])
		for l in range(len(align1)):
			for k in range(l1):
				aligned[l][k+offset] = align1[l][k]
		for l in range(len(align2)):
			for k in range(l2):
				aligned[l+len(align1)][k] = align2[l][k]
	return aligned

#Initalize tree with pairwise alignment of most similar sequence pairs
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

#Recursively merge alignments
cpdef merge_leaves(seqList, lengthList):
	cdef list vals
	cdef list lengths
	cdef int i

	if (len(seqList) == 1):
		return seqList.pop(0)
	else:
		#MERGE
		vals = list()
		lengths = list()
		for i in range(int(len(seqList)/2)):
			if (lengthList[i] >= lengthList[i+1]):
				vals.append(merge_align(seqList[2*i], seqList[2*i+1], lengthList[2*i], lengthList[2*i+1]))
			else:
				vals.append(merge_align(seqList[2*i+1], seqList[2*i], lengthList[2*i+1], lengthList[2*i]))
			lengths.append(len(vals[-1][0]))
		#odd number of sequences case
		if (len(seqList) % 2 == 1): 
			vals.append(seqList[-1])
			lengths.append(len(vals[-1][0]))
		return merge_leaves(vals, lengths)

cpdef entropy(A, freq):
	cdef int w = len(A[0])
	cdef double pseudo = 0.1
	cdef dict count
	cdef int j
	cdef int i
	cdef double total
	cdef int b
	cdef double entropy = 0
	cdef int gapCount

	for j in range(w):
		count = {0:pseudo, 1:pseudo, 2:pseudo, 3:pseudo}
		for i in range(len(A)):
			gapCount = 0
			if A[i][j] in count:
				count[A[i][j]] += 1
			else:
				gapCount += 1
		total = len(A) - gapCount
		for b in range(4):
			entropy += (count[b]/total)*np.log(count[b]/total/freq[b])
	return (entropy/w)

#Look for window maximizing entropy and window size
cpdef infer(multi_align, freq, seqLengths):
	cdef int i
	cdef int j
	cdef double maxVal
	cdef int maxArg
	
	for i in range(len(multi_align[0])-1):
		for j in range(i+1, len(multi_align[0])):
			val = (j-i)*entropy(multi_align[:,i:j], freq)
			if (i == 0 and j == 1):
				maxVal = val
				maxArg = j-i
			if val >= maxVal:
				maxVal = val
				maxArg = j-i
	#low homogeneity in alignment
	if maxArg < 5:
		maxArg = 10 #reasonable estimate
	for i in seqLengths:
		if maxArg > i:
			#if width too large, set it as large as possible
			maxArg = i
	return maxArg

def main(fasta):
	#0-A, 1-C, 2-G, 3-T, 4-GAP
	seqs, freq, seqLengths = cygibbs.parse(np.loadtxt(fasta, dtype="str"))
	dists, oddSeq = get_pairs(seqs, seqLengths)
	seqList, lengthList = build_roots(seqs, dists, seqLengths, oddSeq)
	multi_align = merge_leaves(seqList, lengthList)
	width = infer(multi_align, freq, seqLengths)
	print(width)

if __name__ == "__main__":
	main(sys.argv[1])
