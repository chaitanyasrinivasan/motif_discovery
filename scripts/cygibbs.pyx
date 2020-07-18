import sys
import random
import numpy as np
cimport numpy as np
import pandas as pd
import seqlogo
import matplotlib.pyplot as plt
import pickle
import argparse

#Input: seq - np string array of a fasta sequence
#Output: freq - a dictionary of column-wise frequencies
#Output: intList - a list of quaternarized nucleotide sequences
cdef count(seq):
	cdef dict bases = {"a":0, "c":1, "g":2, "t":3}
	cdef dict freq = {0:0, 1:0, 2:0, 3:0}
	cdef int i
	cdef long[:] intList = np.array([0 for i in range(len(seq))])
	for i in range(len(seq)):
		freq[bases[seq[i]]] += 1
		#Convert string to int for efficiency
		intList[i] = int(bases[seq[i]])
	return (freq, intList)

#Input: stringSeqs - np string matrix of fasta sequences
#Output: quarternarized sequence matrix, background frequency dictionary, 
#			length of sequences (cannot be inferred without O(|seq|) work
#								since all rows equal dimensions)
cdef parse(stringSeqs):
	cdef int maxSize = len(max(stringSeqs, key=len))
	cdef int b = 4 # number of bases
	cdef int i
	cdef int j
	cdef int k = len(stringSeqs)
	cdef double[:] val
	cdef double[:,:] freqs = np.empty(shape=(k, b))
	cdef long[:,:] intSeqs = np.array([[4 for j in range(maxSize)] for i in range(k)])
	cdef long[:] seqLengths = np.empty(shape=(k), dtype="int")
	#0 = a, 1=c, 2=g, 3=t, 4=no base, space filler to maintain dimensions
	for i in range(k):
		freq, intList = count(stringSeqs[i])
		seqLengths[i] = len(stringSeqs[i])
		for j in range(len(intList)):
			intSeqs[i][j] = intList[j]
		for j in range(b):
			freqs[i][j] = freq[j]
	val = (np.sum(freqs, axis=0))
	return np.array(intSeqs), val/np.sum(val), seqLengths

cpdef propensity(A, int w, k, freq):
	cdef int b = 4
	cdef double[:,:] P = np.empty(shape=(b, w))
	cdef double pseudo = np.sqrt(k)*0.25
	cdef int j
	cdef np.ndarray cols = np.array([A[:,j] for j in range(w)])
	cdef int i
	cdef dict counts
	for j in range(len(cols)):
		counts = {0:0, 1:0, 2:0, 3:0}
		for i in range(len(cols[j])):
			if int(cols[j][i]) in counts:
				counts[int(cols[j][i])] += 1
		for i in range(b):
			P[i][j] = (((counts[i] + pseudo)/(k + (pseudo*b))))/freq[i]
	return P

#Complexity: O(kw)
cdef init(seqs, w, freq, seqLengths):
	cdef int k = len(seqs)
	cdef int z = 0
	cdef long[:] tStar = seqs[z]
	cdef int nStar = seqLengths[z]
	cdef int i
	cdef long[:] index = np.array([int(i+1) for i in range(k-1)])
	cdef long[:,:] A = np.empty(shape=(k-1, w), dtype="int")
	cdef int j
	cdef int o
	cdef int m
	#Initialize alignment with random offsets
	for j in range(1, k):
		o = random.randint(0, seqLengths[j]-w-1)
		for m in range(w):
			A[j-1][m] = seqs[j][o+m]
	P =  propensity(A, w, k, freq)
	return (P, A, tStar, nStar, k, index, z)

cdef entropy(A, freq):
	cdef int w = len(A[0])
	cdef double pseudo = 0.1
	cdef dict count
	cdef int j
	cdef int i
	cdef double total
	cdef int b
	cdef double entropy = 0
	for j in range(w):
		count = {0:pseudo, 1:pseudo, 2:pseudo, 3:pseudo}
		for i in range(len(A)):
			if A[i][j] in count:
				count[A[i][j]] += 1
		total = len(A)
		for b in range(4):
			entropy += (count[b]/total)*np.log(count[b]/total/freq[b])
	return entropy/w

cpdef getOptSeq(P, tStar, nStar, w):
	cdef double[:] pdf
	cdef int o
	cdef double productNumerator
	cdef int j
	cdef double sumDenominator
	cdef int i
	cdef double productDenominator
	cdef double[:] cdf
	cdef double sumCDF
	cdef int l
	cdef double oRand
	cdef int oStar

	#Calculate PDF
	pdf = np.empty(shape=nStar-w, dtype=float)
	sumDenominator = 0
	for i in range(nStar-w):
		productDenominator = 1
		for j in range(w):
			productDenominator *= P[int(tStar[i+j])][j]
		sumDenominator += productDenominator
	for o in range(nStar-w):
		productNumerator = 1
		for j in range(w):
			productNumerator *= P[int(tStar[o+j])][j]
		pdf[o] = productNumerator/sumDenominator
	#Calculate CDF
	cdf = np.empty(shape=nStar-w, dtype=float)
	for o in range(len(cdf)-1):
		sumCDF = 0
		for l in range(o+1):
			sumCDF += pdf[l]
		cdf[o] = sumCDF
	cdf[len(cdf)-1] = 1
	oRand = random.random()
	oStar = -1
	for i in range(len(cdf)-1):
		if oRand <= cdf[i]:
			oStar = i
			break
	if oStar == -1: oStar = len(cdf)-1
	return oStar


cpdef search(P, A, tStar, nStar, w, k, index, z, seqs, freq, seqLengths, maxIter):
	cdef int r
	cdef int y
	cdef int c
	cdef long[:,:] Anew
	cdef double[:,:] S
	cdef double currLoss = entropy(A, freq)
	cdef double maxLoss = currLoss
	cdef list loss = [currLoss]

	#Store initial state
	bestA = A
	bestNStar = nStar
	bestTStar = tStar
	bestIndex = index
	#Run through maxIter iterations looking for max entropy
	for c in range(maxIter):
		oStar = getOptSeq(P, tStar, nStar, w)
		#New special sequence
		r = random.randint(0, k-2)
		#Replace new special sequence with tStar
		for i in range(w):
			A[r][i] = tStar[oStar+i]
		currLoss = entropy(A, freq)
		loss.append(currLoss)
		#Store and update parameters
		y = index[r]
		index[r] = z
		z = y
		tStar = seqs[z]
		nStar = seqLengths[z]
		P = propensity(A, w, k, freq)
		#Update max args
		if (currLoss > maxLoss):
			maxLoss = currLoss
			bestA = A
			bestNStar = nStar
			bestTStar = tStar
			bestIndex = index
	
	#Obtain A by one last update
	A = bestA
	nStar = bestNStar
	tStar = bestTStar
	index = bestIndex
	P = propensity(A, w, k, freq)
	oStar = getOptSeq(P, tStar, nStar, w)
	Anew = np.empty(shape=(k, w), dtype="int")
	for i in range(len(A)):
		for j in range(w):
			Anew[i][j] = A[i][j]
	for j in range(w):
		Anew[k-1][j] = int(tStar[oStar+j])
	S = np.log2(propensity(Anew, w, k, freq))
	return S, Anew, loss

def plotLogo(fasta, alignment, w):
	pwm = np.empty(shape=(w, 4))
	pseudo = 0
	for j in range(w):
		count = {0:pseudo, 1:pseudo, 2:pseudo, 3:pseudo}
		for i in range(len(alignment)):
			count[alignment[i][j]] += 1
		total = sum(count.values())
		for b in range(4):
			pwm[j][b] = count[b]/total
	pwmFormatted = seqlogo.CompletePm(pwm)
	seqlogo.seqlogo(pwmFormatted, ic_scale=False, format="png", size="medium", filename=str(fasta)[:-4]+"_motif.png")
	print("Motif logo written to " + str(fasta)[:-4]+"_motif.png")

def plotLoss(fasta, loss):
	iters = np.arange(len(loss))
	loss = np.array(loss)
	fig, ax = plt.subplots()
	ax.plot(iters, loss)
	ax.set(xlabel='Iterations', ylabel='Entropy')
	fig.savefig(str(fasta)[:-4]+"_loss.png")
	print("Entropy function written to " + str(fasta)[:-4]+"_loss.png")


def main(fasta, size, mode):
	w = int(size) #Length of Motif
	print("Parsing from " + fasta)
	seqs, freq, seqLengths = parse(np.loadtxt(fasta, dtype="str"))
	if (mode == "divide"):
		with open(str(fasta)[:-4]+".freq", "wb") as f:
			pickle.dump(freq, f)
			print("Dictionary written to "+str(fasta)[:-4]+".freq")
	print("Initializing")
	P, A, tStar, nStar, k, index, z = init(seqs, w, freq, seqLengths)
	print("Searching")
	maxIter = k*k*np.sum(seqLengths)/len(seqLengths)
	S, Anew, loss = search(P, A, tStar, nStar, w, k, index, z, seqs, freq, seqLengths, maxIter)
	plotLogo(fasta, Anew, w)
	plotLoss(fasta, loss)
	if (mode == "divide"):
		np.savetxt(str(fasta)[:-4]+".matrix", Anew)
		print("Alignment matrix written to "+str(fasta)[:-4]+".matrix")
	print("Done")


if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3])
