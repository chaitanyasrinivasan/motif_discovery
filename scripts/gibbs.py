import sys
import random
import numpy as np
import concurrent.futures
import pandas as pd
import seqlogo
import matplotlib.pyplot as plt

#INPUT: .fasta file of motifs

def count(seq):
	bases = {"a":0, "c":1, "g":2, "t":3}
	freq = {0:0, 1:0, 2:0, 3:0}
	for i in range(len(seq)):
		freq[bases[seq[i]]] += 1
	return freq

def parse(path):
	seqs = np.loadtxt(path, dtype="str")
	b = 4 # number of bases
	freqs = np.empty(shape=(len(seqs), b))
	for i in range(len(seqs)):
		freq = count(seqs[i])
		for j in range(b):
			freqs[i][j] = freq[j]
	val = (freqs.sum(axis=0))
	return seqs, val/val.sum()

def propensity(A, w, k, freq):
	b = 4
	P = np.empty(shape=(b, w))
	pseudo = np.sqrt(k)*0.25
	cols = [A[:,j] for j in range(w)]
	for j in range(len(cols)):
		counts = count(cols[j])
		for i in range(b):
			P[i][j] = (((counts[i] + pseudo)/(k + (pseudo*b))))/freq[i]
	return P

#Complexity: O(kw)
def init(seqs, w, freq):
	k = len(seqs)
	z = 0
	tStar = seqs[z]
	nStar = len(seqs[z])
	index = np.array([i+1 for i in range(k-1)])
	A = np.empty(shape=(k-1, w), dtype="str")
	#Initialize alignment with random offsets
	for j in range(1, k):
		o = random.randint(0, len(seqs[j])-w-1)
		for m in range(w):
			A[j-1][m] = seqs[j][o+m]
	P =  propensity(A, w, k, freq)
	return (P, A, tStar, nStar, k, index, z)

def search(P, A, tStar, nStar, w, k, index, z, seqs, freq):
	counterList = []
	stopList = []
	counter = 0
	stop = 0
	while(counter < 50000):
		#Calculate PDF
		Pindex = {"a":0, "c":1, "g":2, "t":3, "A":0, "C": 1, "G":2, "T":3}
		pdf = np.empty(shape=nStar-w, dtype=float)
		for o in range(nStar-w):
			productNumerator = 1
			for j in range(w):
				productNumerator *= P[Pindex[tStar[o+j]]][j]
			sumDenominator = 0
			for i in range(nStar-w):
				productDenominator = 1
				for j in range(w):
					productDenominator *= P[Pindex[tStar[i+j]]][j]
				sumDenominator += productDenominator
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
		#New special sequence
		r = random.randint(0, k-2)
		#Replace new special sequence with tStar
		for i in range(w):
			A[r][i] = tStar[oStar+i]
		#Store and update parameters
		y = index[r]
		index[r] = z
		z = y
		tStar = seqs[z]
		nStar = len(seqs[z])
		Pnew = propensity(A, w, k, freq)
		counterList.append(counter)
		stopList.append(stop)
		if np.allclose(P, Pnew):
			stop += 1
		else:
			stop = 0
		P = Pnew
		counter += 1

	#Obtain A by one last update
	pdf = np.empty(shape=nStar-w, dtype=float)
	for o in range(nStar-w):
		productNumerator = 1
		for j in range(w):
			productNumerator *= P[Pindex[tStar[o+j]]][j]
		sumDenominator = 0
		for i in range(nStar-w):
			productDenominator = 1
			for j in range(w):
				productDenominator *= P[Pindex[tStar[i+j]]][j]
			sumDenominator += productDenominator
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
	Anew = np.empty(shape=(k, w), dtype="str")
	for i in range(len(A)):
		for j in range(w):
			Anew[i][j] = A[i][j]
	for j in range(w):
		Anew[k-1][j] = tStar[oStar+j]
	S = np.log2(propensity(Anew, w, k, freq))
	return S, Anew, counterList, stopList

def plotLogo(fasta, alignment, w):
	pwm = np.empty(shape=(w, 4))
	pseudo = 0
	for j in range(w):
		count = {"a":pseudo, "c":pseudo, "g":pseudo, "t":pseudo}
		for i in range(len(alignment)):
			count[alignment[i][j]] += 1
		total = sum(count.values())
		pwm[j][0] = count["a"]/total
		pwm[j][1] = count["c"]/total
		pwm[j][2] = count["g"]/total
		pwm[j][3] = count["t"]/total
	pwmFormatted = seqlogo.CompletePm(pwm)
	seqlogo.seqlogo(pwmFormatted, ic_scale=False, format="png", size="medium", filename=str(fasta)[:-4]+"_motif.png")

def main(fasta, size):
	w = int(size) #Length of Motif
	seqs, freq = parse(fasta)
	P, A, tStar, nStar, k, index, z = init(seqs, w, freq)
	S, Anew, counterList, stopList = search(P, A, tStar, nStar, w, k, index, z, seqs, freq)
	plotLogo(fasta, Anew, w)
	fig = plt.figure()
	plt.plot(counterList, stopList)
	plt.savefig("loss.png")



if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])
