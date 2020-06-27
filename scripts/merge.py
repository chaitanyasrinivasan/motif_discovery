import cygibbs
import sys
import numpy as np
import pickle

def merge(jobs, fasta, k, w):
	#Read and merge inputs from partitioned run
	A = np.empty(shape=(k, w), dtype="int")
	seqLengths = np.empty(shape=k, dtype="int")
	freq = {0:0, 1:0, 2:0, 3:0}
	counter = 0
	with open(jobs, "r") as f:
		for line in f.read().splitlines():
			data = np.loadtxt(str(line)[:-4]+".matrix")
			with open (line[:-4]+".freq", "rb") as f:
				splitFreq = pickle.loads(f.read())
			for i in range(len(data)):
				for j in range(w):
					A[counter][j] = int(data[i][j])
					if j in freq:
						freq[j] += splitFreq[j]
				counter += 1
	seqs = np.loadtxt(fasta[:-4]+"_shuffled.txt", dtype="str")
	maxSize = len(max(seqs, key=len))
	intSeqs = np.array([[4 for j in range(maxSize)] for i in range(k)])
	bases = {"a":0, "c":1, "g":2, "t":3}
	for i in range(len(seqs)):
		seqLengths[i] = len(seqs[i])
		for j in range(len(seqs[i])):
			intSeqs[i][j] = bases[seqs[i][j]]
	return A, {key:val/(k) for key,val in freq.items()}, seqLengths, intSeqs

def main(jobs, fasta, k, w):
	print("Merging")
	A, freq, seqLengths, seqs = merge(jobs, fasta, k, w)
	P = np.array(cygibbs.propensity(A, w, k, freq))
	z = 0
	tStar = seqs[0]
	nStar = seqLengths[z]
	index = np.array([int(i+1) for i in range(k-1)])
	maxIter = k*np.sum(seqLengths)/len(seqLengths)
	print("Searching")
	S, Anew, loss = cygibbs.search(P, A, tStar, nStar, w, k, index, z, seqs, freq, seqLengths, maxIter)
	cygibbs.plotLogo(fasta, Anew, w)
	cygibbs.plotLoss(fasta, loss)
	print("Done")


if __name__ == "__main__":
	jobs = sys.argv[1]
	fasta = sys.argv[2]
	k = int(sys.argv[3])
	w = int(sys.argv[4])
	if not isinstance(w, int):
		raise Exception("Width is not an integer, exiting program.")
		sys.exit()
	main(jobs, fasta, k, w)
