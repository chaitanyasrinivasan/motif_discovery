import cygibbs
import sys
import numpy as np
import pickle
import argparse

# read and merge inputs from partitioned run
def merge(jobs, fasta, k, w):
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
	# store sequence lengths
	for i in range(len(seqs)):
		seqLengths[i] = len(seqs[i])
		for j in range(len(seqs[i])):
			intSeqs[i][j] = bases[seqs[i][j]]
	return A, {key:val/(k) for key,val in freq.items()}, seqLengths, intSeqs

def main(jobs, fasta, k, w):
	print("Merging")
	A, freq, seqLengths, seqs = merge(jobs, fasta, k, w)
	# initialize parameters for motif discovery
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
	parser = argparse.ArgumentParser()
	parser.add_argument("-j", "--jobs", type=str, help="File path to list of jobs")
	parser.add_argument("-f", "--fasta", type=str, help="File path to input fasta")
	parser.add_argument("-k", "--size", type=str, help="Number of sequences (integer)")
	parser.add_argument("-w", "--width", type=str, help="Motif width (integer)")

	options, args = parser.parse_known_args()

	if (len(sys.argv)==1):
	    parser.print_help(sys.stderr)
	    sys.exit(1)
	elif (options.input is None):
		parser.print_help(sys.stderr)
		sys.exit(1)
	elif (options.fasta is None):
		parser.print_help(sys.stderr)
		sys.exit(1)
	elif (options.width is None):
		parser.print_help(sys.stderr)
		sys.exit(1)
	elif (int(options.width) < 1):
		raise Exception("Width must be an integer greater than 0")
		sys.exit(1)
	elif (options.size is None):
		parser.print_help(sys.stderr)
		sys.exit(1)
	elif (int(options.size) <= 2):
		raise Exception("Number of sequences must be an integer greater than 2")
		sys.exit(1)
	else:
		main(options.jobs, options.fasta, int(options.size), int(options.width))
