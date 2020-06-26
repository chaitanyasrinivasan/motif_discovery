import cygibbs
import sys
import numpy as np
import json

def merge(jobs, k, w):
	A = np.empty(shape=(k, w), dtype="int")
	freq = {0:0, 1:0, 2:0, 3:0}
	with open(jobs, "r") as f:
		for counter, line in enumerate(f.read().splitlines(), start=0):
			data = np.loadtxt(str(line)[:-4]+".matrix")
			splitFreq = json.loads(line[:-4]+".freq")
			for i in range(len(data)):
				for j in range(4):
					A[counter][j] = int(data[i][j])
					freq[j] += splitFreq[j]
	return A, {k:v/(counter+1) for k,v in freq.iteritems()}



def main(jobs, motif, k, w):
	A, freq = merge(jobs, k, w)
	print(A)
	print(freq)



if __name__ == "__main__":
	jobs = sys.argv[1]
	motif = sys.argv[2]
	k = int(sys.argv[3])
	w = int(sys.argv[4])
	main(jobs, motif, k, w)
