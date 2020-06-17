import sys
import random
import math
from libc.stdlib cimport malloc, free

#INPUT: .fasta file of motifs
#OUTPUT: dictionary of sequences, frequency of nucleotides
cpdef parse(str path) nogil:
	cdef dict seqs = dict()
	cdef int count
	cdef str line
	cdef int start
	cdef str data
	cdef str base
	cdef dict nucleotides = {"a":0, "c":0, "g":0, "t":0}
	cdef dict freq = {"a":0, "c":0, "g":0, "t":0}
	#Read in fasta file
	with open(path, "r") as f:
		for count, line in enumerate(f, start=0):
				data = str(line).strip("\n")
				seqs[count] = [data, len(data)]
				#Get background frequency
				for base in data:
					if base == "a" or base == "A": nucleotides["a"] += 1
					elif base == "c" or base == "C": nucleotides["c"] += 1
					elif base == "g" or base == "G": nucleotides["g"] += 1
					else:
						nucleotides["t"] += 1
	cdef int total = sum(nucleotides.values())
	freq["a"] = nucleotides["a"]/total
	freq["c"] = nucleotides["c"]/total
	freq["g"] = nucleotides["g"]/total
	freq["t"] = nucleotides["t"]/total
	return (seqs, freq)

cpdef propensity(A, w, k, z, freq, inputDict, index):
	cdef double *P = <double *> malloc(4 * w * sizeof(double))
	cdef int pseudo = 1
	cdef int aCount = 0
	cdef int cCount = 0
	cdef int gCount = 0
	cdef int tCount = 0
	cdef int search
	cdef double qA
	cdef double qC
	cdef double qG
	cdef double qT
	cdef str base
	try:
	#Row 0: A, Row1: C, Row 2: G, Row 3: T
		for j in range(w):
			aCount = 0
			cCount = 0
			gCount = 0
			tCount = 0
			for i in range(len(index)):
				search = A[int(j+i*w)]
				base = inputDict[int(index[i])][0][search]
				if base == "a" or base == "A": aCount += 1
				elif base == "c" or base == "C": cCount += 1
				elif base == "g" or base == "G": gCount += 1
				else:
					tCount += 1
			qA = (aCount + pseudo)/(k+(pseudo*4))
			qC = (cCount + pseudo)/(k+(pseudo*4))
			qG = (gCount + pseudo)/(k+(pseudo*4))
			qT = (tCount + pseudo)/(k+(pseudo*4))
			P[0+j*4] = qA/freq["a"]
			P[1+j*4] = qA/freq["c"]
			P[2+j*4] = qA/freq["g"]
			P[3+j*4] = qA/freq["a"]
		return [x for x in P[:4*w]]
	finally:
		free(P) 
		P = NULL

cpdef logOdds(P, w):
	cdef double *S = <double *> malloc(4 * w * sizeof(double))
	cdef int i
	cdef int j
	try:
		for i in range(4):
			for j in range(w): # A C G T
				S[j+i*w] = math.log(P[j+i*4], 2)
		return [x for x in S[:4*w]]
	finally:
		free(S)
		S = NULL

#INPUT: dictionary of nucleotides, motif length, frequency of nucleotides
#OUTPUT: P, t*, n*, k, index, z
cpdef init(dict inputDict, int w, dict freq):
	cdef int z = 0
	cdef str tStar = inputDict[z][0]
	cdef int nStar = inputDict[z][1]
	cdef int k = max(inputDict, key=int)
	cdef double *o = <double *> malloc(k * sizeof(double))
	cdef double *index = <double *> malloc((k-1) * sizeof(double))
	cdef double *A = <double *> malloc((k-1) * w * sizeof(double))
	cdef int j
	cdef int m
	cdef int pseudo 
	#Select random offsets and make a matrix of the coordinate indices
	try:
		for j in range(k):
			if (j != z):
				index[j-1] = j
				o[j] = random.randint(0, inputDict[j][1]-w-1)
				for m in range(w):
					A[m+(j-1)*w] = int(o[j]+m)
		P = propensity([x for x in A[:(k-1)*w]], w, k, z, freq, inputDict, [x for x in index[:(k-1)]])
		return (P, [x for x in A[:(k-1)*w]], tStar, nStar, k, [x for x in index[:(k-1)]], z)
	finally:
		free(o)
		o = NULL
		free(index)
		index = NULL
		free(A)
		A = NULL

cpdef search(P, A, tStar, nStar, w, k, index, z, inputDict, freq):
	cdef int counter = 0
	cdef int randIndex = random.randint(0, 3)
	cdef dict Pindex = {"a":0, "c":1, "g":2, "t":3, "A":0, "C": 1, "G":2, "T":3}
	cdef double *pdf = <double *> malloc((nStar-w)*sizeof(double))
	cdef double *cdf = <double *> malloc((nStar-w)*sizeof(double))
	cdef double *Anew = <double *> malloc(k*w*sizeof(double))
	cdef double productNumerator
	cdef double productDenominator
	cdef double sumDenominator
	cdef double sumCDF
	cdef double oRand
	cdef int oStar
	cdef int o
	cdef int j
	cdef int i
	cdef int l
	cdef int r
	cdef int y

	try:
		while(counter < 10000):
      #Calculate PDF
			for o in range(nStar-w):
				productNumerator = 1
				for j in range(w):
					productNumerator *= P[Pindex[tStar[o+j]]+j*4]
				sumDenominator = 0
				for i in range(nStar-w):
					productDenominator = 1
					for j in range(w):
						productDenominator *= P[Pindex[tStar[i+j]]+j*4]
					sumDenominator += productDenominator
				pdf[o] = productNumerator/sumDenominator
      #Calculate CDF
			for o in range(nStar-w):
				sumCDF = 0
				for l in range(o):
					sumCDF += pdf[l]
				cdf[o] = sumCDF
			#Choose oStar from CDF
			oRand = random.random()
			while(oRand > cdf[nStar-w-1]):
				oRand = random.random()
			oStar = -1
			for o in range(nStar-w):
				if oRand < cdf[o]:
					oStar = o
					break
				else:
					continue
			if oStar == -1: oStar = nStar-w-1
      #New special sequence
			r = random.randint(0, k-2)
      #Replace new special sequence with tStar coordinates
			for i in range(w):
				A[i+r*w] = oStar+i
      #Store and update parameters
			y = index[r]
			index[r] = z
			z = y
			tStar = inputDict[z][0]
			nStar = inputDict[z][1]
			P = [x for x in propensity(A, w, k, z, freq, inputDict, [x for x in index[:(k-1)]])]
			counter += 1
    #Obtain A by one last update
		for i in range(len(index)):
			for j in range(w):
				Anew[int(j+index[i]*w)] = A[j+i*w]
		for i in range(w):
			Anew[i+z*w] = oStar+i
		return logOdds(propensity([x for x in Anew[:(k*w)]], w, k, z, freq, inputDict, [x for x in index[:(k-1)]]), w)
	finally:
		free(pdf)
		pdf = NULL		
		free(cdf)
		cdf = NULL
		free(Anew)
		Anew = NULL

cpdef int main(str fasta):
	w = 19
	print("Parsing data")
	inputDict, freq = parse(fasta)
	print("Initalizing matrices")
	P, A, tStar, nStar, k, index, z = init(inputDict, w, freq)
	print("Searching for motif")
	S = search(P, A, tStar, nStar, w, k, index, z, inputDict, freq)
	for i in range(len(S)/4):
		print("A,"+str(i+1) + " : " + str(S[i]))
		print("C,"+str(i+1) + " : " + str(S[i+1]))
		print("G,"+str(i+1) + " : " + str(S[i+2]))
		print("T,"+str(i+1) + " : " + str(S[i+3]))
	print("Done")
	return 1

if __name__ == "__main__":
	main(sys.argv[1])
    
