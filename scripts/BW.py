import sys
import numpy as np

def parse(path):
	#Load and convert fasta seqs into int seqs
	seqs = np.loadtxt(path, dtype=str)
	k = len(seqs)
	maxSize = len(max(seqs, key=len))
	intSeqs = np.array([[4 for j in range(maxSize)] for i in range(k)])
	seqLengths = np.empty(shape=k, dtype=int)
	bases = {"a":0, "c":1, "g":2, "t":3}
	for i in range(len(seqs)):
		seqLengths[i] = len(seqs[i])
		for j in range(len(seqs[i])):
			intSeqs[i][j] = bases[seqs[i][j]]
	return intSeqs, seqLengths

#L+2=10+2=12 (guess) match states
#11 insertion states
#10 deletion states
# number of states = 12+11+10=33
# double weight : 2/(33+12)
# number of observations = 4
def init(L):
	#Try initializing uniform
	m = L + 2
	i = L + 1
	d = L
	s = m+i+d#number of states
	b = 4 #number of possible observations(nucleotides)
	#index 0 to m-1 rows/cols are match
	#index m to m+i rows/cols are insert
	#index m+i+1 to m+i+d
	#matches transitions twice as likely as indel
	A = np.empty(shape=(s, s), dtype=float)
	E = np.array([[0.25 for _ in range(b)] for _ in range(s)])
	pi = np.zeros(shape=s, dtype=float)
	pi[0] = 1 #start state has certain probability
	for i in range(len(A)):
		for j in range(len(A[i])):
			if (i < m and j < m):
				A[i][j] = 2/45
			else:
				A[i][j] = 1/45
	return A, E, pi


def forward(seq, sLen, A, E, pi):
	t = sLen
	n = len(pi)
	#INITIALIZATION
	alpha = np.empty(shape=(t, n), dtype=float)
	for j in range(n):
		alpha[0][j] = pi[j]*E[j][seq[0]]
	#DYNAMICALLY RECURSE
	for i in range(1,t):
		for j in range(n):
			val = 0
			for k in range(n):
				val += alpha[i-1][k]*A[k][j]*E[j][seq[i]]
			alpha[i][j] = val
	prob = np.sum(alpha[t-1])
	return alpha, prob

def backward(seq, sLen, A, E, pi):
	t = sLen
	n = len(pi)
	beta = np.empty(shape=(t+1, n), dtype=float)
	#INITIALIZATION
	for i in range(n):
		val = 0
		for j in range(n):
			val += A[i][j]*E[j][seq[t-1]]
		beta[t][i] = val
	#RECURSION
	for i in range(t-1, -1, -1):
		for j in range(n):
			val = 0
			for k in range(n):
				val += A[j][k]*E[k][seq[i]]*beta[i+1][k]
			beta[i][j] = val
	prob = 0
	for j in range(n):
		prob += pi[j]*E[j][seq[0]]*beta[1][j]
	return beta, prob


def posterior_decode(alpha, beta, t):
	max, argmax = 0, 0
	for i in range(len(alpha[0])):
		if (alpha[t][i]*beta[t+1][i] > max):
			max = alpha[t][i]*beta[t+1][i]
			argmax = i
	return argmax


def BW(seqs, seqLengths, A, E, pi):
	b = len(pi)
	#Posterior decoding to label data
	labels = np.empty(shape=(len(seqs), len(seqs[0])))
	
	#Iterate through dataset
	for d in range(len(seqs)):
		t = seqLengths[d]
		alpha, prob = forward(seqs[d], seqLengths[d], A, E, pi) 
		beta, _ = backward(seqs[d], seqLengths[d], A, E, pi)
		Anew = np.zeros(shape=(t, t), dtype=float)
		Enew = np.zeros(shape=(t, b), dtype=float)
		pinew = np.zeros(shape=b, dtype=float)
		#Transitions for each seq
		for i in range(len(Anew)):
			for j in range(len(Anew[i])):
				for k in range(t-1):
					Anew[i][j] += alpha[k][i]*A[i][j]*E[j][seqs[d][k+1]]*beta[k+1][j]
				Anew[i][j] = Anew[i][j]/prob
		#Emissions for each seq
		for i in range(len(Enew)):
			for j in range(len(Enew[i])):
				for k in range(t):
					if (seqs[d][k] == j):
						Enew[i][j] += alpha[k][i]*beta[k+1][i]
				Enew[i][j] = Enew[i][j]/prob
		#Initial probability


	
def main(fasta):
	states = {"a": 0, "c": 1, "g": 2, "t":3}
	seqs, seqLengths = parse(fasta)
	L = 10
	A, E, pi = init(L) #guess on size of TF binding site
	BW(seqs, seqLengths, A, E, pi)


if __name__ == "__main__":
	main(sys.argv[1])
