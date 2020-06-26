import sys
import numpy as np

def parse(path):
  return np.loadtxt(path, dtype=str)

def pi(seqs, states):
	#ACGT
    k = len(seqs)
    counts = np.zeros(len(states))
    for d in range(k):
        counts[states[seqs[d]]] += 1/k
    return counts

def transition(seqs, states):
  k = len(seqs)
  a = len(states)
	A = np.empty(shape=(a, a)
	for d in range(k):
		seq = seqs[d]
		for i in range(len(seq)-1):
			A[states[seq[i]]][states[seq[i+1]]] += 1
	for i in range(a):
		for j in range(a):
			A[i][j] = A[i][j]/np.sum(A[i])
	return A

def emission(seqs, states):
    #ACGT
	counts = np.zeros(len(states))
	for d in range(len(seqs)):
		seq = seqs[d]
		for i in range(len(seq)):
			counts[states[seq[i]]] += 1
	return counts/np.sum(counts)

def forward(seq, states, initStates, A, emit):
	t = len(seq)
	n = len(states)
	#INITIALIZATION
	alpha = np.empty(shape=(t, n), dtype=float)
	for j in range(n):
		alpha[0][j] = initStates[states[seq[0]]]*emit[j]
	#DYNAMICALLY RECURSE
	for i in range(1,t):
		for j in range(n):
			val = 0
			for k in range(n):
				val += alpha[i-1][k]*A[k][j]*emit[k]
			alpha[i][j] = val
	prob = np.sum(alpha[t-1])
	return alpha, prob

def backward(seq, states, initStates, A, emit):
	t = len(seq)
	n = len(states)
	beta = np.empty(shape=(t, n), dtype=float)
	#INITIALIZATION
	for i in range(n):
		val = 0
		for j in range(n):
			val += A[i][j]*emit[j]
		beta[t-1][i] = val
    #RECURSION
	for i in range(t-2, -1, -1):
		for j in range(n):
			val = 0
			for k in range(n):
				val += A[j][k]*emit[k]*beta[i+1][k]
			beta[i][j] = val
	prob = 0
	for j in range(n):
		prob += initStates[j]*emit[j]*beta[1][j]
	return beta, prob

def posterior_decode(alpha, beta):
  t = len(alpha


 def BW(seqs):
 	for d in range(len(seqs)):
 		alpha, probA = forward(inputDict[d])
 		beta, probB = backward(inputDict[d])



def main(fasta):
    states = {"a": 0, "c": 1, "g": 2, "t":3}
    seqs = parse(fasta)
    initStates = pi(seqs, states)
    A = transition(seqs, states)
    emit = emission(seqs, states)
    alpha, probA = forward(inputDict[0], states, initStates, A, emit)
    beta, probB = backward(inputDict[0], states, initStates, A, emit)


if __name__ == "__main__":
    main(sys.argv[1])
