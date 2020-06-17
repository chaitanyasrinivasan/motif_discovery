import sys

def parse(path):
    seqs = dict()
    nucleotides = {"a":0, "c":0, "g":0, "t":0}
    freq = {"a":0, "c":0, "g":0, "t":0}
    with open(path, "r") as f:
        for count, line in enumerate(f, start=0):
            data = str(line).strip("\n")
            seqs[count] = data
    return seqs

def pi(inputDict, states):
	#ACGT
    k = len(inputDict)
    counts = {0:0, 1:0, 2:0, 3:0}
    for d in range(len(inputDict)):
        counts[states[inputDict[d][0]]] += 1/k
    return counts

def transition(inputDict, states):
	A = [[0 for i in range(4)] for j in range(4)]
	for d in range(len(inputDict)):
		seq = inputDict[d]
		for i in range(len(seq)-1):
			A[states[seq[i]]][states[seq[i+1]]] += 1
	for i in range(4):
		observed = sum(A[i])
		for j in range(4):
			A[i][j] = A[i][j]/observed
	return A

def emission(inputDict, states):
    #ACGT
	counts = [0, 0, 0, 0]
	for d in range(len(inputDict)):
		seq = inputDict[d]
		for i in range(len(seq)):
			counts[states[seq[i]]] += 1
	total = sum(counts)
	for i in range(len(counts)):
		counts[i] = counts[i]/total
	return counts


def forward(seq, states, initStates, A, emit):
	t = len(seq)
	n = 4
	#INITIALIZATION
	alpha = [[0 for j in range(n)] for i in range(t)]
	for j in range(n):
		alpha[0][j] = initStates[states[seq[0]]]*emit[j]
	#DYNAMICALLY RECURSE
	for i in range(1,t):
		for j in range(n):
			val = 0
			for k in range(n):
				val += alpha[i-1][k]*A[k][j]*emit[k]
			alpha[i][j] = val
	prob = sum(alpha[t-1])
	return alpha, prob

def backward(seq, states, initStates, A, emit):
	t = len(seq)
	n = 4
	beta = [[0 for i in range(n)] for j in range(t)]
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

 def BW(inputDict):
 	for d in range(len(inputDict)):
 		alpha, probA = forward(inputDict[d])
 		beta, probB = backward(inputDict[d])



def main(fasta):
    states = {"a": 0, "c": 1, "g": 2, "t":3, "A":0, "C": 1, "G": 2, "T": 3}
    inputDict = parse(fasta)
    initStates = pi(inputDict, states)
    A = transition(inputDict, states)
    emit = emission(inputDict, states)
    alpha, probA = forward(inputDict[0], states, initStates, A, emit)
    beta, probB = backward(inputDict[0], states, initStates, A, emit)


if __name__ == "__main__":
    main(sys.argv[1])