import sys
import random
import numpy as np
import concurrent.futures
import pandas as pd
import seqlogo

#INPUT: .fasta file of motifs

def count(seq):
    bases = {"a":0, "c":1, "g":2, "t":3}
    freq = {0:0, 1:0, 2:0, 3:0}
    for i in range(len(seq)):
        freq[bases[seq[i]]] += 1
    return freq

def parse(path):
    #Compute nucleotide frequencies in parallel
    with concurrent.futures.ProcessPoolExecutor() as executor:
        seqs = np.loadtxt(path, dtype="str")
        counter = 0
        b = 4 # number of bases
        freqs = np.empty(len(seqs)*4).reshape(len(seqs), b)
        for freq in executor.map(count, seqs):
            for j in range(b):
                freqs[counter][j] = freq[j]
            counter += 1
    val = (freqs.sum(axis=0))
    return seqs, val/val.sum()

def propensity(A, w, k, freq):
    b = 4
    P = np.empty(shape=(b, w))
    pseudo = np.sqrt(k)*0.25
    with concurrent.futures.ProcessPoolExecutor() as executor:
        cols = [A[:,j] for j in range(w)]
        counter = 0
        for counts in executor.map(count, cols):
            for i in range(b):
                P[i][counter] = (((counts[i] + pseudo)/(k + (pseudo*b))))/freq[i]
            counter += 1
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

def pdfGenerate(o, w, P, Pindex, tStar, nStar):
    productNumerator = 1
    for j in range(w):
        productNumerator *= P[Pindex[tStar[o+j]]][j]
    sumDenominator = 0
    for i in range(nStar-w):
        productDenominator = 1
        for j in range(w):
            productDenominator *= P[Pindex[tStar[i+j]]][j]
        sumDenominator += productDenominator
    return productNumerator/sumDenominator

def search(P, A, tStar, nStar, w, k, index, z, seqs, freq):
    counter = 0
    stop = 0
    while(stop < len(seqs)):
        #Calculate PDF
        Pindex = {"a":0, "c":1, "g":2, "t":3, "A":0, "C": 1, "G":2, "T":3}
        pdf = np.empty(shape=nStar-w, dtype=float)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            i = 0
            for val in executor.map(pdfGenerate, [o for o in range(nStar-w)], [w for i in range(nStar-w)], [P for i in range(nStar-w)], [Pindex for i in range(nStar-w)], [tStar for i in range(nStar-w)], [nStar for i in range(nStar-w)]):
                pdf[i] = val
                i += 1
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
        if np.allclose(P, Pnew):
            stop += 1
        else:
            stop = 0
        P = Pnew
        counter += 1

    #Obtain A by one last update
    pdf = np.empty(shape=nStar-w, dtype=float)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        i = 0
        for val in executor.map(pdfGenerate, [o for o in range(nStar-w)], [w for i in range(nStar-w)], [P for i in range(nStar-w)], [Pindex for i in range(nStar-w)], [tStar for i in range(nStar-w)], [nStar for i in range(nStar-w)]):
            pdf[i] = val
            i += 1
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
    return S, Anew

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
    S, Anew = search(P, A, tStar, nStar, w, k, index, z, seqs, freq)
    plotLogo(fasta, Anew, w)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2])
