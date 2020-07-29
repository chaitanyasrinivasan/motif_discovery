import numpy as np
import matplotlib.pyplot as plt

'''
Chaitanya Srinivasan

This script plots the sequential and parallel performance metrics for
find_motif.sh
'''

if __name__ == "__main__":
	sequential= np.loadtxt("sequential_metrics.txt", dtype=float)
	parallel= np.loadtxt("parallel_metrics.txt", dtype=float)
	size = np.arange(5, len(sequential)+5, dtype=int)
	fig, ax1 = plt.subplots()
	ax1.plot(size, sequential, 'g-', label="Sequential")
	ax1.plot(size, parallel, 'b-', label="Parallel")
	ax1.set_xlabel("Number of sequences")
	ax1.set_ylabel("Time (seconds)")
	fig.suptitle("Gibb's Sampler Runtime")
	ax1.legend()
	fig.savefig("performance.png")


