import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
	sequential= np.loadtxt("sequential_metrics.txt", dtype=float)
	parallel= np.loadtxt("parallel_metrics.txt", dtype=float)
	size = np.arange(len(sequential), dtype=int)
	fig, ax1 = plt.subplots()
	ax1.plot(size, sequential, 'g-', label="Sequential Time (seconds)")
	ax1.plot(size, parallel, 'b-', label="Parallel Time (seconds)")
	fig.savefig("performance.png")


