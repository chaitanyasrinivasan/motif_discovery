import gibbs
import parallel_gibbs
import time
import numpy as np
import sys

iters=3

print("Testing sequential implementation...")
sequential_time = np.empty(5)
for i in range(iters):
	start = time.time()
	gibbs.main(sys.argv[1], sys.argv[2])
	end = time.time()

	sequential_time[i] = (end-start)
	print("Sequential time = {}".format(sequential_time[i]))
print("Sequential average time = {}".format(np.mean(sequential_time)))
print("\n")
print("Testing parallel implementation...")
parallel_time = np.empty(5)
for i in range(iters):
	start = time.time()
	parallel_gibbs.main(sys.argv[1], sys.argv[2])
	end = time.time()

	parallel_time[i] = end-start
	print("Parallel time = {}".format(parallel_time[i]))
print("Parallel average time = {}".format(np.mean(parallel_time)))
