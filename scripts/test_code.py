import cygibbs #SEQUENTIAL CYTHON IMPLEMENTATION
import time
import numpy as np
import sys


iters=1

print("Testing sequential implementation...")

sequential_time = np.empty(iters)
for i in range(iters):
	start = time.time()
	cygibbs.main(sys.argv[1], sys.argv[2])
	end = time.time()
	sequential_time[i] = (end-start)
	print("Sequential time = {}".format(sequential_time[i]))
print("Sequential average time = {}".format(np.mean(sequential_time)))

'''
short_options = "hi:n:o:v"
	long_options = ["help", "input", "length", "out", "verbose"]
	argument_list = sys.argv[1:]
	try:
		arguments, values = getopt.getopt(argument_list, short_options, long_options)
		for arg,val in arguments:
			if arg in ("-h", "--help"):
				print("Gibbs Sampler will run motif discovery on an modified input fasta file.\n\n")
				print("Usage:\t ")

	except getopt.error as err:
		print(str(err))
		sys.exit(2)
'''