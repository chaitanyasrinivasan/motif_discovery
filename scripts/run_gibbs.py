import cygibbs #SEQUENTIAL CYTHON IMPLEMENTATION
import time
import numpy as np
import sys
import argparse

#Run the sequential implementation of gibbs sampler

def main(input, width):
	cygibbs.main(input, width)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type=str, help="Input set of sequences")
	parser.add_argument("-w", "--width", type=str, help="Integer of motif width")

	options, args = parser.parse_known_args()

	if (len(sys.argv)==1):
	    parser.print_help(sys.stderr)
	    sys.exit(1)
	elif (options.input is None):
		parser.print_help(sys.stderr)
		sys.exit(1)
	elif (options.width is None):
		parser.print_help(sys.stderr)
		sys.exit(1)
	elif (int(options.width) < 1):
		raise Exception("Width must be an integer greater than 0")
		sys.exit(1)
	else:
		main(options.input, int(options.width))