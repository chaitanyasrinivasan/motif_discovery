import align #CYTHON ALIGNMENT
import time
import numpy as np
import sys
import argparse

#Run the cython alignment for width inference

def main(input):
	align.main(input)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type=str, help="Input set of sequences")

	options, args = parser.parse_known_args()

	if (len(sys.argv)==1):
	    parser.print_help(sys.stderr)
	    sys.exit(1)
	elif (options.input is None):
		parser.print_help(sys.stderr)
		sys.exit(1)
	else:
		main(options.input)