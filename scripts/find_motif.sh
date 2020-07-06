#!/bin/bash
helpFunction()
{
	echo -e "Usage: $0 -i [/path/to/fasta] -w [motif size] -t [BED/FASTA/GENES]\n"
	echo -e "Required arguments:"
	echo -e "\t-i, --input\tFile path to the processed fasta file"
	echo -e "\t-w, --width\tPositive integer of motif width"
	echo -e "\t-t --type\t Supported types: BED/FASTA/GENES"
	echo -e "Optional arguments:"
	echo -e "\t-s, --sequential\tRun sequentially"
	echo -e "\t-p, --parallel\tRun in parallel"
	echo -e "\t-v, --verbose"
	echo -e "\t-h, --help\n"
	echo "Example run call below:"
	echo ""
	echo "$0 -i myfasta.txt -w 10"
	exit 1 # Exit script after printing help
}
# check for args
if [ $# -eq 0 ]; then
    helpFunction
    exit 1
fi

WIDTH=-1
VERBOSE=0
PARALLEL=0
SEQUENTIAL=0

# read in command line arguments
while [ "$1" != "" ]; do
	case $1 in
		-i | --input )			shift
								INPUT=$1
								;;
		-w | --width )			shift
								WIDTH=$1
								;;
		-t | --type )			shift
								TYPE=$1
								;;
		-v | --verbose )		shift
								VERBOSE=1
								;;
		-s | --sequential )		shift
								SEQUENTIAL=1
								;;
		-p | --parallel )		shift
								PARALLEL=1
								;;
		-h | --help )           helpFunction
								exit 1
								;;
		*)
			helpFunction
			exit 1
	esac
	shift

done
############## INPUT COMPATABILITY CHECKS ################

# check if input file exists
if [ $(ls ${INPUT}| wc -l) -eq 0 ]
	then
	echo "Error: input fasta ${INPUT} does not exist"
	exit 1
fi
# check if more than 1 sequence
if [ $(wc -l <${INPUT}) -lt 2 ]
	then
	echo "Error: input fasta ${INPUT} needs at least 2 sequences"
	exit 1
fi
# check input width
if (($WIDTH < 1))
	then
	echo "Error: width must be an integer greater than 0"
	exit 1
fi
#check width is greater than sequence lengths
while IFS= read -r line; 
do 
	if ((${#line} < $WIDTH))
	then
	echo "Error: width ${WIDTH} is smaller than a sequence length"
	exit 1
fi; 
done < ${INPUT}
if [[ $TYPE != "BED" && $TYPE != "FASTA" && $TYPE != "GENES" ]]
then
	echo "Error: Provided file type not supported"
	helpFunction
	exit 1
fi
############## PARALLEL ALGORITHM ################

parallelRun()
{
	mkdir -p "${INPUT::-4}_splits"
	#SHUFFLE
	if (( VERBOSE )); then
		echo "Splitting data"
	fi
	#Randomize sequence order in case nearby sequences are similar
	shuf ${INPUT} > "${INPUT::-4}_shuffled.txt"
	#PARTITION
	lines_per_file=5
	split -d --lines=${lines_per_file} --additional-suffix=.txt "${INPUT::-4}_shuffled.txt" "${INPUT::-4}_splits"/split
	#If there is a split file with length 1, add it to the first file
	last_file=$(ls -1 "${INPUT::-4}_splits"/split* | tail -n 1)
	first_file=$(ls -1 "${INPUT::-4}_splits"/split* | head -n 1)
	if [ $(wc -l <$last_file) -eq 1 ]
	then
		cat $last_file >> $first_file 
		rm $last_file
	fi;
	#CREATE JOB ARRAY OF PARTITIONS
	if (( VERBOSE )); then
		echo "Creating job array"
	fi
	ls -d "${INPUT::-4}_splits"/*.txt >> jobs.txt
	#SUBMIT JOBS
	if (( VERBOSE )); then
		echo "Submitting jobs..."
	fi
	NUM_JOBS=$(wc -l <jobs.txt)
	sed "s/REPLACE/${NUM_JOBS}/g" template.sb > jobs.sb
	sbatch --wait jobs.sb ${WIDTH}
	if (( VERBOSE )); then
		echo "All jobs ended, logs written to logs/"
	fi
	#MERGE OUTPUTS AND RUN SEQUENTIAL GIBBS
	if (( VERBOSE )); then
		echo "Running merge..."
	fi
	total_lines=$(wc -l <${INPUT})
	python merge.py jobs.txt ${INPUT} ${total_lines} ${WIDTH}
	rm jobs.txt
	rm jobs.sb
	rm "${INPUT::-4}_shuffled.txt"
	rm -r "${INPUT::-4}_splits"
}

############### DETERMINE MODE ##########################
GRAN=20 #empirically determined
if [[ $TYPE = "FASTA" ]]
then
	#Automatically infer if parallel or sequential is faster
	if [ $SEQUENTIAL -eq $PARALLEL ]
	then
		if [ $(wc -l <${INPUT}) -lt $GRAN ]
		then
			python run_gibbs.py -i $INPUT -w $WIDTH
			exit 1
		else
			parallelRun
			exit 1
		fi
	fi
	if (( SEQUENTIAL ))
	then
		python run_gibbs.py -i $INPUT -w $WIDTH
		exit 1
	fi
	if (( PARALLEL ))
	then
		parallelRun
		exit 1
	fi
fi
if [[ $TYPE = "BED" ]]
then
	#Check BED is correctly formatted using bedtools quick command

	#Download hg38 fasta
	wget -nc ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	bedtools getfasta -fi hg38.fa.gz -bed ${INPUT} > "${INPUT::-4}.fa"
	sh preprocess.sh "${INPUT::-4}.fa"
	#Redirect input var to processed fasta
	INPUT="${INPUT::-4}.txt" 
	#Automatically infer if parallel or sequential is faster
	if [ $SEQUENTIAL -eq $PARALLEL ]
	then
		if [ $(wc -l <${INPUT}) -lt $GRAN ]
		then
			python run_gibbs.py -i $INPUT -w $WIDTH
			exit 1
		else
			parallelRun
			exit 1
		fi
	fi
	if (( SEQUENTIAL ))
	then
		python run_gibbs.py -i $INPUT -w $WIDTH
		exit 1
	fi
	if (( PARALLEL ))
	then
		parallelRun
		exit 1
	fi
fi
if [[ $TYPE = "GENE" ]]
then
	sh gene_map.sh ${INPUT}
fi