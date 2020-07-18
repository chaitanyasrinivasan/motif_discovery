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
scanSeqs() {
	while IFS= read -r line; 
	do 
		if ((${#line} < $WIDTH))
		then
		echo "Error: width ${WIDTH} is greater than the length of a sequence"
		exit 1
	fi; 
	done < ${INPUT}
}

if [[ $TYPE = "FASTA" ]]
then
	scanSeqs
fi

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
	echo "Splitting data..."
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
	echo "Creating job array..."
	ls -d "${INPUT::-4}_splits"/*.txt >> jobs.txt
	#SUBMIT JOBS
	echo "Submitting jobs..."
	NUM_JOBS=$(wc -l <jobs.txt)
	sed "s/REPLACE/${NUM_JOBS}/g" template.sb > jobs.sb
	#Clear old error logs
	if [ $(grep err logs/* | wc -l) -gt 0 ]
	then
		rm logs/*err*
	fi
	sbatch --wait jobs.sb ${WIDTH}
	#Check that all jobs completed
	for file in logs/*err*;
	do
		if [ $(wc -l <$file) -gt 0 ]
		then
			echo "A job did not complete. Error logs written to logs/"
			exit 1
		fi
	done
	rm logs/*err*
	echo "All jobs ended, logs written to logs/"
	#MERGE OUTPUTS AND RUN SEQUENTIAL GIBBS
	echo "Running merge..."
	total_lines=$(wc -l <${INPUT})
	python merge.py jobs.txt ${INPUT} ${total_lines} ${WIDTH}
	#Program compelete, clean-up
	rm jobs.txt jobs.sb "${INPUT::-4}_shuffled.txt"
	rm -r "${INPUT::-4}_splits"
}

############### DETERMINE MODE ##########################

startAnalysis() {
	GRAN=20 #empirically determined, see motif_discovery/performance
	#Automatically infer if parallel or sequential is faster
	if [ $SEQUENTIAL -eq $PARALLEL ]
	then
		if [ $(wc -l <${INPUT}) -lt $GRAN ]
		then
			python run_gibbs.py -i $INPUT -w $WIDTH -m conquer
			exit 0
		else
			parallelRun
			exit 0
		fi
	fi
	if (( SEQUENTIAL ))
	then
		python run_gibbs.py -i $INPUT -w $WIDTH -m conquer
		exit 0
	fi
	if (( PARALLEL ))
	then
		parallelRun
		exit 0
	fi
}

############# BED TO FASTA ####################################

bedToFasta() {
	#Check BED is correctly formatted using bedtools quick command
	echo "Merging bed coordinates..."
	sort -k 1,1 -k 2,2n ${INPUT} | bedtools merge -i stdin > merged.bed
	#Download hg38 fasta
	if [ ! -f hg38.fa ]
	then
		wget -nc ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
		gunzip hg38.fa.gz
	fi
	bedtools getfasta -fi hg38.fa -bed merged.bed > "${INPUT::-4}.fa"
	#Clean up
	rm merged.bed
	gzip hg38.fa
	sh preprocess.sh "${INPUT::-4}.fa"
	#Redirect input var to processed fasta
	INPUT="${INPUT::-4}.txt" 	
}

###################  MAIN ####################################

if [[ $TYPE = "FASTA" ]]
then
	startAnalysis
fi

if [[ $TYPE = "BED" ]]
then
	bedToFasta
	#Quality check
	scanSeqs
	#Run motif discovery
	startAnalysis
fi

if [[ $TYPE = "GENE" ]]
then
	#Map genes to regulatory coordinates in hg38
	sh gene_map.sh ${INPUT} #outputs ../data/introns_and_flanks/
	#Map bed to fasta
	INPUT="${INPUT::-4}.bed"
	bedToFasta
	#Quality check
	scanSeqs
	#Run motif discovery
	startAnalysis
fi