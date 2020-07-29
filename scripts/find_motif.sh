#!/bin/bash

##### Chaitanya Srinivasan #####
##### Parallel de novo regulatory motif discovery tool #####

helpFunction()
{
	echo -e "Usage: $0 -i [/path/to/data] -w [motif size] -t [BED/FASTA/GENES]\n"
	echo -e "Required arguments:"
	echo -e "\t-i, --input\tFile path to the sequence, genomic coordinates, or genes list data."
	echo -e "Optional arguments:"
	echo -e "\t-w, --width\tPositive integer of motif width"
	echo -e "\t-s, --sequential\tRun sequentially"
	echo -e "\t-p, --parallel\tRun in parallel"
	echo -e "\t-h, --help\n"
	echo "Example run calls below:"
	echo ""
	echo "$0 -i myfasta.fa -w 10 -p"
	echo "$0 -i mybed.bed -s"
	echo "$0 -i mygenes.txt"
	exit 1 # Exit script after printing help
}
# check for args
if [ $# -eq 0 ]; then
    helpFunction
    exit 1
fi

# default glob vars
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
	echo "Error: input ${INPUT} does not exist"
	exit 1
fi
# check type is compatabile
if [[ ${INPUT: -4} != ".bed" && ${INPUT: -3} != ".fa" && ${INPUT: -4} != ".txt" ]]
then
	echo "Error: The file must have extension .fa, .bed, or .txt"
	helpFunction
	exit 1
fi

scanSeqs() {
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
	# check width is not greater than a sequence length
	while IFS= read -r line; 
	do 
		if ((${#line} < $WIDTH))
		then
		echo "Error: width ${WIDTH} is greater than the length of a sequence"
		exit 1
	fi; 
	done < ${INPUT}
}

#Remove sequences from FASTA that are not compatabile
preProcessing() {
	echo "Preprocessing fasta..."
	awk 'NR % 2 == 0 {print}' $INPUT > "${INPUT::-3}_seqs.txt"
	grep -vE "(X)" "${INPUT::-3}_seqs.txt" | grep -vE "(N)" | grep -vE "(n)" | tr '[:upper:]' '[:lower:]' > "${INPUT::-3}_filtered.txt"
	rm "${INPUT::-3}_seqs.txt"
	mv "${INPUT::-3}_filtered.txt" "${INPUT::-3}_seqs.txt"
	INPUT="${INPUT::-3}_seqs.txt"
}

############## WIDTH INFERENCE ######################

alignSeqs() {
	# run width inference if width not provided
	if [ -z "${WIDTH}" ];
	then
		echo "Inferring width from sequence alignment..."
		WIDTH=`python run_align.py -i ${INPUT}`
		echo "Setting motif width as ${WIDTH}..."
	fi	
}

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
	# clear old error logs
	if [ $(grep err logs/* | wc -l) -gt 0 ]
	then
		rm logs/*err*
	fi
	sbatch --wait jobs.sb ${WIDTH}
	# check that all jobs completed
	for file in logs/*err*;
	do
		if [ $(wc -l <$file) -gt 0 ]
		then
			echo "Error : Job ${file::-8} did not complete. Error logs written to ${file}"
			exit 1
		fi
	done
	rm logs/*err*
	echo "All jobs ended, logs written to logs/"
	#MERGE OUTPUTS AND RUN SEQUENTIAL GIBBS
	echo "Running merge..."
	total_lines=$(wc -l <${INPUT})
	python merge.py -j jobs.txt -f ${INPUT} -k ${total_lines} -w ${WIDTH}
	# program compelete, clean-up
	rm jobs.txt jobs.sb "${INPUT::-4}_shuffled.txt"
	rm -r "${INPUT::-4}_splits"
}

############### DETERMINE MODE ##########################

startAnalysis() {
	GRAN=20 #empirically determined, see motif_discovery/performance
	# automatically infer if parallel or sequential is faster
	if [ $SEQUENTIAL -eq $PARALLEL ]
	then
		if [ $(wc -l <${INPUT}) -lt $GRAN ]
		then
			python run_gibbs.py -i $INPUT -w $WIDTH -m conquer
			exit 0
		else
			# run parallel if slurm can be used
			if [ -x "$(command -v sbatch)" ]
			then
				parallelRun
				exit 0
			else
				python run_gibbs.py -i $INPUT -w $WIDTH -m conquer
				exit 0	
			fi
		fi
	fi
	# sequential run
	if (( SEQUENTIAL ))
	then
		python run_gibbs.py -i $INPUT -w $WIDTH -m conquer
		exit 0
	fi
	# parallel run
	if (( PARALLEL ))
	then
		if ! [ -x "$(command -v sbatch)" ] 
		then
			echo "Error : parallel job submission requires Slurm, use -s instead of -p"
			exit 1
		else
			parallelRun
			exit 0
		fi
	fi
}

############# BED TO FASTA ####################################

bedToFasta() {
	if [ -x "$(command -v bedtools)" ] 
	then
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

		INPUT="${INPUT::-4}.fa"
	else
		echo "Error: bedtools is not installed or is not executable from your path."
		exit 0
	fi
}

################### GENES TO BED #############################

genesToBed() {
	if [ -x "$(command -v bedtools)" ] 
	then
		#Download GENCODE v33 and hg38 chrom sizes
		if [ ! -f gencode.v33.annotation.gff3 ]
		then
			wget -nc ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz
			gunzip gencode.v33.annotation.gff3.gz
		fi
		wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

		#GET GENE AND EXON COORDINATES
		echo "Mapping gene names to gene and exon coordinates..."
		# returns gene and exon coordinates (2 files) for each gene in the list
		python ../scripts/get_gene_and_exon_coordinates.py $INPUT
		GENEBED="${INPUT::-4}_genes.bed"
		EXONBED="${INPUT::-4}_exons.bed"

		#MERGE, SUBTRACT EXONS FROM GENE, AND FLANK GENE 20KB TO CAPTURE REGULATORY ACTIVITY
		echo "Merging exons..."
		sort -k1,1 -k2,2n $EXONBED | bedtools merge -i stdin > "${EXONBED::-4}_merged.bed"
		EXONMERGEDBED="${EXONBED::-4}_merged.bed"
		echo "Subtracting merged exons from gene and getting gene flanks..."
		sort -k 1,1 -k 2,2n $GENEBED > "${GENEBED::-4}_sorted.bed"
		GENESORTEDBED="${GENEBED::-4}_sorted.bed"
		bedtools subtract -a $GENESORTEDBED -b $EXONMERGEDBED > "${INPUT::-4}_introns.bed"
		INTRONBED="${INPUT::-4}_introns.bed"
		bedtools flank -i $GENESORTEDBED -g hg38.chrom.sizes -b 20000 | bedtools subtract -a stdin -b $GENESORTEDBED > "${GENEBED::-4}_20KBflank.bed"
		GENEFLANKBED="${GENEBED::-4}_20KBflank.bed"

		#ADD GENE FLANKS TO INTRONS
		cat $GENEFLANKBED $INTRONBED | sort -k 1,1 -k 2,2n | bedtools merge -i stdin > "${INTRONBED::-4}_and_intergenics.bed"
		
		mv "${INPUT::-4}_introns_and_intergenics.bed" "${INPUT::-4}.bed"
		# set glob var to new format
		INPUT="${INPUT::-4}.bed"

		#CLEANUP AND COMPLETE
		gzip ../scripts/gencode.v33.annotation.gff3
		rm ${GENEFLANKBED} ${INTRONBED} ${GENESORTEDBED} ${EXONMERGEDBED} ${EXONBED} ${GENEBED}
		echo "Done mapping genes to BED."
	else
		echo "Error: bedtools is not installed or is not executable from your path."
		exit 0
	fi
}

###################  MAIN ####################################
#FASTA
if [[ ${INPUT: -3} = ".fa" ]]
then
	# preprocess sequences
	preProcessing
	# infer width if not provided
	alignSeqs
	# quality check
	scanSeqs
	# run motif discovery
	startAnalysis
fi
#BED
if [ ${INPUT: -4} = ".bed" ]
then
	# convert BED to FASTA
	bedToFasta
	# same as FASTA from here
	preProcessing 
	alignSeqs
	scanSeqs
	startAnalysis
fi
#GENES
if [ ${INPUT: -4} = ".txt" ]
then
	#Map genes to regulatory coordinates in hg38
	genesToBed
	# same as for BED from here
	bedToFasta
	preProcessing
	alignSeqs
	scanSeqs
	startAnalysis
fi