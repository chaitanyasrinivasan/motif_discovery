#!/bin/bash

#Download GENCODE v33 hg38 gene coordinates
if [ ! -f gencode.v33.annotation.gff3 ]
then
	wget -nc ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz
	gunzip gencode.v33.annotation.gff3.gz
fi

cd ../data
mkdir -p genes exons exons_merged introns genes_20KBflank

#GET GENE AND EXON COORDINATES
file=$1
echo "Mapping gene names to gene and exon coordinates..."

python ../scripts/get_gene_and_exon_coordinates.py $file

mv *genes.bed genes/
mv *exons.bed exons/
#MERGE, SUBTRACT EXONS FROM GENE, AND FLANK GENE 20KB
#TO CAPTURE REGULATORY ACTIVITY
echo "Merging exons..."

for file in exons/*_exons.bed;
do
	sort -k1,1 -k2,2n $file | bedtools merge -i stdin > "${file::-4}_merged.bed"
done
mv exons/*_merged.bed exons_merged/

echo "Subtracting merged exons from gene and getting gene flanks..."

i=0
for file in genes/*_genes.bed;
do
	genes[$i]=$file
	((i++))
done
i=0
for file in exons_merged/*_exons_merged.bed;
do
	exons[$i]=$file
	((i++))
done
arraylength=${#genes[@]}
for ((i=0; i<${arraylength}; i++));
do
	sort -k 1,1 -k 2,2n ${genes[i]} > "${genes[i]::-4}_sorted.bed"
	bedtools subtract -a "${genes[i]::-4}_sorted.bed" -b ${exons[i]}  > "${genes[i]::-9}introns.bed"
	bedtools flank -i "${genes[i]::-4}_sorted.bed" -g ../scripts/hg38.chrom.sizes.sorted -b 20000 | bedtools subtract -a stdin -b "${genes[i]::-4}_sorted.bed" > "${genes[i]::-4}_20KBflank.bed"
done
rm genes/*genes_sorted.bed
mv genes/*introns.bed introns/
mv genes/*_20KBflank.bed genes_20KBflank/

#ADD GENE FLANKS TO INTRONS

echo "Adding gene flanks to introns..."

i=0
for file in genes_20KBflank/*_20KBflank.bed;
do
	flanks[$i]=$file
	((i++))
done
i=0
for file in introns/*_introns.bed;
do
	introns[$i]=$file
	((i++))
done
arraylength=${#flanks[@]}
for ((i=0; i<${arraylength}; i++));
do
	cat ${flanks[i]} ${introns[i]} | sort -k 1,1 -k 2,2n | bedtools merge -i stdin > "${introns[i]::-4}_and_flanks.bed"
done
## FIX THIS BUG
mv introns/*_and_flanks.bed "${file::-4}.bed"

#CLEANUP AND COMPLETE
gzip ../scripts/gencode.v33.annotation.gff3
rm -r genes exons exons_merged introns genes_20KBflank
echo "Done."
exit 0

