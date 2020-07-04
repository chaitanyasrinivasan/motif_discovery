import sys
import gzip

if __name__ == '__main__':
	geneMarker = sys.argv[1]
	geneoutfile = geneMarker[:-4]+"_genes.bed"
	exonoutfile = geneMarker[:-4]+"_exons.bed"
    #Read in hg38 gene coordinates
	with gzip.open("../scripts/gencode.v33.annotation.gff3.gz", 'r') as f, open(geneMarker, 'r') as g, open(geneoutfile, 'w') as h, open(exonoutfile, 'w') as e:
		geneCoordinates = dict()
		exonCoordinates = dict()
		for linef in f:
			if (linef[0] != "#"):
				genef = linef.split("\t")
				if (genef[2] == "gene"):
					chrf = genef[0]
					startf = genef[3]
					stopf = genef[4]
					geneInfof = genef[8].split(";")
					geneNamef = geneInfof[3].split("=")[1]
					geneCoordinates[str(geneNamef)] = (str(chrf),str(startf),str(stopf))
				elif (genef[2] == "exon"):
					chrf = genef[0]
					startf = genef[3]
					stopf = genef[4]
					geneInfof = genef[8].split(";")
					geneNamef = geneInfof[5].split("=")[1]
					if (geneNamef in exonCoordinates): #build a list of exonic coordinates for each gene
						exonCoordinates[str(geneNamef)].append([str(chrf),str(startf),str(stopf)])
					else:
						exonCoordinates[str(geneNamef)] = [[str(chrf),str(startf),str(stopf)]]
		lines = g.readlines()
		for lineg in lines:
			geneNameg = lineg.strip("\n")
			if (geneNameg in geneCoordinates):
				info = geneCoordinates[str(geneNameg)]
				h.write(info[0]+"\t"+info[1]+"\t"+info[2]+"\t"+geneNameg+"\n")
			if (geneNameg in exonCoordinates):
				infoList = exonCoordinates[str(geneNameg)]
				for i in range(len(infoList)):
					e.write(infoList[i][0]+"\t"+infoList[i][1]+"\t"+infoList[i][2]+"\t"+geneNameg+"\n")
f.close()
g.close()
h.close()
e.close()
