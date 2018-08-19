import sys
import csv
from collections import defaultdict

## input is a gtf file
## make a dictionary of lists where keys are gene names and values are start-stop pairs

def genecoverage(row):
	if row[2] == "exon":
		gene = row[8]
		start = row[3]
		stop = row[4]

		if gene in genedict:
			oldstart = genedict[gene][1]
			oldstop = genedict[gene][2]

			if start < oldstart:
				del genedict[gene][1]
				genedict[gene].insert(1,start) #adds the new start location
				# print(genedict[gene]) #DEBUG
			if stop > oldstop:
				del genedict[gene][2]
				genedict[gene].insert(2, stop) #adds the new stop location
				# print(genedict[gene]) #DEBUG

		else:
			chrom = row[0]
			strand = row[6]
			geneinfo = [chrom, start, stop, strand]
			genedict[gene] = geneinfo	#will attach the start and stop values the first time a gene is encountered
			# print(genedict[gene]) #DEBUG
	# if the row isn't an exon it just does nothing lmao
	print(genedict) #DEBUG
	return genedict


infile = sys.argv[1]
outpath = infile.replace(".gtf","_parsed.gtf")
outfile = open(outpath, "w+")

genedict = {}

## field 0 is chr, field 2 read type, field 3 and 4 range of read, field 6 strand, field 8 gene name.
with open(infile) as tsvfile: 
	reader = csv.reader(tsvfile, delimiter ='\t')
	for row in reader:
		genecoverage(row) #after running this we have a full dictionary of genes with their biggest starts and stops


for gene in genedict: # tabs after the first three for that good good tsv formatting
	chrom = str(genedict[gene][0]) + "\t"
	start = str(genedict[gene][1]) + "\t"
	stop = str(genedict[gene][2]) + "\t"
	strand = str(genedict[gene][3]) + "\n"
	line = gene + chrom + start + stop + strand
	outfile.write(line)
outfile.close()

#this works
