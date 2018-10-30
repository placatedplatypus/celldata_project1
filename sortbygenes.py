import sys
import csv
import re


def reversecomplement(seq):
    seq = seq.upper()
    seq = seq.replace('A', 't')
    seq = seq.replace('T', 'a')
    seq = seq.replace('C', 'g')
    seq = seq.replace('G', 'c')
    seq = seq.upper()
    seq = seq[::-1]
    return seq


# takes a fasta followed by a bed file followed by desired output file as inputs.
infasta = open(sys.argv[1], 'r')
inbed = sys.argv[2]
outfile = sys.argv[3]

genedict = {}


for line in infasta:
    if re.match('^>', line):
        gene = line.replace('>', '').rstrip()
        if gene not in genedict:    # check if it's a new gene
            genedict[gene] = []     # initializes an empty gene info list
    else:
        if genedict[gene]:          # if gene code is already started, adds next part of sequence
            genedict[gene][0].append(line.rstrip())
        else:                       # if gene code is empty, starts the gene code.
            genedict[gene] = [[line.rstrip()]]
        # appends line to current genecode reading if it's not a header line.
# genedict is now an entire dictionary of gene lists, currently only with their constituent codes as first entries.
infasta.close()

with open(inbed, 'rt') as tsvin:
    bedinfo = csv.reader(tsvin, delimiter='\t')
    for row in bedinfo:
        start = row[1]
        stop = row[2]
        gene = row[3]
        if len(genedict[gene]) == 1:   # checks if the gene has any entries other than its code
            chr = row[0]
            strand = row[5]
            genedict[gene].append(chr)    # adds one-off info to the gene if it's the first encounter
            genedict[gene].append(strand)
        genedict[gene].append(start)   # adds start and stop pairs for every entry in the gene
        genedict[gene].append(stop)
# genedict is now a dictionary of genes, with code, chr, strand, and then all start-stop pairs.
tsvin.close()

with open(outfile, 'w') as out:
    for gene in genedict:
        sequence = genedict[gene][0]
        for i in range(0, (len(genedict[gene][3:])-2), 2):  # checks all starts and stops for overlap
            stop = genedict[gene][i+4]
            start = genedict[gene][i+5]
            if stop == start:
                sequence[i//2] = sequence[i//2][:-1]        # clips overlap nucleotide
                genedict[gene][i+4] = str(int(stop) - 1)         # corrects read range for future algorithms
        if genedict[gene][2] == "+":    # if positive strand, forward-splices exons
            genedict[gene][0] = "".join(sequence)
        else:                           # if negative strand, reverse-splices exons and inverts sequence.
            sequence = "".join(list(reversed(sequence)))
            sequence = reversecomplement(sequence)
            genedict[gene][0] = "".join(sequence)
        geneinfo = '\t'.join(genedict[gene][1:])    # writes chr strand starts and stops as TSV
        out.write('>' + gene + '\t' + geneinfo + '\n')
        out.write(genedict[gene][0] + '\n')    # writes entire genetic code to next line
out.close()
