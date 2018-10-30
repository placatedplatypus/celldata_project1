library(plyr)

geneinfo <- read.table("mm10_longerUTR.gtf")

# split the gtf into by gene
corrected <- ddply(geneinfo, ~ V9, function(gene) {
	# only consider exons, as we're only modifying them
	exons <- subset(gene, V3!="five_prime_UTR" & V3!="three_prime_UTR")
	if (gene[1,7]=="+"){
		localmax <- max(exons[,5])
		gene[gene==localmax] <- (localmax + 1)
	}
	# adds one to the last exon of positive strands
	else {
		localmin <- min(exons[,4])
		gene[gene==localmin] <- (localmin - 1)
	}
	# adds one to the last exon of negative strands
	gene
})

write.table(corrected, file='mm10_longerUTR_c.gtf', quote=FALSE, sep='\t', col.names = NA)
