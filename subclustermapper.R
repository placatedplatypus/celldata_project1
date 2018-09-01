library(EnrichedHeatmap)
library(circlize)
library(ggplot2)

# input is a single cluster bedgraph file passed to the R script
args = commandArgs(trailingOnly=TRUE)
raw1 = gsub("bigwig/","",args[1])
raw = gsub(".bedGraph","",raw1)  ## vectorize this later
# raw is now just the name of the cluster


filenames <- list.files(path = "bigwig/", pattern = raw) # returns a list of all subcluster files
filenames <- filenames[filenames != raw1] # removes overcluster bam from filenames list
# filenames is now a list of cluster.subcluster.bedGraph files.

n <- length(filenames) 	 # number of subclusters, used later.
c = as.character(1:19)
c <- append(c,c("X","Y"))
c <- sapply(c,function(x) paste0("chr",x))
names(c) <- c 										# why this was necessary we may never know
# c is now a character vector of approved chromosomes in df, we'll use this with the apply below us

df2 <- read.table("mousegenome.bed") # make sure you have this file lol
# can use IRanges(df2[[3]],df2[[4]]) to get TSS and TES. As it stands it only grabs TES
ref = GRanges(seqnames = df2[[2]], ranges = df2[[4]], strand = df2[[5]])
# this gets our read index
	    
plotname <- paste0(raw, ".png")
png(plotname)

range <- seq(-5000, 4950, by=50)
xrange <- range(range)
yrange <- range(seq(0,1, by=0.2)) 
# set up the plot 
plot(xrange, yrange, type="n", xlab="Distance From TES",
  	ylab="Read Density") 
colors <- rainbow(n) 

title(raw)
legend(xrange[1], yrange[2], 1:n, cex=0.8, col=colors,
  	pch=plotchar, lty=linetype, title="Subcluster")

overlist <- lapply(filenames, function(x) {
	a = read.table(paste0("bigwig/",x))			# paste function gives correct path
	subraw = gsub(".bedGraph","",x)  			# var is now "cluster.subcluster"
	subraw = gsub(paste0(raw,"."),"", subraw)	# var is now "subcluster"
	# There must be a better way to do this

	df <- a[a$V1 %in% c,]						# sanitizes chromosome names and input

	gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]),
			     coverage = log10(df[[4]]+1))	# converts our frame to a granges object
	
	mat = normalizeToMatrix(gr, ref, value_column = "coverage", extend = 5000, 
	mean_mode = "w0", w = 50)					# makes a matrix of binned reads

	readlist <- as.numeric(lapply(1:200, function(x) sum(mat[,x])))
	normalizedreads <- (readlist / max(readlist))
})
# overlist is a list of listed reads as numberic vectors

for (i in 1:n) {
	lines(range, overlist[[i]], type="l", lwd=1.5, col=colors[i])
}

dev.off()
