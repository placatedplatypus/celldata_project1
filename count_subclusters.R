# takes input of a csv of cell clustering data in a folder named "data"
# gives output of a TSV txt file of subcluster frequency and a histogram of such data.
library(plyr)

pathcsv <- "data/df_cell_cluster_1.csv"

df <- data.frame(read.csv(file = pathcsv, header = TRUE, colClasses = c("NULL","character","character")))
# returns two line data frame of cell types and subclusters.
counts <- ddply(df, .(df$cell_name, df$sub_Cluster), nrow)
names(counts) <- c("cell_name", "sub_Cluster", "count") ##makes the counts prettier.
write.table(counts, file = "subcluster_counts.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# this gives us a table of cluster, subcluster, and frequency in a TSV table.

df$overcluster <- paste(df$cell_name,df$sub_Cluster, sep = "_") ## sep="_" inserts an underscore to the column.
# makes a new column called "overcluster" with the full cluster and sub-cluster
pdf("subcluster_plot.pdf",width = 200, height = 12)
par(mar = c(18, 4.1, 4.1, 2.1)) ## sets bottom margin to larger than default to hold data
barplot(table(df$overcluster), las = 2) ## las=2 for vertical labels
mtext(text="Sub Cluster Counts", side=1, line=15) ## line argument puts the main text somewhere below the labels
dev.off()