# our working directory is /home/sereno/projects/cell_clustering
# we need to be in /home/vagar/projects/forSereno/ [files]
library(data.table)
PATH <- "../../../vagar/projects/forSereno"

pathcsv <- paste(PATH,"/df_cell_cluster_1.csv",sep = "") ## need sep = "" to stop it from inserting space

DT <- data.table(read.csv(file = pathcsv, header = TRUE, colClasses = c("NULL","character","NULL")))
# returns one line data table of cell types
counts <- DT[,.N,by = "cell_name"] 
write.csv(counts, file = "cluster_counts.csv", quote = FALSE) ## quote = false to take out r weirdness

pdf("cluster_plot.pdf",width = 20, height = 10)
par(mar = c(15, 4.1, 4.1, 2.1)) ## sets bottom margin to larger than default to hold data
barplot(table(DT$cell_name), las = 2) ## las=2 for vertical labels
mtext(text="Cell Cluster Counts", side=1, line=12) ## line argument puts the main text somewhere below the labels
dev.off()