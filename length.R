#!/home/nfs/SOFTWARE/bin/Rscript 
library("ggplot2")
require(ggthemes)
argv <- commandArgs(TRUE)

countsTable <- read.delim( argv[1], header=TRUE, stringsAsFactors=TRUE )
pdf(argv[2])
ggplot(countsTable, aes(Length))+ geom_bar(alpha=0.2)+ylab("Unigene Number")+theme_few()
dev.off()
tiff(argv[3])
ggplot(countsTable, aes(Length))+ geom_bar(alpha=0.2)+ylab("Unigene Number")+theme_few() 
dev.off()
