#!/home/nfs/SOFTWARE/bin/Rscript
library("ggplot2")
require(ggthemes)
argv <- commandArgs(TRUE)

countsTable <- read.delim( argv[1], header=TRUE, stringsAsFactors=TRUE )
pdf(argv[2])
ggplot(countsTable, aes(Range,NO))+ geom_bar(alpha=0.2,stat = "identity")+ylab("Unigene Number")+theme_few()+xlab("Unigene Length")+ theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+scale_x_continuous(breaks=seq(0, 8000, 200))
dev.off()
tiff(argv[3])
ggplot(countsTable, aes(Range,NO))+ geom_bar(alpha=0.2,stat = "identity")+ylab("Unigene Number")+theme_few()+xlab("Unigene Length")+ theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+scale_x_continuous(breaks=seq(0, 8000, 200))
dev.off()
