#!/home/nfs/SOFTWARE/bin/Rscript
# functions
colors <- c("blue", "red", "green", "yellow", "pink")
nucDistrFromFastx <- function(df){
  xlim <- c(0, nrow(df)+1)
  #ylim <- c(0, max(df$count))
  ylim<-c(0,100)
  par(mar=c(5,4,2,2)+0.1, xaxs="i", yaxs="i")
  plot(xlim, ylim, type="n", xlab="Read position", ylab="The percentage of nucleotides components(%)")
  for (i in 1:nrow(df)){
    sum=df$A_Count[i]+df$C_Count[i]+df$G_Count[i]+df$T_Count[i]+df$N_Count[i]
    ys <- c(0,
            (df$A_Count[i])*100/sum,
            (df$A_Count[i]+df$C_Count[i])*100/sum,
            (df$A_Count[i]+df$C_Count[i]+df$G_Count[i])*100/sum,
            (df$A_Count[i]+df$C_Count[i]+df$G_Count[i]+df$T_Count[i])*100/sum,
            (df$A_Count[i]+df$C_Count[i]+df$G_Count[i]+df$T_Count[i]+df$N_Count[i])*100/sum
            )
    for (j in 1:5){
      plot_rect(i, ys[j], ys[j+1], colors[j])
    }
    legend("topright", legend=c("A", "C", "G", "T", "N"), fill=colors)
  }
}
plot_rect <- function(x, y1, y2, col){
  polygon(x=c(x-0.4, x-0.4, x+0.4, x+0.4, x-0.4),
        y=c(y1, y2, y2, y1, y1),
        col=col)
}
# read arguments
args <- (commandArgs(TRUE))
data <- read.table(args[1], head=TRUE, as.is=TRUE)
#data <- read.table("BBb.qualstats", head=TRUE, as.is=TRUE)
if (length(args) > 1){
  outfile <- args[2]
} else {
  outfile <- paste(args[1], ".png", sep="")
}
png(file=outfile, w=1000, h=500)
nucDistrFromFastx(data)
dev.off()
