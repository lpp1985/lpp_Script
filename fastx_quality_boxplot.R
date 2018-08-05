#!/home/nfs/SOFTWARE/bin/Rscript 
# functions
boxplotFromFastx <- function(df){
  xlim <- c(0, nrow(df))
  ylim <- c(min(0, min(df$lW)), max(df$rW)+1)
  par(mar=c(5,4,2,2)+0.1)
  plot(xlim, ylim, type="n", xlab="Read position", ylab="Average quality")
  grid1 <- floor(ylim[1]/2)*2
  grid2 <- ceiling(ylim[2]/2)*2
  abline(h=seq(grid1, grid2, by=2), lty=2, col=grey(0.3))
  for (i in 1:nrow(df)){
    plot_whisk(i, df$Q1[i], df$lW[i])
    plot_whisk(i, df$Q3[i], df$rW[i])
    plot_box(i, df$Q3[i], df$Q1[i])
    plot_med(i, df$med[i])
  }
}
plot_whisk <- function(x, y1, y2){
  segments(x0=x, x1=x, y0=y1, y1=y2, col="red")
  segments(x0=x-0.4, x1=x+0.4, y0=y2, y1=y2, col="red")
}
plot_box <- function(x, y1, y2){
  lines(x=c(x-0.4, x-0.4, x+0.4, x+0.4, x-0.4),
        y=c(y1, y2, y2, y1, y1),
        col="red")
}
plot_med <- function(x, y){
  segments(x0=x-0.4, x1=x+0.4, y0=y, y1=y, lwd=2)
}
# read arguments
args <- (commandArgs(TRUE))
data <- read.table(args[1], head=TRUE, as.is=TRUE)
if (length(args) > 1){
  outfile <- args[2]
} else {
  outfile <- paste(args[1], ".png", sep="")
}
png(file=outfile, w=1000, h=500)
boxplotFromFastx(data)
dev.off()

