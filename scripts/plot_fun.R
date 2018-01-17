#!/usr/bin/env Rscript
# input is *.total.txt file from the regression

require(gplots)
options("scipen"=100)

plotreg <- function(fl, verbose=FALSE) {
    print(fl)
    x=read.table(fl, head=T)
    #sub = subset(x, Chr=="chr7" & position=="6065804")
    ## p-value = 1.76362e-15
    corr = round(cor(x$prechipfreq, x$POSTfreq, method="pearson"), 2)

    y <- densCols(x$prechipfreq, x$POSTfreq, colramp=colorRampPalette(c("black", "white")))
    x$dens <- col2rgb(y)[1,] + 1L

    ## Map densities to colors
    cols <-  colorRampPalette(c("#054A91", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
    x$col <- cols[x$dens]

    ## Plot it, reordering rows so that densest points are plotted on top
    png(paste(fl, ".prepost.png", sep=""), width=1200, height=1200, pointsize=40)
    plot(x$POSTfreq~x$prechipfreq, data=x[order(x$dens),], pch=20, col=col, cex=0.5, xlim=c(0,1), ylim=c(0,1),las=1, cex.axis=1.5, main=paste(fl, corr, sep="\t"), xlab="Pre-ATAC frequency", ylab="Post-ATAC frequency")
    #par(new=T)
    #plot(sub$prechipfreq, sub$POSTfreq, xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n' ,xlab="", ylab="", col="red", pch=16, cex=1.5)
    dev.off()
}

args <- commandArgs(trailingOnly = TRUE)

for (i in args) {
    plotreg(i)
}

