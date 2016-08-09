#!/usr/bin/Rscript

arg <- commandArgs(trailingOnly=TRUE)

rtt <- function(filepath)
{

require(TeachingDemos, quietly=TRUE, warn.conflicts=FALSE)
require(gplots, quietly=TRUE, warn.conflicts=FALSE)

cat("  reading log-file...\n")

tdat <- read.table(filepath, header=T)
headernames <- names(tdat)

if (length(grep("birth.rate", headernames)) == 0) {
	stop ("rates not found. please check if the input file is correct.")
}

nrates <- length(tdat) - 1
y.val <- NULL
upper <- NULL
lower <- NULL

for (n in 2:(nrates+1)) {
	y.val <- append(y.val,mean(tdat[,n]), after=length(y.val))
	upper <- append(upper, emp.hpd(tdat[,n], conf=0.95)[2], after=length(upper))
	lower <- append(lower, emp.hpd(tdat[,n], conf=0.95)[1], after=length(lower))
}

x.val <- seq(from=1, to=nrates)
x.val <- (x.val-0.5)*-1
dat <- cbind(x.val, y.val, upper, lower)
colnames(dat) <- c("Time", "Mean rate", "95%HPD upper", "95%HPD lower")

output <- gsub ("\\.log$", ".pdf", filepath)
pdf (file=output)
	plotCI(x=dat[,1], y=dat[,2], ui=dat[,3], li=dat[,4], err="y", col=NULL, barcol="black", lwd=1, sfrac=0.01, gap=0, xlim=c((nrates*-1),0), ylim=c(floor(min(dat[,4]*10))/10, ceiling(max(dat[,3]*10))/10), xlab="Time [myr]", ylab="Net diversification rate [r]", main="Rates-through-time plot", cex.main=0.9)
	par(new=TRUE)
	plot(spline(dat[,1], dat[,2]), type="l", xlim=c((nrates*-1),0), ylim=c(floor(min(dat[,4]*10))/10,ceiling(max(dat[,3]*10))/10), axes=FALSE, xlab="", ylab="")
dev.off()

output2 <- gsub ("\\.log$", ".txt", filepath)
write.table(round(dat, digits=4), file=output2, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
cat("  results are saved in: ", sprintf("%s", output), " and ", sprintf("%s", output2), "\n\n")
}

rtt(as.character(arg[1]))