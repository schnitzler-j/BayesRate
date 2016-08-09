#!/usr/bin/Rscript

arg <- commandArgs(trailingOnly=TRUE)

rate.shift <- function(filepath)
{

require(hdrcde, quietly=TRUE, warn.conflicts=FALSE)
require(TeachingDemos, quietly=TRUE, warn.conflicts=FALSE)

cat("  reading log-file...\n")


tdat <- read.table(filepath, header=T)
headernames <- names(tdat)

if (length(grep("shift", headernames)) == 0) {
	stop ("rate-shifts not found. please check if the input file is correct.")
}

shifts <- length(tdat) - 1
if (round(length(tdat[,1])/500) < 100){
	histbreaks <- 100
	} else {
	histbreaks <- round(length(tdat[,1])/500)	
	}
	
tdat2 <- tdat
tdat2[,2:length(tdat2)] <- tdat2[,2:length(tdat2)]*-1

maxX <- 0
minX <- floor(min(tdat2[,2]))

output <- gsub ("\\.log$", ".pdf", filepath)
pdf (file=output)

par(mfcol=c(shifts,1))

for (n in 2:(shifts+1)) {
	sp <- round(hdr(tdat2[,n], all.modes=F)$mode, digits=3)
	cred_int <- emp.hpd (tdat[,n], conf=0.95)
	ahist <- hist(tdat2[,n], breaks=histbreaks, plot=F)
	if (floor(minX/5) == 0) {
		ticks <- 1
		} else {
			ticks <- floor(minX/5)*-1
		}
	hist(tdat2[,n], xlab = paste("Time"), sub = paste(sprintf("MAP: %s", sp), paste(", 95% HPD:", sprintf("%.3f-%.3f", cred_int[1], cred_int[2]))), main = "Time of rate-shift \n (maximum-a-posteriori estimate)", cex.main=0.9, cex.axis=0.9, cex.lab=0.9, cex.sub=0.9, xlim = c(minX, maxX), col="black", freq=FALSE, breaks=histbreaks, xaxp=c(minX, maxX, ticks))
	d <- (density(tdat2[,n]))
	new.dy <- d$y*((max(ahist$density)/max(d$y))*0.66)
	lines(d$x,new.dy, col="red", lwd=1.5)
	}
dev.off()
cat("  results are saved in: ", sprintf("%s", output), "\n\n")
}

rate.shift(as.character(arg[1]))