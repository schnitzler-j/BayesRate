#!/usr/bin/Rscript

arg <- commandArgs(trailingOnly=TRUE)

rate_summary <- function(filepath)
{

require(MASS, quietly=TRUE, warn.conflicts=FALSE)
require(TeachingDemos, quietly=TRUE, warn.conflicts=FALSE)
require(hdrcde, quietly=TRUE, warn.conflicts=FALSE)

cat("  reading log-file...\n")

tdat <- read.table(filepath, header=T)
headernames <- names(tdat)

if (any(headernames == "net.diversification")) {
	a <- grep("speciation.rate", headernames)
	b <- grep("extinction.rate", headernames)	
	c <- grep("net.diversification", headernames)
	d <- grep("extinction.fraction", headernames)
	} else if (any(headernames == "sp.rate_1")) {
		a <- grep("^sp.rate_", headernames)
		b <- rep(0, length(a))
		c <- rep(0, length(a))
		d <- rep(grep("ext.frac", headernames), length(a))
		} else {
			a <- grep("^speciation_", headernames)
			b <- grep("^extinction_", headernames)
			c <- grep("^net_div_", headernames)
			d <- grep("^ext_frac_", headernames)
	}

if (any(headernames == "net.diversification" | headernames == "sp.rate_1")) {
	clade.name <- NULL
	} else {
	clade.name <- headernames[c]
	clade.name <- gsub("net_div_", "", clade.name)
	}

df <- data.frame(a,b,c,d)

l <- c("speciation rate", "extinction rate", "net diversification", "extinction fraction")
cat("  computing HPDs and printing pdf...\n")

output <- gsub ("\\.log$", ".pdf", filepath)

pdf (file=output)

# check if all models are pb
if (sum(tdat[,d]) == 0) {
	for (n in 1:length(a)) {
			var_name <- l[1]
			cred_int <- emp.hpd (tdat[,a[n]], conf=0.95)
			map <- round(hdr(tdat[,a[n]], all.modes=F)$mode, digits=3)
			minX <- floor(10*min(tdat[,a[n]]))/10
			maxX <- ceiling(10*max(tdat[,a[n]]))/10
			ahist <- hist(tdat[,a[n]], plot=F)
			rate <- round(mean(tdat[,a[n]]), digits=3)			
			if (length(a) > 1) {
				hist(tdat[,a[n]], xlab = paste(sprintf("%s", var_name), "- mean:", rate[1], "/ MAP:", map), sub = paste("95% HPD:", sprintf("%.3f-%.3f", cred_int[1], cred_int[2])), main = sprintf("rate %s", n), ylim = c(0, max(ahist$counts+(ahist$counts/3))), xlim = c(minX, maxX))
			} else {
				hist(tdat[,a[n]], xlab = paste(sprintf("%s", var_name), "- mean:", rate[1], "/ MAP:", map), sub = paste("95% HPD:", sprintf("%.3f-%.3f", cred_int[1], cred_int[2])), main = "speciation rate", ylim = c(0, max(ahist$counts+(ahist$counts/3))), xlim = c(minX, maxX))
			}
		}
		
	} else {

par(mfrow=c(2,2))
for (n in 1:length(a)) {
	if (sum(tdat[,d[n]]) == 0) {
		var_name <- l[1]
		cred_int <- emp.hpd (tdat[,a[n]], conf=0.95)
		map <- round(hdr(tdat[,a[n]], all.modes=F)$mode, digits=3)
		minX <- floor(10*min(tdat[,a[n]]))/10
		maxX <- ceiling(10*max(tdat[,a[n]]))/10
		ahist <- hist(tdat[,a[n]], plot=F)
		rate <- round(mean(tdat[,a[n]]), digits=3)
		hist(tdat[,a[n]], xlab = paste(sprintf("%s", var_name), "- mean:", rate[1], "/ MAP:", map), sub = paste("95% HPD:", sprintf("%.3f-%.3f", cred_int[1], cred_int[2])), main = sprintf("clade: %s", clade.name[n]), ylim = c(0, max(ahist$counts+(ahist$counts/3))), xlim = c(minX, maxX))
	for (m in 1:3){
		frame()
		}
		
	} else {
		
	j = 1
	coln <- as.vector(df[n,])
		for(i in coln){
			var_name <- l[j]
			cred_int <- emp.hpd (tdat[,i], conf=0.95)
			map <- round(hdr(tdat[,i], all.modes=F)$mode, digits=3)
			minX <- floor(10*min(tdat[,i]))/10
			maxX <- ceiling(10*max(tdat[,i]))/10
				if(j==4){
					minX=0
					maxX=1}
			ahist <- hist(tdat[,i], plot=F)
			rate <- round(mean(tdat[,i]), digits=3)
			hist(tdat[,i], xlab = paste(sprintf("%s", var_name), "- mean:", rate[1], "/ MAP:", map), sub = paste("95% HPD:", sprintf("%.3f-%.3f", cred_int[1], cred_int[2])), main = sprintf("clade: %s", clade.name[n]), ylim = c(0, max(ahist$counts+(ahist$counts/3))), xlim = c(minX, maxX), cex.lab=0.8, cex.sub=0.8, cex.main=0.9, cex.axis=0.9)
	j <- j+1
			}
		}
	}
}

dev.off()
cat("  results are saved in: ", sprintf("%s", output), "\n\n")
}

rate_summary(as.character(arg[1]))