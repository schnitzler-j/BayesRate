#!/usr/bin/Rscript

arg <- commandArgs(trailingOnly=TRUE)

subsample <- function(file = NULL, burnin = 0.25, samples = 10, random=TRUE)
{
	
if (!is.null(file)) 
    filey <- paste(file, sep="")
else stop("you must enter a filename or a character string\n")

if (is.null(samples))
	stop ("please enter the number of trees you would like to sample")

bi<-burnin
no<-samples

x <- scan(filey, what=character(), sep="\n")
tree.subsamp <- NULL

range <- grep("tree STATE|tree rep.|tree gen.|TREE \\* UNTITLED|TREE UNTITLED", x)    
end <- grep("End;|END;|end;", x)

if (bi < 1){
	bi <- bi * 100
	burn <- ceiling((((max(end)) - min(range))*bi)/100)
	} else
	{ burn <- bi } 

if (no >= length(range)-burn)
	stop ("error: there are not enough trees available")

if (random == FALSE) {
	freq <- floor((length(range)-burn)/no)
	tree <- seq(from=min(range)+burn+1, to=max(end)-1, by=freq)
	} else 
	{ tree <- sample((min(range)+burn+1):(max(end)-1), no, replace=FALSE) }

header <- x[1:(min(range)-1)]
footer <- x[c(max(end):length(x))]

tree.subsamp <- header
tree.subsamp <- c(tree.subsamp,x[tree])
tree.subsamp <- c(tree.subsamp,footer)

if (grepl("\\.tre", file))
	{ file2 <- gsub ("\\.tre", "_mod.tre", file) 
} else if (grepl("\\.t$", file))
		{ file2 <- gsub ("\\.t$", "_mod.t", file) 
	} else if (grepl("\\.nex$", file))
		{ file2 <- gsub ("\\.nex$", "_mod.nex", file) 
		} else 
		{ file2 <- paste(file,"mod") }
	
write(tree.subsamp, file=file2)
cat("  The resampled trees are saved as: ", sprintf("%s", file2), "\n\n")

range <- NULL
end <- NULL
tree <- NULL
file2 <- NULL
}

subsample(as.character(arg[1]), as.double(arg[2]), as.integer(arg[3]), as.logical(arg[4]))

