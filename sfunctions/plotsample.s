"plotsample" <-
function(file, psname, type = "l", ...)
{
	if(!missing(psname)) {
		if(file == psname)
			stop("Argument \"file\" and argument \"psname\" are identical"
				)
		else postscript(psname,horizontal=F)
	}
	if(exists("is.R") && is.function(is.R) && is.R())
		old.par <- par(no.readonly=T)	
	else
		old.par <- par()
	data <- read.table(file, header = T)
	i <- 1
	count <- 0
	while(i < length(data)) {
		if(count %% 6 == 0) {
			count <- 0
			if(exists("is.R") && is.function(is.R) && is.R() && missing(psname))
				windows()
			par(mfrow = c(3, 2), omi = c(1.5, 1, 0, 1.5)/2.54)
		}
		plot(data[, 1], data[, i + 1], type = type, xlab = "iteration", 
			ylab = dimnames(data)[[2]][i + 1], ...)
		count <- count + 1
		i <- i + 1
	}
	par(old.par)
	if(!missing(psname)) {
		dev.off()
	}
	else {
		if(exists("is.R") && is.function(is.R) && is.R())
			cat("Use 'graphics.off()' to close all graphic windows\n"
				)
	}
	return(invisible())
}
