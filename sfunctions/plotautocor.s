"plotautocor" <-
function(file, psname, type = "l", mean.autocor = F, ...)
{
	limit <- (-0.144)
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
	if(mean.autocor == F) {
		data <- read.table(file, header = T)
		print <- seq(1, ncol(data), 1)
		print1 <- is.na(match(substring(names(data),nchar(names(data))-3),c("_min","_max")))
		print2 <- is.na(match(substring(names(data),nchar(names(data))-4),c("_mean")))
		print <- print[print1&print2]
		i <- 2
		count <- 0
		while(i <= length(print)) {
			if(count %% 6 == 0 || (print[i] - print[i - 1]) > 1) {
				count <- 0
				if(exists("is.R") && is.function(is.R) && is.R(
				  ) && missing(psname) && i>2)
				  windows()
				par(mfrow = c(3, 2))
				par(omi = c(1.5, 1, 0, 1.5)/2.54)
			}
			data[, print[i]][data[, print[i]] < limit] <- limit
			if(sum(!is.na(data[, print[i]])))
				plot(data[, 1], data[, print[i]], type = type, 
				  xlab = "lag", ylab = dimnames(data)[[2]][
				  print[i]], ylim = c(-0.1, 1), axes = F, ...)
			else plot(data[, 1], rep(1, dim(data)[1]), xlab = "lag",
				  ylab = dimnames(data)[[2]][print[i]], ylim = 
				  c(-0.1, 1), axes = F, pch = " ", ...)
			box()
			axis(1)
			axis(side = 2, at = seq(0, 1, 0.2), label = seq(0, 1, 
				0.2))
			abline(h = 0.1, lty = 2)
			count <- count + 1
			i <- i + 1
		}
	}
	else {
		data <- read.table(file, header = T)
		print <- seq(1, ncol(data), 1)
		print1 <- is.na(match(substring(names(data),nchar(names(data))-3),c("_min","_max")))
		print2 <- is.na(match(substring(names(data),nchar(names(data))-4),c("_mean")))
		print <- print[!(print1&print2)]
		i <- 1
		count <- 0
		while(i <= length(print)) {
			if(count %% 3 == 0) {
				count <- 0
				if(exists("is.R") && is.function(is.R) && is.R(
				  ) && missing(psname) && i>1)
				  windows()
				par(omi = c(2, 2, 0.5, 2.5)/2.54)
				par(mfrow = c(3, 1))
			}
			data[, print[i]][data[, print[i]] < limit] <- limit
			if(sum(!is.na(data[, print[i]])))
				plot(data[, 1], data[, print[i]], type = type, 
				  xlab = "lag", ylab = dimnames(data)[[2]][
				  print[i]], ylim = c(-0.1, 1), axes = F, ...)
			else plot(data[, 1], rep(1, dim(data)[1]), xlab = "lag",
				  ylab = dimnames(data)[[2]][print[i]], ylim = 
				  c(-0.1, 1), axes = F, pch = " ", ...)
			box()
			axis(1)
			axis(side = 2, at = seq(0, 1, 0.2), label = seq(0, 1, 
				0.2))
			abline(h = 0.1, lty = 2)
			count <- count + 1
			i <- i + 1
		}
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
