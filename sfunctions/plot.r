"plotnonp" <-
function(file, psname, month, year, step = 12, ylimtop, ylimbottom, xlab, ylab, 
	maintitle = "", subtitle = "", type = "l", linetype = 1, linecol = 3, 
	level = 0, mfrow, ...)
{
	if(!missing(psname)) {
		if(file == psname)
			stop("Argument \"file\" and argument \"psname\" are identical"
				)
		else postscript(psname,horizontal=F)
		par(mfrow = c(2, 1))
	}
	block <- read.table(file, header = T)
	block <- block[, 2:8]
	if(missing(xlab))
		xlab <- dimnames(block)[[2]][1]
	if(missing(ylab))
		ylab <- paste("f(", xlab, ")")
	if(!missing(mfrow))
		par(mfrow = mfrow)
	par(omi = c(2.5, 2, 2.5, 2.5)/2.54)
	if(missing(ylimbottom)) {
		if(level != 1 && level != 2)
			ylimbottom <- min(apply(block[, 2:7], 2, min, na.rm = T
				))
		else ylimbottom <- min(block[, 5 - level], na.rm = T)
	}
	if(missing(ylimtop)) {
		if(level != 1 && level != 2)
			ylimtop <- max(apply(block[, 2:7], 2, max, na.rm = T))
		else ylimtop <- max(block[, 5 + level], na.rm = T)
	}
	plot(block[, 1], block[, 2], type = type, xlab = xlab, ylab = ylab, 
		ylim = c(ylimbottom, ylimtop), main = maintitle, sub = subtitle,
		axes = F, ...)
	box()
	axis(2)
	if(!missing(month) & !missing(year)) {
		start <- block[1, 1] - month + 1
		stop <- max(block[, 1] + 1, na.rm = T)
		pos <- seq(start, stop, step)
		label <- (pos - pos[1])/step + year
		if(nrow(block) <= 24) {
			if(step == 12) {
				label2 <- c("Jan", "Feb", "Mar", "Apr", "May", 
				  "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", 
				  "Dec")
			}
			else if(step == 4) {
				label2 <- c("Jan", "Apr", "Jul", "Oct")
			}
			else if(step == 2) {
				label2 <- c("Jan", "Jul")
			}
			else {
				label2 <- F
			}
			label2 <- rep(label2, length.out = nrow(block) + month - 
				1)
			label2 <- label2[month:(nrow(block) + month - 1)]
			start2 <- block[1, 1]
			stop2 <- max(block[, 1], na.rm = T)
			pos2 <- seq(start2, stop2, 1)
			axis(side = 1, at = pos2, label = label2, mgp = c(3, 
				0.5, 0))
			axis(side = 1, at = pos, label = label, mgp = c(3, 1.5, 
				0))
		}
		else axis(side = 1, at = pos, label = label)
	}
	else axis(1)
	if(level != 1) {
		lines(block[, 1], block[, 3], col = linecol, lty = linetype)
		lines(block[, 1], block[, 7], col = linecol, lty = linetype)
	}
	if(level != 2) {
		lines(block[, 1], block[, 4], col = linecol, lty = linetype)
		lines(block[, 1], block[, 6], col = linecol, lty = linetype)
	}
	if(!missing(psname))
		dev.off()
	return(invisible())
}
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
		not.to.print <- c(grep("*min", dimnames(data)[[2]]), grep(
			"*max", dimnames(data)[[2]]), grep("*mean", dimnames(
			data)[[2]]))
		not.to.print <- sort(not.to.print)
		print <- seq(1, ncol(data), 1)
		print <- print[ - not.to.print]
		i <- 2
		count <- 0
		while(i <= length(print)) {
			if(count %% 6 == 0 || (print[i] - print[i - 1]) > 1) {
				count <- 0
				if(exists("is.R") && is.function(is.R) && is.R(
				  ) && missing(psname))
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
		print <- c(grep("*min", dimnames(data)[[2]]), grep("*max", 
			dimnames(data)[[2]]), grep("*mean", dimnames(data)[[2]]
			))
		print <- sort(print)
		i <- 1
		count <- 0
		while(i <= length(print)) {
			if(count %% 3 == 0) {
				count <- 0
				if(exists("is.R") && is.function(is.R) && is.R(
				  ) && missing(psname))
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
"drawmap" <-
function (map, dfile, outfile, regionvar, 
    plotvar, lowerlimit, upperlimit, nrcolors = 100, pstitle = "", 
    color = F, legend = T, drawnames = F, swapcolors = F, pcat = F) 
{
    if (!missing(outfile) && !missing(dfile)) {
        if (dfile == outfile) 
            stop("Argument \"dfile\" and argument \"outfile\" are identical")
    }
    if(pcat) 
        nrcolors <- 3
    black <- grey(0)
    white <- grey(1)
    xmin <- 1:length(map)
    xmax <- 1:length(map)
    ymin <- 1:length(map)
    ymax <- 1:length(map)
    for (i in 1:length(map)) {
        xmin[i] <- min(map[[i]][, 1], na.rm = T)
        xmax[i] <- max(map[[i]][, 1], na.rm = T)
        ymin[i] <- min(map[[i]][, 2], na.rm = T)
        ymax[i] <- max(map[[i]][, 2], na.rm = T)
    }
    xlimits <- c(min(xmin), max(xmax))
    ylimits <- c(min(ymin) - (max(ymax) - min(ymin)) * 0.1, max(ymax))
    ratio.xy <- diff(range(ylimits))/diff(range(xlimits))
    if (missing(dfile) && missing(regionvar) && missing(plotvar)) {
        if (missing(outfile)) 
            x11()
        else postscript(outfile, horizontal = F, width = 8, height = 9.5 * 
            ratio.xy)
        plot(xlimits, ylimits, type = "n", axes = F, col = white, 
            xlab = "", ylab = "")
        for (k in 1:length(map)) polygon(map[[k]][, 1], map[[k]][, 
            2], lwd = 0.3, col = white, border = black)
    }
    else {
        if (missing(dfile)) {
            ord <- order(regionvar)
            plotvar <- plotvar[ord]
            regionvar <- regionvar[ord]
        }
        else {
            data <- read.table(dfile, header = T, row.names = NULL)
            ord <- order(data[, regionvar])
            plotvar <- data[, plotvar]
            plotvar <- plotvar[ord]
            regionvar <- data[, regionvar]
            regionvar <- regionvar[ord]
        }
        if (missing(outfile)) {
            x11()
        }
        else {
            postscript(outfile, horizontal = F, width = 8, height = 9.5 * 
                ratio.xy)
        }
        maxim <- max(plotvar, na.rm = T)
        minim <- min(plotvar, na.rm = T)
	if(pcat) {
		upperlimit <- 1
		lowerlimit <- -1
	}
	else {
	        if (missing(lowerlimit)) {
        	    lowerlimit <- minim
		}
	        else if (lowerlimit > minim) {
        	    plotvar[plotvar < lowerlimit] <- lowerlimit
	            cat(paste("Note: lowerlimit is above minimum value (", 
        	        lowerlimit, " > ", minim, ")\n"))
	        }
	        if (missing(upperlimit)) {
        	    upperlimit <- maxim
	        }
        	else if (upperlimit < maxim) {
	            plotvar[plotvar > upperlimit] <- upperlimit
        	    cat(paste("Note: upperlimit is below maximum value (", 
	                upperlimit, " < ", maxim, ")\n"))
        	}
	}
        fill.colors <- cut(c(lowerlimit, plotvar, upperlimit), 
            nrcolors)
        fill.colors <- fill.colors[c(-1, -length(fill.colors))]
        fill.colors <- as.vector(fill.colors, mode = "numeric")
        if (color == T) {
            fill.colors <- (fill.colors - 1)/(3 * (nrcolors - 
                1))
            if (swapcolors == T) 
                fill.colors <- sum(range(fill.colors)) - fill.colors
            fill.colors <- hsv(h = fill.colors)
            legend.colors <- hsv(h = (0:(nrcolors - 1))/(3 * 
                (nrcolors - 1)))
        }
        else {
            fill.colors <- (fill.colors - 1)/(nrcolors - 1)
            if (swapcolors == T) 
                fill.colors <- max(fill.colors) - fill.colors + 
                  min(fill.colors)
            fill.colors <- grey(fill.colors)
            legend.colors <- grey((0:(nrcolors - 1))/(nrcolors - 
                1))
        }
        if (swapcolors == T) {
            legend.colors <- legend.colors[length(legend.colors):1]
        }
        plot(xlimits, ylimits, type = "n", axes = F, col = white, 
            xlab = "", ylab = "")
        if (sum(!is.na(match(names(map), regionvar))) == 0) 
            warning("map probably doesn't match datafile")
        block1 <- c()
        block2 <- c()
        for (k in 1:length(map)) {
            if (is.na(map[[k]][1, 1]) && is.na(map[[k]][1, 2])) 
                block2 <- c(block2, k)
            else block1 <- c(block1, k)
        }
        m <- match(names(map), regionvar)
        for (k in block1) {
            if (is.na(m[k])) {
	    polygon(map[[k]][, 1], map[[k]][, 2], col = white, border = F)
	    polygon(map[[k]][, 1], map[[k]][, 2], density  = 15, lwd = 0.3, col = black)
            }
            else {
                polygon(map[[k]][, 1], map[[k]][, 2], col = fill.colors[m[k]], 
                  border = black)
            }
        }
        for (k in block2) {
            if (is.na(m[k])) {
         	    polygon(map[[k]][-1, 1], map[[k]][-1, 2], col = white, border = F)
                polygon(map[[k]][-1, 1], map[[k]][-1, 2], density = 15, lwd = 0.3, col = black)
            }
            else {
                polygon(map[[k]][-1, 1], map[[k]][-1, 2], col = fill.colors[m[k]], 
                  border = black)
            }
        }
        if (legend == T) {
            ylo <- yro <- ylimits[1] + (0.7 * (ylimits[2] - ylimits[1]))/11
            ylu <- yru <- ylimits[1] + (0.3 * (ylimits[2] - ylimits[1]))/11
            tylu <- tyru <- ylimits[1]
            xlu <- xlo <- xlimits[1]
            xru <- xro <- xlimits[1] + 0.3 * (xlimits[2] - xlimits[1])
            step <- (xru - xlu)/nrcolors
            for (i in 0:(nrcolors - 1)) {
                polygon(c(xlo + step * i, xlo + step * (i + 1), 
                  xlu + step * (i + 1), xlu + step * i), c(ylo, 
                  yro, yru, ylu), col = legend.colors[i + 1], 
                  border = legend.colors[i + 1])
            }
            lines(c(xlo, xro, xru, xlu, xlo), c(ylo, yro, yru, 
                ylu, ylo), col = black)
            text(xlu + 0.5 * step, tylu, lowerlimit, cex = 0.7, 
                col = black)
            text(xru - 0.5 * step, tyru, upperlimit, cex = 0.7, 
                col = black)
            if (lowerlimit + (upperlimit - lowerlimit)/3 < 0 && 
                0 < upperlimit - (upperlimit - lowerlimit)/3) {
                help <- cut(c(0, lowerlimit, upperlimit), nrcolors)
                help <- as.vector(help, mode = "numeric")
                text(xlu + step * (help[1] - 0.5), tylu, "0", 
                  cex = 0.7, col = black)
            }
        }
    }
    if (drawnames == T) {
        xpos <- (xmin + xmax)/2
        ypos <- (ymin + ymax)/2
        text(xpos, ypos, labels = names(map), cex = 0.7, col = black)
    }
    title(main = pstitle, col = black)
    if (!missing(outfile)) 
        dev.off()
    return(invisible())
}
"readbndfile" <-
function(path = path, name = name)
{
	options(warn = -1)
	data.raw <- scan(path, what = list("", ""), sep = ",", quote = "")
	data.numeric <- list(as.numeric(data.raw[[1]]), as.numeric(data.raw[[2
		]]))
	anzkreise <- sum(is.na(data.numeric[[1]])) - sum(data.raw[[1]] == 
		"is.in")
	cat("Note: map consists of", anzkreise, "polygons\n")
	cat("Reading map ...\n")
	map <- list()
	i <- 1
	for(k in 1:anzkreise) {
		j <- 1
		npoints <- data.numeric[[2]][i]
		if(is.na(data.numeric[[1]][i + 1]) && is.na(data.numeric[[2]][i + 1])) {
			npoints <- npoints + 1
		}
		elem1 <- data.numeric[[1]][i+1:npoints]
		elem2 <- data.numeric[[2]][i+1:npoints]
		map[[k]] <- matrix(c(elem1, elem2), ncol = 2)
		names(map)[k] <- substring(data.raw[[1]][i], 2, nchar(data.raw[[
			1]][i]) - 1)
		i <- i + npoints + 1
	}
	if(sum(is.na(as.numeric(names(map)))) == 0) {
		map <- map[order(as.numeric(names(map)))]
		cat("Note: regions sorted by number\n")
	}
	else {
		map <- map[order(names(map))]
		cat("Note: regions sorted by name\n")
	}
	count<-length(unique(names(map)))
 	cat("Note: map consists of ", count,"regions\n")
	assign(name, map, pos=1)
	options(warn = 1)
	return(invisible())
}
