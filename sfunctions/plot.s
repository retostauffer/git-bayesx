"plotnonp"<-
function(file, psname, month, year, step = 12, ylimtop, ylimbottom, xlab, ylab, 
	maintitle = "", subtitle = "", type = "l", linetype = 1, linecol = 3, 
	level = 0, mfrow, ...)
{
	if(!missing(psname)) {
		if(file == psname)
			stop("Argument \"file\" and argument \"psname\" are identical"
				)
		else postscript(psname, horizontal = F)
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
"plotsample"<-
function(file, psname, type = "l", ...)
{
	if(!missing(psname)) {
		if(file == psname)
			stop("Argument \"file\" and argument \"psname\" are identical"
				)
		else postscript(psname, horizontal = F)
	}
	if(exists("is.R") && is.function(is.R) && is.R())
		old.par <- par(no.readonly = T)
	else old.par <- par()
	data <- read.table(file, header = T)
	i <- 1
	count <- 0
	while(i < length(data)) {
		if(count %% 6 == 0) {
			count <- 0
			if(exists("is.R") && is.function(is.R) && is.R() && 
				missing(psname))
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
"drawmap"<-
function(map, dfile, outfile, regionvar, plotvar, lowerlimit, upperlimit, 
	nrcolors = 100, pstitle = "", color = F, legend = T, drawnames = F, 
	swapcolors = F, pcat = F)
{
	if(!missing(outfile) && !missing(dfile)) {
		if(dfile == outfile)
			stop("Argument \"dfile\" and argument \"outfile\" are identical"
				)
	}
	if(pcat)
		nrcolors <- 3
	black <- 1
	xmin <- 1:length(map)
	xmax <- 1:length(map)
	ymin <- 1:length(map)
	ymax <- 1:length(map)
	for(i in 1:length(map)) {
		xmin[i] <- min(map[[i]][, 1], na.rm = T)
		xmax[i] <- max(map[[i]][, 1], na.rm = T)
		ymin[i] <- min(map[[i]][, 2], na.rm = T)
		ymax[i] <- max(map[[i]][, 2], na.rm = T)
	}
	xlimits <- c(min(xmin), max(xmax))
	ylimits <- c(min(ymin) - (max(ymax) - min(ymin)) * 0.1, max(ymax))
	ratio.xy <- diff(range(ylimits))/diff(range(xlimits))
	if(missing(dfile) && missing(regionvar) && missing(plotvar)) {
# draw only borders
		if(missing(outfile)) win.graph() else postscript(outfile, 
				horizontal = F, colors = 0, width = 8, height
				 = 9.5 * ratio.xy)
		plot(xlimits, ylimits, type = "n", axes = F, col = 0, xlab = "",
			ylab = "")
		for(k in 1:length(map))
			polygon(map[[k]][, 1], map[[k]][, 2], density = 0, lwd
				 = 0.3, col = 1)
	}
	else {
# draw and fill borders
		if(missing(dfile)) {
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
		if(!missing(outfile)) {
# define colors
			black <- nrcolors + 1
			white <- nrcolors + 2
			if(color == T) {
				ps.colors <- matrix(1, nrow = white, ncol = 3)
				if(swapcolors == F)
				  ps.colors[1:nrcolors, 1] <- ((1:nrcolors) - 1
				    )/(3 * (nrcolors - 1))
				else ps.colors[1:nrcolors, 1] <- ((nrcolors:1) - 
				    1)/(3 * (nrcolors - 1))
				ps.colors[black,  ] <- c(0, 0, 0)
				ps.colors[white,  ] <- c(0, 0, 1)
			}
			else {
				if(swapcolors == F)
				  ps.colors <- ((1:nrcolors) - 1)/(nrcolors - 1
				    )
				else ps.colors <- ((nrcolors:1) - 1)/(nrcolors - 
				    1)
				ps.colors[black] <- 0
				ps.colors[white] <- 1
			}
			postscript(outfile, horizontal = F, colors = ps.colors, 
				width = 8, height = 9.5 * ratio.xy)
			legend.colors <- 1:nrcolors
		}
		else {
			screen.colors <- as.integer(seq(0, 255, length = 2))
			screen.colors <- matrix(screen.colors, ncol = 3, nrow
				 = 2)
			black <- 1
			white <- 2
			if(color == T) {
				if(nrcolors %% 2 == 1) {
				  if(swapcolors == F)
				    graphsheet(color.table = screen.colors, 
				      num.image.colors = 3, num.image.shades = 
				      paste((nrcolors - 3)/2, ",", (nrcolors - 
				      3)/2), image.color.table = 
				      "255,0,0|255,255,0|0,255,0")
				  else graphsheet(color.table = screen.colors, 
				      num.image.colors = 3, num.image.shades = 
				      paste((nrcolors - 3)/2, ",", (nrcolors - 
				      3)/2), image.color.table = 
				      "0,255,0|255,255,0|255,0,0")
				}
				else {
				  shades <- paste(nrcolors/2 - 2, ",0,", 
				    nrcolors/2 - 2)
				  if(swapcolors == F)
				    colors <- paste("255,0,0|255,", as.integer((
				      nrcolors - 1)/nrcolors * 255), ",0|", 
				      as.integer((nrcolors - 1)/nrcolors * 255),
				      ",255,0|0,255,0", sep = "")
				  else colors <- paste("0,255,0|", as.integer((
				      nrcolors - 1)/nrcolors * 255), 
				      ",255,0|255,", as.integer((nrcolors - 1)/
				      nrcolors * 255), ",0|255,0,0", sep = "")
				  graphsheet(color.table = screen.colors, 
				    num.image.colors = 4, num.image.shades = 
				    shades, image.color.table = colors)
				}
			}
			else {
				if(swapcolors == F)
				  graphsheet(color.table = screen.colors, 
				    num.image.colors = 2, num.image.shades = 
				    paste(nrcolors - 2), image.color.table = 
				    "0,0,0|255,255,255")
				else graphsheet(color.table = screen.colors, 
				    num.image.colors = 2, num.image.shades = 
				    paste(nrcolors - 2), image.color.table = 
				    "255,255,255|0,0,0")
			}
			legend.colors <- 1:nrcolors + 16
		}
		maxim <- max(plotvar, na.rm = T)
		minim <- min(plotvar, na.rm = T)
		if(pcat) {
			upperlimit <- 1
			lowerlimit <- -1
		}
		else {
			if(missing(lowerlimit)) {
				lowerlimit <- minim
			}
			else if(lowerlimit > minim) {
				plotvar[plotvar < lowerlimit] <- lowerlimit
				cat(paste(
				  "Note: lowerlimit is above minimum value (", 
				  lowerlimit, " > ", minim, ")\n"))
			}
			if(missing(upperlimit)) {
				upperlimit <- maxim
			}
			else if(upperlimit < maxim) {
				plotvar[plotvar > upperlimit] <- upperlimit
				cat(paste(
				  "Note: upperlimit is below maximum value (", 
				  upperlimit, " < ", maxim, ")\n"))
			}
		}
		fill.colors <- cut(c(lowerlimit, plotvar, upperlimit), nrcolors
			)
		fill.colors <- fill.colors[c(-1,  - length(fill.colors))]
		if(missing(outfile))
			fill.colors <- fill.colors + 16
		plot(xlimits, ylimits, type = "n", axes = F, col = white, xlab
			 = "", ylab = "")
		if(sum(!is.na(match(names(map), regionvar))) == 0)
			warning("map probably doesn't match datafile")
		block1 <- c()
		block2 <- c()
		for(k in 1:length(map)) {
			if(is.na(map[[k]][1, 1]) && is.na(map[[k]][1, 2]))
				block2 <- c(block2, k)
			else block1 <- c(block1, k)
		}
		m <- match(names(map), regionvar)
		for(k in block1) {
			if(is.na(m[k])) {
				polygon(map[[k]][, 1], map[[k]][, 2], col = 
				  white, border = F)
				polygon(map[[k]][, 1], map[[k]][, 2], density
				   = 15, lwd = 0.3, col = black)
			}
			else {
				polygon(map[[k]][, 1], map[[k]][, 2], col = 
				  fill.colors[m[k]], border = F)
				polygon(map[[k]][, 1], map[[k]][, 2], col = 
				  black, density = 0)
			}
		}
		for(k in block2) {
			if(is.na(m[k])) {
				polygon(map[[k]][-1, 1], map[[k]][-1, 2], col
				   = white, border = F)
				polygon(map[[k]][-1, 1], map[[k]][-1, 2], 
				  density = 15, lwd = 0.3, col = black)
			}
			else {
				polygon(map[[k]][-1, 1], map[[k]][-1, 2], col
				   = fill.colors[m[k]], border = F)
				polygon(map[[k]][-1, 1], map[[k]][-1, 2], col
				   = black, density = 0)
			}
		}
		if(legend == T) {
			ylo <- yro <- ylimits[1] + (0.7 * (ylimits[2] - ylimits[
				1]))/11
			ylu <- yru <- ylimits[1] + (0.3 * (ylimits[2] - ylimits[
				1]))/11
			tylu <- tyru <- ylimits[1]
			xlu <- xlo <- xlimits[1]
			xru <- xro <- xlimits[1] + 0.3 * (xlimits[2] - xlimits[
				1])
			step <- (xru - xlu)/nrcolors
			for(i in 0:(nrcolors - 1)) {
				polygon(c(xlo + step * i, xlo + step * (i + 1), 
				  xlu + step * (i + 1), xlu + step * i), c(ylo, 
				  yro, yru, ylu), col = legend.colors[i + 1], 
				  border = F)
			}
			polygon(c(xlo, xro, xru, xlu), c(ylo, yro, yru, ylu), 
				density = 0, col = black, lwd = 0.3)
			text(xlu + 0.5 * step, tylu, lowerlimit, cex = 0.7, col
				 = black)
			text(xru - 0.5 * step, tyru, upperlimit, cex = 0.7, col
				 = black)
			if(lowerlimit + (upperlimit - lowerlimit)/3 < 0 && 0 < 
				upperlimit - (upperlimit - lowerlimit)/3) {
				help <- cut(c(0, lowerlimit, upperlimit), 
				  nrcolors)
				text(xlu + step * (help[1] - 0.5), tylu, "0", 
				  cex = 0.7, col = black)
			}
		}
	}
	if(drawnames == T) {
		xpos <- (xmin + xmax)/2
		ypos <- (ymin + ymax)/2
		text(xpos, ypos, labels = names(map), cex = 0.7, col = black)
	}
	title(main = pstitle, col = black)
	if(!missing(outfile))
		dev.off()
	return(invisible())
}
"readbndfile"<-
function(path = path, name = name)
{
	options(warn = -1)
	data.raw <- scan(path, what = list("", ""), sep = ",")
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
		if(is.na(data.numeric[[1]][i + 1]) && is.na(data.numeric[[2]][i +
			1])) {
			npoints <- npoints + 1
		}
		elem1 <- data.numeric[[1]][i + 1:npoints]
		elem2 <- data.numeric[[2]][i + 1:npoints]
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
	count <- length(unique(names(map)))
	cat("Note: map consists of ", count, "regions\n")
	assign(name, map, where = 1)
	options(warn = 1)
	return(invisible())
}
"plotsurf"<-
function(data, outfile, x, y, z, mode = 1, cols = c(2, 3, 4), zlim, ...)
{
	if(!missing(data)) {
		if(is.character(data))
			data <- read.table(data, header = T)
		x <- data[, cols[1]]
		y <- data[, cols[2]]
		z <- data[, cols[3]]
	}
	else {
		x <- x
		y <- y
		z <- z
	}
	if(mode == 1) {
		guiCreate("SurfacePlot", xValues = x, yValues = y, zValues = z, 
			NumXOutputGrids = "DataGrid", NumYOutputGrids = 
			"DataGrid", FillSurface = F, DrawMesh = T, Extrapolate
			 = F, ...)
	}
	else if(mode == 2) {
		if(!missing(outfile))
			postscript(outfile, colors = 0:100/100)
		data <- interp(x, y, z)
		if(missing(zlim))
			persp(data, ...)
		else persp(data, zlim = zlim, ...)
		if(!missing(outfile))
			dev.off()
	}
	else if(mode == 3) {
		if(!missing(outfile))
			postscript(outfile)
		data <- interp(x, y, z)
		if(missing(zlim))
			image(data, ...)
		else image(data, zlim = zlim, ...)
		if(!missing(outfile))
			dev.off()
	}
	else if(mode == 4) {
		guiCreate("ContourPlot", xValues = x, yValues = y, zValues = z, 
			...)
	}
	else if(mode == 5) {
		guiPlot("16 Color Surface", DataSetValues = data[, cols], ...)
	}
	else if(mode == 0) {
	}
	else {
		cat("ERROR\n")
	}
	return(invisible())
}
