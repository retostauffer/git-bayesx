"drawmap"<-
function(map, dfile, outfile, regionvar, plotvar, lowerlimit, upperlimit, 
	nrcolors = 100, pstitle = "", color = F, legend = T, drawnames = F, 
	swapcolors = F)
{
	if(!missing(outfile) && !missing(dfile)) {
		if(dfile == outfile)
			stop("Argument \"dfile\" and argument \"outfile\" are identical"
				)
	}
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
		if(missing(outfile)) motif() else postscript(outfile, 
				horizontal = F, colors = 0, width = 8, height
				 = 8 * ratio.xy)
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
			postscript(outfile, horizontal = F, ps.colors = ps.hsb2rgb(ps.colors), 
				width = 8, height = 8 * ratio.xy)
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
		if(missing(lowerlimit)) {
			lowerlimit <- minim
		}
		else if(lowerlimit > minim) {
			plotvar[plotvar < lowerlimit] <- lowerlimit
			cat(paste("Note: lowerlimit is above minimum value (", 
				lowerlimit, " > ", minim, ")\n"))
		}
		if(missing(upperlimit)) {
			upperlimit <- maxim
		}
		else if(upperlimit < maxim) {
			plotvar[plotvar > upperlimit] <- upperlimit
			cat(paste("Note: upperlimit is below maximum value (", 
				upperlimit, " < ", maxim, ")\n"))
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
