"drawmap" <-
function (map, dfile, outfile, regionvar, 
    plotvar, lowerlimit, upperlimit, nrcolors = 100, pstitle = "", 
    color = F, legend = T, drawnames = F, swapcolors = F, pcat = F, cex.legend=0.7, h = c(25, 130), c = 100, l = c(90, 70), hcl=F) 
{
    hcl.available <- hcl && require("vcd", character.only=TRUE)

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
            postscript(outfile, horizontal = F, width = 8, height = 8.5 * 
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
            if(hcl.available) {
                if (swapcolors == T) 
                    h <- rev(h)
                fill.colors <- diverge_hcl(nrcolors,h=h,c=c,l=l)[fill.colors]
                legend.colors <- diverge_hcl(nrcolors,h=h,c=c,l=l)
		}
            else {
                fill.colors <- (fill.colors-1)/(3*(nrcolors-1))
                if (swapcolors == T) 
                    fill.colors <- 1/3 - fill.colors
                fill.colors <- hsv(h = fill.colors)
                legend.colors <- hsv(h = (0:(nrcolors-1))/(3*(nrcolors-1)))
                if (swapcolors == T)
                   legend.colors <- rev(legend.colors)
                }
        }
        else {
            fill.colors <- (fill.colors-1)/(nrcolors-1)
            if (swapcolors == T) 
                fill.colors <- 1 - fill.colors
            fill.colors <- grey(fill.colors)
            legend.colors <- grey((0:(nrcolors-1))/(nrcolors-1))
            if (swapcolors == T)
               legend.colors <- rev(legend.colors)
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
            text(xlu + 0.5 * step, tylu, lowerlimit, cex = cex.legend, 
                col = black)
            text(xru - 0.5 * step, tyru, upperlimit, cex = cex.legend, 
                col = black)
            if (lowerlimit + (upperlimit - lowerlimit)/3 < 0 && 
                0 < upperlimit - (upperlimit - lowerlimit)/3) {
                help <- cut(c(0, lowerlimit, upperlimit), nrcolors)
                help <- as.vector(help, mode = "numeric")
                text(xlu + step * (help[1] - 0.5), tylu, "0", 
                  cex = cex.legend, col = black)
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

