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
