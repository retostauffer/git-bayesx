plotsurf<-function(data, outfile, x, y, z, mode = 1, cols = c(2, 3, 4), ticktype="detailed", expand=0.5, d=200, theta=-40, phi=20, ...)
{
	require(akima)
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

	data <- interp(x, y, z)

	if(!missing(outfile))
		postscript(outfile,horizontal=F)

	if(mode==1)
		persp(data, ticktype=ticktype, expand=expand, d=d, theta=theta, phi=phi, ...)
	else if(mode==2)
		contour(data,...)
	else if(mode==3)
		image(data,...)

	if(!missing(outfile))
		dev.off()

	return(invisible())
}

