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
