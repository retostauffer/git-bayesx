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
