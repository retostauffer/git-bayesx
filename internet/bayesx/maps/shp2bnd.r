shp2bnd <- function(shpname = "c:\\temp\\gm", regionnames = 1:16, bndname = "c:\\temp\\test.bnd", check.is.in = TRUE, replace = FALSE)
{
# TASK:
# 
# convert a shape-file into a boundary-file as expected by BayesX
# (see Ch. 5 of the Reference Manual for a detailled description of the expected file format)
#
#
# ARGUMENTS:
#
# shpname     = base filename of the shapefile (including path)
# regionnames = either a vector of region names or the name of the variable in the dbf-file reperesenting these names
# bndname     = target file to store the bnd-file in
# check.is.in = test whether some regions are surrounded by other regions (FALSE speeds up the execution time but may result in a corrupted bnd-file)
# replace     = replace an existing version of the bnd-file?
# 
# DEPENDS ON:
#
# shapefiles
# maptools
#
# EXAMPLES:
# 
# shp2bnd(shpname = "c:\\temp\\gm", regionnames = 1:16, bndname = "c:\\temp\\test.bnd")
# shp2bnd(shpname = "c:\\temp\\gm", regionnames = "ADMIN_NAME", bndname = "c:\\temp\\test.bnd")

# load library 'shapefiles'
require(shapefiles)
require(maptools)

# read the shapefile information
shp<-read.shapefile(shpname)
dbf<-read.dbf(paste(shpname,".dbf",sep=""))

# extract names of the regions

if(is.character(regionnames))
  {
  regions<-dbf[regionnames][[1]]
  }
else
  {
  regions <- regionnames
  }

# split data into closed polygons

regionsnew <- numeric(0)
origind <- numeric(0)
polynew <- list()
k <- 1

for(i in 1:length(regions))
  {
  temppoly <- shp$shp$shp[[i]]$points
  n <- nrow(temppoly)
  i1 <- 1
  i2 <- which(temppoly[,1]==temppoly[i1,1] & temppoly[,2]==temppoly[i1,2])[2]

  while(i2 < n)
    {
    tempname<-paste("\"",regions[i],"\",",i2,sep="")
    regionsnew <- c(regionsnew,tempname)
    origind <- c(origind,i)
    polynew[[k]] <- temppoly[i1:i2,]
    k <- k+1

    temppoly <- temppoly[-(i1:i2),]
    n <- nrow(temppoly)
    i1 <- 1
    i2 <- which(temppoly[,1]==temppoly[i1,1] & temppoly[,2]==temppoly[i1,2])[2]
    }

  tempname<-paste("\"",regions[i],"\",",i2,sep="")
  regionsnew <- c(regionsnew,tempname)
  origind <- c(origind,i)
  polynew[[k]] <- temppoly[i1:i2,]
  k <- k+1
  }

# check for regions contained in another region

if(check.is.in)
  {
  dims <- rep(0,length(regionsnew))
  rmcheck <- rep(FALSE,length(regionsnew))
  for(i in 1:length(regionsnew))
    {
    dims[i] <- nrow(polynew[[i]])
    }
  for(i in 1:length(regionsnew))
    {
    dimshelp <- which(dims==dims[i])
    if(length(dimshelp)>1)
      {
      for(j in 2:length(dimshelp))
        {
        test <- (sum((polynew[[i]]-polynew[[dimshelp[j]]][nrow(polynew[[dimshelp[j]]]):1,])^2) < 10e-6)
        if(test)
          {
          if(maptools:::.ringDirxy(polynew[[dimshelp[j]]])<0)
            {
            rmcheck[dimshelp[j]] <- TRUE
            regionsnew[i] <- paste(regionsnew[i],"\nis.in,\"",regions[origind[dimshelp[j]]],"\"",sep="")
            }
          else
            {
            rmcheck[i] <- TRUE
            regionsnew[dimshelp[j]] <- paste(regionsnew[dimshelp[j]],"\nis.in,\"",regions[origind[i]],"\"",sep="")
            }
          }
        }
      }
    }

  regionsnew <- regionsnew[-which(rmcheck)]
  polynew <- polynew[-which(rmcheck)]
  }

# write to the specified file

if(replace & file.exists(bndname))
  {
  test <- file.remove(bndname)
  }
if(!replace & file.exists(bndname))
  {
  stop("Specified bnd-file already exists!")
  }
for(i in 1:length(regionsnew))
  {
  write(regionsnew[i],bndname,append=T)
  write.table(polynew[[i]],bndname,append=T,col.names=F,row.names=F,sep=",",quote=F)
  }

return(invisible())

}

