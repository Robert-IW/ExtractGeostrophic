# Extract the geostrophic velocities from 25 km AVISO netcdf files

# open netcdf as raster stack and extract data using the wind transects
# create a direction and speed from the u and v vectors

# for local transect, only the nearest value is retained


library(raster)
library(rasterVis)
library(ncdf4)
library(Imap)
library(dplyr)
library(data.table)

setwd("/media/robert/KELP-HDD-Portable/Altimetry/AVISO/data/")
filelist <- list.files()

# path to wind transects
trans.path <- ("~/R/myFunctions/Wind_transects/")

N <- length(filelist)

vgosa.stack <- stack()
ugosa.stack <- stack()

vgos.stack <- stack()
ugos.stack <- stack()

dates <- vector()

for (i in 1:N){
  
  vgosa <- stack(filelist[i],varname="vgosa")
  dates <- c(dates, format(getZ(vgosa), format="%Y%m%d"))

  ugosa <- stack(filelist[i],varname="ugosa")

  vgos <- stack(filelist[i],varname="vgos")
  ugos <- stack(filelist[i],varname="ugos")

  vgosa.stack <- stack(vgosa.stack, vgosa)
  ugosa.stack <- stack(ugosa.stack, ugosa)

  vgos.stack <- stack(vgos.stack, vgos)
  ugos.stack <- stack(ugos.stack, ugos)

  gc()
}

# open the station names and locations
station.loc <- read.table("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/station_locations.Rdata",
                          stringsAsFactors=F)

# create empty list for all the days data
mylist <- sapply(station.loc[,3],function(x) NULL)

# For each station
for (h in 60:nrow(station.loc)){
  
  stat.name <- station.loc$station[h]
  st.lon <- station.loc$lon[h]
  st.lat <- station.loc$lat[h]
  
  cat(paste("Starting on",stat.name,"transects\n"))
  
  # this produces 'station.list' with each element a transect
  file.load <- paste0(trans.path,stat.name,"_transects.Rdata")
  load(file.load)
  rm(file.load)
  
  # get the individual transects from the list
  trans.along <- station.list$alongshore
  trans.offshore <- station.list$offshore
  region.wind <- station.list$local
  rm(station.list)

  # ------------------------------------------------------------------ ALONGSHORE
  cat("Working on ALONGSHORE transects\n")
  
  # extract the data from the stack
  vgosa.along <- extract(vgosa.stack, trans.along, cellnumbers=TRUE)
  ugosa.along <- extract(ugosa.stack, trans.along, cellnumbers=TRUE)
  
  vgos.along <- extract(vgos.stack, trans.along, cellnumbers=TRUE)
  ugos.along <- extract(ugos.stack, trans.along, cellnumbers=TRUE)
  
  # get the coordinates of the extracted points
  along.loc <- xyFromCell(vgosa.stack[[1]], vgosa.along[[1]][,'cell'])   # a matrix
  colnames(along.loc) <- c("lon","lat")
  
  # remove first column from matrix, which contains the cell numbbers
  vgosa.along <- vgosa.along[[1]][,-1]
  ugosa.along <- ugosa.along[[1]][,-1]
  
  vgos.along <- vgos.along[[1]][,-1]
  ugos.along <- ugos.along[[1]][,-1]
  
  # create empty list for each daily transect
  vgosa.data <- vector("list",ncol(vgosa.along))
  ugosa.data <- vector("list",ncol(ugosa.along))
  
  vgos.data <- vector("list",ncol(vgos.along))
  ugos.data <- vector("list",ncol(ugos.along))
  
  for (j in 1:ncol(vgosa.along)){
    vgosa.data[[j]] <- cbind(dates[j],along.loc, vgosa.along[,j])
    ugosa.data[[j]] <- cbind(dates[j],along.loc, ugosa.along[,j])
    
    vgos.data[[j]] <- cbind(dates[j],along.loc, vgos.along[,j])
    ugos.data[[j]] <- cbind(dates[j],along.loc, ugos.along[,j])
  }
  rm(vgosa.along, ugosa.along,vgos.along, ugos.along)
  
  # create matrix with all the data
  vgosa.matrix <- do.call(rbind, vgosa.data)
  colnames(vgosa.matrix) <- c("date","lon","lat","vanom")
  vgosa.matrix <- apply(vgosa.matrix, 2, as.numeric)
  
  ugosa.matrix <- do.call(rbind, ugosa.data)
  colnames(ugosa.matrix) <- c("date","lon","lat","uanom")
  ugosa.matrix <- apply(ugosa.matrix, 2, as.numeric)
  
  vgos.matrix <- do.call(rbind, vgos.data)
  colnames(vgos.matrix) <- c("date","lon","lat","vgeos")
  vgos.matrix <- apply(vgos.matrix, 2, as.numeric)
  
  ugos.matrix <- do.call(rbind, ugos.data)
  colnames(ugos.matrix) <- c("date","lon","lat","ugeos")
  ugos.matrix <- apply(ugos.matrix, 2, as.numeric)
  
  rm(vgosa.data, ugosa.data, vgos.data, ugos.data)
  
  # coerce to data frame and change character to numeric and add 'station' and 'class'
  vgosa.df.along <- as.data.frame(vgosa.matrix, stringsAsFactors=F)
  vgosa.df.along <- cbind("station"=stat.name, "class"="alongshore", vgosa.df.along)
  ugosa.df.along <- as.data.frame(ugosa.matrix, stringsAsFactors=F)
  ugosa.df.along <- cbind("station"=stat.name, "class"="alongshore", ugosa.df.along)
  
  vgos.df.along <- as.data.frame(vgos.matrix, stringsAsFactors=F)
  vgos.df.along <- cbind("station"=stat.name, "class"="alongshore", vgos.df.along)
  ugos.df.along <- as.data.frame(ugos.matrix, stringsAsFactors=F)
  ugos.df.along <- cbind("station"=stat.name, "class"="alongshore", ugos.df.along)
  
  rm(vgosa.matrix, ugosa.matrix, vgos.matrix, ugos.matrix)
  
  # --------------------------------------------------------------------------------- OFFSHORE
  cat("Working on OFFSHORE transects\n")
  
  # extract the data from stack
  vgosa.off <- extract(vgosa.stack, trans.offshore, cellnumbers=TRUE)
  ugosa.off <- extract(ugosa.stack, trans.offshore, cellnumbers=TRUE)
  
  vgos.off <- extract(vgos.stack, trans.offshore, cellnumbers=TRUE)
  ugos.off <- extract(ugos.stack, trans.offshore, cellnumbers=TRUE)
  
  # get the coordinates of the extracted points
  off.loc <- xyFromCell(vgosa.stack[[1]], vgosa.off[[1]][,'cell'])   # a matrix
  colnames(off.loc) <- c("lon","lat")
  
  # remove first column from matrix, which contains the cell numbbers
  vgosa.off <- vgosa.off[[1]][,-1]
  ugosa.off <- ugosa.off[[1]][,-1]
  
  vgos.off <- vgos.off[[1]][,-1]
  ugos.off <- ugos.off[[1]][,-1]
  
  # create empty list for each daily transect
  vgosa.data <- vector("list",ncol(vgosa.off))
  ugosa.data <- vector("list",ncol(ugosa.off))
  
  vgos.data <- vector("list",ncol(vgos.off))
  ugos.data <- vector("list",ncol(ugos.off))
  
  for (j in 1:ncol(vgosa.off)){
    vgosa.data[[j]] <- cbind(dates[j],off.loc, vgosa.off[,j])
    ugosa.data[[j]] <- cbind(dates[j],off.loc, ugosa.off[,j])
    
    vgos.data[[j]] <- cbind(dates[j],off.loc, vgos.off[,j])
    ugos.data[[j]] <- cbind(dates[j],off.loc, ugos.off[,j])
  }
  rm(vgosa.off, ugosa.off,vgos.off, ugos.off)
  
  # create matrix with all the data
  vgosa.matrix <- do.call(rbind, vgosa.data)
  colnames(vgosa.matrix) <- c("date","lon","lat","vanom")
  vgosa.matrix <- apply(vgosa.matrix, 2, as.numeric)
  
  ugosa.matrix <- do.call(rbind, ugosa.data)
  colnames(ugosa.matrix) <- c("date","lon","lat","uanom")
  ugosa.matrix <- apply(ugosa.matrix, 2, as.numeric)
  
  vgos.matrix <- do.call(rbind, vgos.data)
  colnames(vgos.matrix) <- c("date","lon","lat","vgeos")
  vgos.matrix <- apply(vgos.matrix, 2, as.numeric)
  
  ugos.matrix <- do.call(rbind, ugos.data)
  colnames(ugos.matrix) <- c("date","lon","lat","ugeos")
  ugos.matrix <- apply(ugos.matrix, 2, as.numeric)
  
  rm(vgosa.data,ugosa.data,vgos.data,ugos.data)
  
  # coerce to data frame and change character to numeric and add 'station' and 'class'
  vgosa.df.off <- as.data.frame(vgosa.matrix, stringsAsFactors=F)
  vgosa.df.off <- cbind("station"=stat.name, "class"="offshore", vgosa.df.off)
  
  ugosa.df.off <- as.data.frame(ugosa.matrix, stringsAsFactors=F)
  ugosa.df.off <- cbind("station"=stat.name, "class"="offshore", ugosa.df.off)
  
  vgos.df.off <- as.data.frame(vgos.matrix, stringsAsFactors=F)
  vgos.df.off <- cbind("station"=stat.name, "class"="offshore", vgos.df.off)
  
  ugos.df.off <- as.data.frame(ugos.matrix, stringsAsFactors=F)
  ugos.df.off <- cbind("station"=stat.name, "class"="offshore", ugos.df.off)
  
  rm(vgosa.matrix, ugosa.matrix, vgos.matrix, ugos.matrix)
  
  # --------------------------------------------------------------------------------- LOCAL
  # extract the data from the stack: produces list with matrix 3 rows
  cat("Working on LOCAL transects\n")
  
  vgosa.local <- extract(vgosa.stack, region.wind, cellnumbers=TRUE)
  
  # get the coordinates of the extracted points and find the nearest
  local.loc <- xyFromCell(vgosa.stack[[1]], vgosa.local[[1]][,'cell'])   # a matrix
  
  # remove first column from matrix, which contains the cell numbbers: produces 3 x length matrix

  vgosa.local <- extract(vgosa.stack, region.wind, cellnumbers=TRUE)
  vgosa.local <- vgosa.local[[1]][,-1]
  
  ind.na <- which(is.na(rowSums(vgosa.local)))  # get index of rows all na
  
  # if all rows are NA extent the bbox of region.wind
  if (nrow(vgosa.local) == length(ind.na)){
    cnt <- 1
    poly.ext <- region.wind@bbox
    poly.lon <- c(poly.ext[1,1],poly.ext[1,2],poly.ext[1,2],poly.ext[1,1],poly.ext[1,1])
    poly.lat <- c(poly.ext[2,2],poly.ext[2,2],poly.ext[2,1],poly.ext[2,1],poly.ext[2,2])
    
    poly.mat <- matrix(c(poly.lon,poly.lat),ncol=2)
    poly.adj <- matrix(c(-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5),ncol=2)

    poly.poly <- Polygon(poly.mat)
    poly.sp <- SpatialPolygons(list(Polygons(list(poly.poly),ID=1)))

    while (nrow(vgosa.local) == length(ind.na)){
      cat(paste("Working on new bbox, expanded=",cnt,"\n"))
      
      vgosa.local <- extract(vgosa.stack, poly.sp, cellnumbers=TRUE)
      local.loc <- xyFromCell(vgosa.stack[[1]], vgosa.local[[1]][,'cell'])   # a matrix
      vgosa.local <- vgosa.local[[1]][,-1]
      ind.na <- which(is.na(rowSums(vgosa.local)))  # get index of rows all na
      cnt <- cnt+1
      poly.poly <- Polygon(poly.mat+poly.adj)
      poly.sp <- SpatialPolygons(list(Polygons(list(poly.poly),ID=1)))
    }
    # extract remaining variables with poly.sp
    ugosa.local <- extract(ugosa.stack, poly.sp, cellnumbers=TRUE)
    vgos.local <- extract(vgos.stack, poly.sp, cellnumbers=TRUE)
    ugos.local <- extract(ugos.stack, poly.sp, cellnumbers=TRUE)
  } else {
    # extract remaining variables with region.wind
    ugosa.local <- extract(ugosa.stack, region.wind, cellnumbers=TRUE)
    vgos.local <- extract(vgos.stack, region.wind, cellnumbers=TRUE)
    ugos.local <- extract(ugos.stack, region.wind, cellnumbers=TRUE)
  }
  
  ugosa.local <- ugosa.local[[1]][,-1]
  vgos.local <- vgos.local[[1]][,-1]
  ugos.local <- ugos.local[[1]][,-1]
  
  # estimate the distance from the station coords
  local.dist <- vector()
  for (k in 1:nrow(local.loc)){
    local.dist[k] <- gdist(lat.1=st.lat,
                           lon.1=st.lon,
                           lat.2=local.loc[k,2],
                           lon.2=local.loc[k,1],
                           units = "km")
  }
  
  # if there are NA rows but not all NA, remove
  if (length(ind.na)){
    local.dist <- local.dist[-ind.na]             # remove distance corresponding to na
    local.loc <- local.loc[-ind.na, ,drop=F]
    
    vgosa.local <- vgosa.local[-ind.na, ,drop=F]          # remove rows all na
    ugosa.local <- ugosa.local[-ind.na, ,drop=F]          # remove rows all na
    vgos.local <- vgos.local[-ind.na, ,drop=F]            # remove rows all na
    ugos.local <- ugos.local[-ind.na, ,drop=F]            # remove rows all na
  }
  
  ind.near <- which(local.dist==min(local.dist))  # get index of nearest dist
  local.loc <- local.loc[ind.near, ,drop=F]
  colnames(local.loc) <- c("lon","lat")
  
  vgosa.local <- vgosa.local[ind.near,,drop=F]         # keep nearest
  ugosa.local <- ugosa.local[ind.near,,drop=F]         # keep nearest
  vgos.local <- vgos.local[ind.near,,drop=F]           # keep nearest
  ugos.local <- ugos.local[ind.near,,drop=F]           # keep nearest
  
  # create empty list for each daily transect
  vgosa.data <- vector("list",length(vgosa.local))
  ugosa.data <- vector("list",length(ugosa.local))
  
  vgos.data <- vector("list",length(vgos.local))
  ugos.data <- vector("list",length(ugos.local))
  
  for (j in 1:length(vgosa.local)){
    vgosa.data[[j]] <- cbind(dates[j],local.loc, vgosa.local[j])
    ugosa.data[[j]] <- cbind(dates[j],local.loc, ugosa.local[j])
    
    vgos.data[[j]] <- cbind(dates[j],local.loc, vgos.local[j])
    ugos.data[[j]] <- cbind(dates[j],local.loc, ugos.local[j])
  }
  
  rm(vgosa.local, ugosa.local,vgos.local, ugos.local)
  
  # create matrix with all the data and covert data to numeric
  vgosa.matrix <- do.call(rbind, vgosa.data)
  colnames(vgosa.matrix) <- c("date","lon","lat","vanom")
  vgosa.matrix <- apply(vgosa.matrix, 2, as.numeric)
  
  ugosa.matrix <- do.call(rbind, ugosa.data)
  colnames(ugosa.matrix) <- c("date","lon","lat","uanom")
  ugosa.matrix <- apply(ugosa.matrix, 2, as.numeric)
  
  vgos.matrix <- do.call(rbind, vgos.data)
  colnames(vgos.matrix) <- c("date","lon","lat","vgeos")
  vgos.matrix <- apply(vgos.matrix, 2, as.numeric)
  
  ugos.matrix <- do.call(rbind, ugos.data)
  colnames(ugos.matrix) <- c("date","lon","lat","ugeos")
  ugos.matrix <- apply(ugos.matrix, 2, as.numeric)
  
  rm(vgosa.data,ugosa.data,vgos.data,ugos.data)
  
  cat("Creating data frame for LOCAL\n")
  
  # coerce to data frame and change character to numeric and add 'station' and 'class'
  vgosa.df.local <- as.data.frame(vgosa.matrix, stringsAsFactors=F)
  vgosa.df.local <- cbind("station"=stat.name, "class"="local", vgosa.df.local)
  ugosa.df.local <- as.data.frame(ugosa.matrix, stringsAsFactors=F)
  ugosa.df.local <- cbind("station"=stat.name, "class"="local", ugosa.df.local)
  
  vgos.df.local <- as.data.frame(vgos.matrix, stringsAsFactors=F)
  vgos.df.local <- cbind("station"=stat.name, "class"="local", vgos.df.local)
  ugos.df.local <- as.data.frame(ugos.matrix, stringsAsFactors=F)
  ugos.df.local <- cbind("station"=stat.name, "class"="local", ugos.df.local)
  
  rm(vgosa.matrix, ugosa.matrix, vgos.matrix, ugos.matrix)
  
  # combine all transects for each variable and save in list of stations
  A <- bind_rows(vgosa.df.along,vgosa.df.off,vgosa.df.local)
  B <- bind_rows(ugosa.df.along,ugosa.df.off,ugosa.df.local)
  C <- bind_rows(vgos.df.along,vgos.df.off,vgos.df.local)
  D <- bind_rows(ugos.df.along,ugos.df.off,ugos.df.local)
  
  df <- bind_rows(A,B,C,D)
  rm(A,B,C,D)
  
  mylist[[h]] <- df
  rm(df)
  gc()
}

# for each in the list save
for (k in 1:length(mylist)){
  
  stat <- station.loc[k,3]
  
  cat(paste("Working on saving", stat, "station number", k, "\n"))
  
  df <- mylist[[k]]
  dir.save <- "/media/robert/KELP-HDD-Portable/Altimetry/AVISO/ESdata/"
  filename1 <- paste0(stat, "_geostropic.csv",sep='')
  
  data.table::fwrite(df, file = filename1)
  rm(df, stat, filename1)
  gc()
}
  
  
