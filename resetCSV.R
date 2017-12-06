# Adjust the csv output of extractBrick.R

# Geostrophic
dir.main <- "/media/robert/KELP-HDD-Portable/"
dir.geos <- paste0(dir.main, "Altimetry/AVISO/ESdata/")

file.list <- list.files(dir.geos, full.names = T)

reset <- function(temp.sub){
  s1 <- (nrow(temp.sub)/4)
  s2 <- (nrow(temp.sub)/4)*2
  s3 <- (nrow(temp.sub)/4)*3
  
  temp.subA <- temp.sub[1:s1,]
  temp.subB <- temp.sub[(s1+1):s2,]
  temp.subC <- temp.sub[(s2+1):s3,]
  temp.subD <- temp.sub[(s3+1):nrow(temp.sub),]
  
  temp.new <- cbind(temp.subA[,1:6],temp.subB[,7],temp.subC[,8],temp.subD[,9])
  return(temp.new)
}

for (i in 1:length(file.list)){
  temp <- fread(file.list[i])
  
  temp.sub <- subset(temp, temp$class=="alongshore")
  al <- reset(temp.sub)
  rm(temp.sub)

  temp.sub <- subset(temp, temp$class=="offshore")
  of <- reset(temp.sub)
  rm(temp.sub)
  
  temp.sub <- subset(temp, temp$class=="local")
  lo <- reset(temp.sub)
  rm(temp.sub)
  
  new.df <- rbind(al,of,lo)
  
  cat(paste("Writing file",file.list[i],"\n"))
  data.table::fwrite(new.df, file = file.list[i])
}
