library(org.Hs.eg.db)
library(jsonlite)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

commands <- c("mdsScaled", "geneScaled")
args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

#--------------------------------------------------------------#

mds_manifest <- "../manifests/os.mds.network.2016-06-23.manifest.json"
hg19_manifest <- "../manifests/os.hg19.manifest.json"

mds_scaled_dir <- "../data/networks/mds/scaled"
hg19_dir <- "../data/molecular/hg19"

#---------------------------------------------------------
get.new.collection.index <- function(datasetName, dataTypeName){
  
  if(nrow(Manifest) == 0) return(1)
  
  dataObj <- subset(Manifest, dataset == datasetName && dataType == dataTypeName)
  if(nrow(dataObj) == 0) return(1)
  
  return(nrow(dataObj$collections[[1]]) +1)
}
#---------------------------------------------------------
add.new.collection <- function(datasetName, dataTypeName, collection){
  
  if(nrow(Manifest) == 0){	
    newCollection <- data.frame(dataset=datasetName, dataType=dataTypeName)
    newCollection$collections <- list(collection)
    Manifest <<- newCollection
    return()
  }
  
  dataObj <- subset(Manifest, dataset == datasetName & dataType == dataTypeName)
  if(nrow(dataObj) == 1){
    Manifest[Manifest$dataset==datasetName & Manifest$dataType ==dataTypeName,"collections"] <<- list(rbind(dataObj$collections[[1]],collection))
    return()
  }
  if(nrow(dataObj) == 0){	
    newCollection <- data.frame(dataset=datasetName, dataType=dataTypeName)
    newCollection$collections <- list(collection)
    Manifest <<- rbind(Manifest, newCollection)
    return()
  }
  stop(printf("add.new.collection found %d instances of dataset %s and dataType %s", length(dataObj), datasetName, dataTypeName))
  
}
#----------------------------------------------------------------------------------------------------
### Save Function Takes A matrix/data.frame + Base File Path (w/o extension) & Writes to Disk In Multiple (optionally specified) Formats
os.data.save <- function(df, directory, file, format = c("tsv", "csv", "RData", "JSON")){

  if(!dir.exists(directory))
    dir.create(file.path(directory), recursive=TRUE)
  
  if(!grepl("/$", directory)) directory <- paste(directory, "/", sep="")
  outFile = paste(directory, file, sep="")

  
  # Write Tab Delimited
  if("tsv" %in% format)
    write.table(df, file=paste(outFile,".tsv", sep = ""), quote=F, sep="\t")
  
  # Write CSV Delimited
  if("csv" %in% format)
    write.csv(df, file=paste(outFile,".csv",sep = ""), quote = F)
  
  # Write RData File
  if("RData" %in% format)
    save(df, file=paste(outFile,".RData", sep = "") )
  
  # Write JSON File
  if("JSON" %in% format)
  write(toJSON(df, pretty=TRUE, digits=I(8)), file=paste(outFile,".json", sep = "") )
  
  # Return DataFrame For Chaining
  return(df)
}

#--------------------------------------------------------------#
getChromosomeOffsets <- function(chromosomes, chrLengths, pLength, scaleFactor=1000){
	
  pLength = pLength/scaleFactor
  chrLengths = chrLengths/scaleFactor
  
	yCent <- max(pLength)
	yOffset <- sapply(chromosomes, function(chr) {  yCent - pLength[chr] })
	names(yOffset) <- chromosomes
	yChrLengths <- sapply(chromosomes, function(chr) { yOffset[chr] + chrLengths[chr] })
	names(yChrLengths) <- chromosomes
	
	yHeight <- max(unlist(yChrLengths))
	xWidth <- yHeight * (4/3)  # 3x4 (y,x) ratio
	numChrs <- length(chromosomes)
	chrWidth <- xWidth/numChrs
	xChrOffset <- sapply(1:numChrs, function(i) { chrWidth*i })
	
	chrCoordinates <- data.frame(name=chromosomes,length=unlist(chrLengths[chromosomes]),centromere=unlist(pLength), yOffset=unlist(yOffset), xOffset = unlist(xChrOffset))
	rownames(chrCoordinates) <- chromosomes
	
	return(list(chrCoordinates =chrCoordinates, dim=c(xWidth, yHeight)))
}
#--------------------------------------------------------------#
getChromosomePositions <- function(chromosomes, chrCoordinates){

	chrPos <- lapply(chromosomes, function(chr){
		offset = chrCoordinates[chr,"yOffset"]
		data.frame(x=chrCoordinates[chr, "xOffset"], p=offset, c=offset+chrCoordinates[chr,"centromere"], q=offset+chrCoordinates[chr,"length"])
	})
	names(chrPos) <- chromosomes
	return(chrPos)

}
#--------------------------------------------------------------#
scaleSamplesToChromosomes <- function(mtx.xy, chrDim){
	
	xOffset <- -1* min(mtx.xy[,"x"])
	yOffset <- -1* min(mtx.xy[,"y"])
	
	mtx.xy[,"x"] = mtx.xy[,"x"] + xOffset; mtx.xy[,"y"] = mtx.xy[,"y"] + yOffset
		# offset mtx so min val is 0,0
	
	mtxDim <- c(max(mtx.xy[,"x"]), max(mtx.xy[,"y"]))
	
	r2Chr <- sum(chrDim*chrDim)
	r2Mtx <- sum(mtxDim*mtxDim)	
	scale <- sqrt(r2Chr/r2Mtx)
		# make diagonal of drawing regions equal
		
	mtx.xy <- mtx.xy * scale

	return(mtx.xy)
}
#--------------------------------------------------------------#
scaleGenesToChromosomes <- function(genePos, chrCoordinates, scaleFactor=1000){
	
	genePos_xy <- lapply(genePos, function(gene){
		x <- chrCoordinates[gene[1], "xOffset"]
		y <- chrCoordinates[gene[1], "yOffset"] + as.numeric(gene[2])/scaleFactor
		c(x,y)
	})

	return(genePos_xy)	
	
}

#--------------------------------------------------------------#
run.batch <- function(scaleFactor=10000){

	hg19_data <- fromJSON(hg19_manifest)

	hg19 <- subset(hg19_data, dataset=="hg19")
  
	chromObj <- hg19[hg19$dataType=="chromosome", "collections"][[1]]
	geneObj  <- hg19[hg19$dataType=="genes",      "collections"][[1]]
	centObj  <- hg19[hg19$dataType=="centromere", "collections"][[1]]
	genesets <- hg19[hg19$dataType=="genesets", "collections"][[1]]
	
	chrLengths <- fromJSON(paste(chromObj$directory, chromObj$file, ".json", sep=""))
	pLength    <- fromJSON(paste(centObj$directory,  centObj$file,  ".json", sep=""))
	genePos    <- fromJSON(paste(geneObj$directory,  geneObj$file,  ".json", sep=""))
	genesets   <- fromJSON(paste(genesets$directory, genesets$file, ".json", sep=""))
	
	chromosomes <- c(seq(1:22), "X", "Y")

	chrSpecs <- getChromosomeOffsets(chromosomes, chrLengths, pLength, scaleFactor=scaleFactor)
	chrPos <- getChromosomePositions(chromosomes, chrSpecs$chrCoordinates)
	save.json(chrPos, hg19_dir, "chromosome_coordinates_scaled")
	
	genePos_scaled <- scaleGenesToChromosomes(genePos, chrSpecs$chrCoordinates, scaleFactor=scaleFactor)
	save.json(genePos_scaled, hg19_dir, paste(genepos_file, "scaled", sep="_"))

	if("mdsScaled" %in% commands) {

		mds_Files<- fromJSON(mds_manifest)
		mds_Files <- subset(mds_Files, dataType=="mds")

		for(mdsFile in mds_Files){
			datasetName <- mdsFile$dataset
			dataType <- mdsFile$dataType
			collections <- mdsFile$collections[[1]]
			
			for(i in 1:nrow(collections)){
				dataObj <- collections[i,]
				process= c(dataObj$process, "scale10k")
				processName <- paste(process, collapse="-")

				mtx <- fromJSON(paste(dataObj$directory, dataObj$file,".json", sep=""))
				mtx <- t(as.data.frame(mtx)); 
				colnames(mtx) <- c("x", "y")
				mtx_scaled <- scaleSamplesToChromosomes(mtx, chrSpecs$dim)
				
				index <- get.new.collection.index(datasetName, dataType)
				file= paste(datasetName, dataType, index, processName, sep="_")
				parent <- list(c(datasetName, dataType, dataObj$id))
				collection <- data.frame(id=index, date=date,directory=mds_scaled_dir, file=file)
				collection$parent <- list(parent)
				collection$process <- list(process)

				add.new.collection(datasetName, dataType, collection)
				os.data.save(mtx_scaled, mds_scaled_dir, file, format="JSON")
			}
		}
	}
	if("geneScaled" %in% commands) {

		datasetName = "genesets"
		datatype = "position"

		genePos_Files<- fromJSON(mds_manifest)
		mds_Files <- subset(mds_Files, dataType=="mds")


		for (genesetName in names(genesets)){	
			genes <- genesets[[genesetName]]
			map_genes <- intersect(genes, names(genePos_scaled))
			genesetPos <- genePos_scaled[map_genes]
			process <- c(genesetName,"scale10k")
			processName <- paste(process, collapse="-")
			
			index <- get.new.collection.index(datasetName, dataType)
			file= paste(datasetName, dataType, index, processName, sep="_")
			parent <- list(c(datasetName, dataType, dataObj$id))
			collection <- data.frame(id=index, date=date,directory=mds_scaled_dir, file=file)
			collection$parent <- list(parent)
			collection$process <- list(process)

			add.new.collection(datasetName, dataType, collection)
			
			file= paste("network_chrPos_", genesetName, sep="")
			os.data.save(genesetPos, hg19_dir, file, format="JSON")
		}	
	}
}

#--------------------------------------------------------------#

	run.batch()
