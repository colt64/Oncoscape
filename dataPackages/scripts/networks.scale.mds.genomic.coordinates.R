library(org.Hs.eg.db)
library(jsonlite)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

commands <- c("chromosomeScaled", "mdsScaled", "pcaScaled", "geneScaled")
#commands <- c("chromosomeScaled", "geneScaled")
#commands <- c("chromosomeScaled", "pcaScaled", "mdsScaled")
scaleFactor <- 100000 
date <- as.character(Sys.Date())

source("common.R")

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

#--------------------------------------------------------------#

mds_manifest_file <- "../manifests/os.mds.manifest.json"
hg19_manifest_file <- "../manifests/os.hg19.manifest.json"

mds_scaled_dir <- "../data/molecular/mds/scaled/"
hg19_dir <- "../data/molecular/hg19/"

Manifest_hg19 <- fromJSON(hg19_manifest_file)
Manifest_mds  <- fromJSON(mds_manifest_file)


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
scaleSamplesToChromosomes <- function(mtx, chrDim, dim.names=c("x", "y", "z")){
  
  mtx <- apply(mtx, 2, function(col){ -1* min(col) + col})
    # offset mtx so min val is 0,0
  mtx.max <- apply(mtx, 2, max)
  
  r2Chr <- sum(chrDim*chrDim)
  r2Mtx <- sum(mtx.max*mtx.max)	
  scale <- sqrt(r2Chr/r2Mtx)
  # make diagonal of drawing regions equal
  
  mtx <- mtx * scale
  mtx <- round(mtx)
  
  list.coord <- lapply(rownames(mtx), function(name){
    df <- data.frame(mtx[name,dim.names])
    colnames(df) <- dim.names
    df
  })
  names(list.coord) <- rownames(mtx)
  return(list.coord)
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
	mtx.xy <- round(mtx.xy)

	list.xy <- lapply(rownames(mtx.xy), function(name){
	  x<- mtx.xy[name,"x"]
	  y<- -1*  mtx.xy[name,"y"]
	  data.frame(x=x,y=y)
	})
	names(list.xy) <- rownames(mtx.xy)
	return(list.xy)
}
#--------------------------------------------------------------#
scaleGenesToChromosomes <- function(genePos, chrCoordinates, scaleFactor=1000){
	
	genePos_xy <- lapply(genePos, function(gene){
		x <- chrCoordinates[gene[1], "xOffset"]
		y <- chrCoordinates[gene[1], "yOffset"] + as.numeric(gene[2])/scaleFactor
		data.frame(x=round(x),y=round(y))
	})

	return(genePos_xy)	
	
}
#--------------------------------------------------------------#
get.chromosome.dimensions <- function(manifest, scaleFactor=1000){
  
  chrPosObj <- subset(manifest, dataType=="chromosome")$collections[[1]]
  chrPosScaledObj <- chrPosObj[which(sapply(chrPosObj$process, function(proc){proc["calculation"] == "scaled" & proc["input"] == scaleFactor})),]

  chrCoord <- fromJSON(paste(chrPosScaledObj$directory, chrPosScaledObj$file, ".json", sep=""))
  chrPos_xy <-t(sapply(chrCoord, function(chr){ return(c(chr$x, chr$q))}))
  chrDim <- c(max(chrPos_xy[,1]), max(chrPos_xy[,2]))

  return(chrDim)
}
#--------------------------------------------------------------#
#creates 2 files: scaled positions of chromosomes and all genes with offsets to align centromeres
run.scale.chr.genes <- function(manifest, scaleFactor=10000){

  # define data objects
	chromObj <- subset(manifest, dataType=="chromosome")
	geneObj  <- subset(manifest, dataType=="genes")
	centObj  <- subset(manifest, dataType=="centromere")
	
	chrLenObj <- subset(chromObj$collections[[1]], process=="length")
	centPosObj <- subset(centObj$collections[[1]], process=="position")
	genePosObj <- subset(geneObj$collections[[1]], all(c("position", "min", "abs", "start") %in% unlist(process)))
	  
	chrLengths <- fromJSON(paste(chrLenObj$directory,   chrLenObj$file,   ".json", sep=""))
	pLength    <- fromJSON(paste(centPosObj$directory,  centPosObj$file,  ".json", sep=""))
	genePos    <- fromJSON(paste(genePosObj$directory,  genePosObj$file,  ".json", sep=""))

	chromosomes <- c(seq(1:22), "X", "Y")

	## calculate Chromosome & Gene positions with scaling
	chrSpecs <- getChromosomeOffsets(chromosomes, chrLengths, pLength, scaleFactor=scaleFactor)
	chrPos <- getChromosomePositions(chromosomes, chrSpecs$chrCoordinates)
	genePos_scaled <- scaleGenesToChromosomes(genePos, chrSpecs$chrCoordinates, scaleFactor=scaleFactor)
	
	## create collections
	process <- data.frame(calculation="scaled", input=scaleFactor); processName <- paste(process, collapse="-")
	## Chr Positions
	parent <- list(c(chromObj$dataset, chromObj$dataType, chrLenObj$id), c(centObj$dataset, centObj$dataType, centPosObj$id))
	manifest <- save.collection(Manifest=manifest, dataset=chromObj$dataset, dataType=chromObj$dataType, result=chrPos,
	                            parent=list(parent), process=list(process),processName=processName, outputDirectory=hg19_dir)
	## Gene Positions
	parent <- list(c(geneObj$dataset, geneObj$dataType, genePosObj$id), c(chrLenObj$dataset, chrLenObj$dataType, chrLenObj$id))
	manifest <- save.collection(Manifest=manifest, dataset=geneObj$dataset, dataType=geneObj$dataType, result=genePos_scaled,
	                            parent=list(parent), process=list(process),processName=processName, outputDirectory=hg19_dir)

	return(manifest)
}	
#--------------------------------------------------------------#
run.batch.pca <- function(manifest,chrDim, scaleFactor=10000){
  
  pca_Files <- subset(manifest, dataType=="pcaScores")
  
  for(i in 1:nrow(pca_Files)){
    pcaFile  <- pca_Files[i,]
    datasetName <- pcaFile$dataset
    dataType <- pcaFile$dataType
    collections <- pcaFile$collections[[1]]
    
    orig_pca <- collections[sapply(collections$process, function(proc){"prcomp" %in% proc$calculation }),]
    
    for(j in 1:nrow(orig_pca)){
      dataObj <- orig_pca[j,]
      
      scores <- fromJSON(paste(dataObj$directory, dataObj$file,".json", sep=""))

      #3d
      mtx.3d <- scores$data[,1:3]
      colnames(mtx.3d) <- c("x", "y", "z")
      rownames(mtx.3d) <- scores$rows
      list3d_scaled <- scaleSamplesToChromosomes(mtx, chrDim)
      
      mtx <- scores$data[,1:2]
      colnames(mtx) <- c("x", "y")
      rownames(mtx) <- scores$rows
      list_scaled <- scaleSamplesToChromosomes(mtx, chrDim, dim.names=c("x", "y"))

      parent <- list(c(datasetName, dataType, dataObj$id))
      geneset <- dataObj$process[[1]]$geneset[[1]]
      if(is.null(geneset)) geneset <- NA
      process= data.frame(calculation="scaled", input=scaleFactor, geneset=geneset, dim="2d");
      processName <- paste(process, collapse="-")
      manifest <- save.collection(Manifest=manifest, dataset=datasetName, dataType=dataType, result=list_scaled,
                                  parent=parent, process=list(process),processName=processName, outputDirectory=mds_scaled_dir)

      process= data.frame(calculation="scaled", input=scaleFactor, geneset=geneset, dim="3d");
      processName <- paste(process, collapse="-")
      manifest <- save.collection(Manifest=manifest, dataset=datasetName, dataType=dataType, result=list3d_scaled,
                                  parent=parent, process=list(process),processName=processName, outputDirectory=mds_scaled_dir)
      
      
      }
  }
  
  return(manifest)
}

#--------------------------------------------------------------#
run.batch.mds <- function(manifest,chrDim, scaleFactor=10000){

		mds_Files <- subset(manifest, dataType=="mds")
		
		
		for(i in 1:nrow(mds_Files)){
		  mdsFile  <- mds_Files[i,]
			datasetName <- mdsFile$dataset
			dataType <- mdsFile$dataType
			collections <- mdsFile$collections[[1]]
			
			orig_mds <- collections[sapply(collections$process, function(proc){"mds" %in% proc$calculation }),]
			
			for(j in 1:nrow(orig_mds)){
				dataObj <- orig_mds[j,]

				ptList <- fromJSON(paste(dataObj$directory, dataObj$file,".json", sep=""))
				mtx <- t(sapply(ptList, function(id){ return(c(id$x, id$y))}))
				colnames(mtx) <- c("x", "y")
				mtx[,"y"] <- -1 * mtx[,"y"]
				list_scaled <- scaleSamplesToChromosomes(mtx, chrDim)

				process= data.frame(calculation="scaled", scaleFactor=scaleFactor, geneset=dataObj$process[[1]]$geneset[[1]], dim="2d");				
				processName <- paste(process, collapse="-")
				parent <- list(c(datasetName, dataType, dataObj$id))
				manifest <- save.collection(Manifest=manifest, dataset=datasetName, dataType=dataType, result=list_scaled,
				                            parent=parent, process=list(process),processName=processName, outputDirectory=mds_scaled_dir)
			}
		}
		
		return(manifest)
}

#--------------------------------------------------------------#
run.batch.genesets <- function(manifest, scaleFactor=10000){
 
		geneObj<- subset(manifest, dataType=="genes")
		geneColl <- geneObj$collections[[1]]
		
		genePosObj <- geneColl[sapply(geneColl$process, function(proc){all(proc == c("scaled", scaleFactor))}),]
		genePos_scaled <- fromJSON(paste(genePosObj$directory, genePosObj$file,".json", sep=""))
		
		genesetObj <-  subset(manifest, dataType=="genesets")
		genesetColl <- genesetObj$collections[[1]]
		
		genesetSymObj <- subset(genesetColl, process == "hgnc")
		genesets <- fromJSON(paste(genesetSymObj$directory, genesetSymObj$file, ".json", sep=""))

		for (genesetName in names(genesets)){	
			genes <- genesets[[genesetName]]
			map_genes <- intersect(genes, names(genePos_scaled))
			genesetPos <- genePos_scaled[map_genes]
			
			process <- data.frame(calculation="scaled", scaleFactor=scaleFactor, geneset=genesetName); 
			processName <- paste(process, collapse="-")
			parent <- list(c(geneObj$dataset, geneObj$dataType, genePosObj$id),c(genesetObj$dataset, genesetObj$dataType, genesetSymObj$id))
			manifest <- save.collection(Manifest=manifest, dataset=geneObj$dataset, dataType=geneObj$dataType, result=genesetPos,
			                            parent=list(parent), process=list(process),processName=processName, outputDirectory=hg19_dir)
		}	
		
		return(manifest)
}


#--------------------------------------------------------------#
if("chromosomeScaled" %in% commands) {
  Manifest_hg19<- run.scale.chr.genes(Manifest_hg19,scaleFactor)
  os.data.save(Manifest_hg19,"../manifests/", "os.hg19.scaled.manifest", format="JSON")
}
if("mdsScaled" %in% commands){
   chrDim <- get.chromosome.dimensions(Manifest_hg19, scaleFactor) 
   Manifest_mds <- run.batch.mds(Manifest_mds,chrDim, scaleFactor)
   os.data.save(Manifest_mds,"../manifests/", "os.mds.scaled.manifest", format="JSON")
}
if("pcaScaled" %in% commands){
  chrDim <- get.chromosome.dimensions(Manifest_hg19, scaleFactor) 
  Manifest_pca <- run.batch.pca(Manifest_mds,chrDim, scaleFactor)
  os.data.save(Manifest_pca,"../manifests/", "os.pca.scaled.manifest", format="JSON")
}

if("geneScaled" %in% commands) {
  Manifest_hg19<- run.batch.genesets(Manifest_hg19, scaleFactor)
  os.data.save(Manifest_hg19,"../manifests/", "os.hg19.scaled.manifest", format="JSON")
}

