library(org.Hs.eg.db)
library(jsonlite)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

commands <- c("chromosomeScaled", "mdsScaled", "geneScaled")
scaleFactor <- 10000
date <- Sys.Date()

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

#---------------------------------------------------------
get.new.collection.index <- function(Manifest, datasetName, dataTypeName){
  
  if(nrow(Manifest) == 0) return(1)
  
  dataObj <- subset(Manifest, dataset == datasetName & dataType == dataTypeName)
  if(nrow(dataObj) == 0) return(1)
  
  return(nrow(dataObj$collections[[1]]) +1)
}
#---------------------------------------------------------
add.new.collection <- function(Manifest, datasetName, dataTypeName, collection){
  
  if(nrow(Manifest) == 0){	
    newCollection <- data.frame(dataset=datasetName, dataType=dataTypeName)
    newCollection$collections <- list(collection)
    Manifest <- newCollection
    return(Manifest)
  }
  
  dataObj <- subset(Manifest, dataset == datasetName & dataType == dataTypeName)
  if(nrow(dataObj) == 1){
    newCollection <- list(rbind(dataObj$collections[[1]],collection))
    Manifest[Manifest$dataset==datasetName & Manifest$dataType ==dataTypeName,"collections"] <- list(newCollection)
    return(Manifest)
  }
  if(nrow(dataObj) == 0){	
    newCollection <- data.frame(dataset=datasetName, dataType=dataTypeName)
    newCollection$collections <- list(collection)
    Manifest <- rbind(Manifest, newCollection)
    return(Manifest)
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
	mtx.xy <- round(mtx.xy)

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
get.chromosome.dimensions <- function(manifest, scaleFactor=1000){
  
  chrPosObj <- subset(manifest, dataType=="chromosome")$collections[[1]]
  chrPosScaledObj <- chrPosObj[sapply(chrPosObj$process, function(proc){all(proc == c("scaled", scaleFactor))}),]

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
	
	## create collection: Chromosome Positions
	process <- c("scaled", scaleFactor); processName <- paste(process, collapse="-")
	## Chromosome Positions
	index_chr <- get.new.collection.index(manifest, chromObj$dataset, chromObj$dataType)
	outfile_chr <- paste(chromObj$dataset, chromObj$dataType, index_chr, processName, sep="_")
	parent <- list(c(chromObj$dataset, chromObj$dataType, chrLenObj$id), c(centObj$dataset, centObj$dataType, centPosObj$id))
	collection_chr <- data.frame(id=index_chr, date=date, directory=hg19_dir, file=outfile_chr)
	collection_chr$parent <- list(parent)
	collection_chr$process <- list(process)
	## Gene Positions
	index_gene <- get.new.collection.index(manifest, geneObj$dataset, geneObj$dataType)
	outfile_gene <- paste(geneObj$dataset, geneObj$dataType, index_gene, processName, sep="_")
	parent <- list(c(geneObj$dataset, geneObj$dataType, genePosObj$id), c(chrLenObj$dataset, chrLenObj$dataType, chrLenObj$id))
	collection_gene <- data.frame(id=index_gene, date=date, directory=hg19_dir, file=outfile_gene)
	collection_gene$parent <- list(parent)
	collection_gene$process <- list(process)
	
	#save and update manifest
	manifest <- add.new.collection(manifest,chromObj$dataset, chromObj$dataType, collection_chr)
	os.data.save(chrPos, hg19_dir, outfile_chr, format="JSON")
	manifest <- add.new.collection(manifest, geneObj$dataset, geneObj$dataType, collection_gene)
	os.data.save(genePos_scaled, hg19_dir, outfile_gene, format="JSON")

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
			
			orig_mds <- collections[sapply(collections$process, function(proc){"mds" %in% proc }),]
			
			for(j in 1:nrow(orig_mds)){
				dataObj <- orig_mds[j,]

				ptList <- fromJSON(paste(dataObj$directory, dataObj$file,".json", sep=""))
				mtx <- t(sapply(ptList, function(id){ return(c(id$x, id$y))}))
				colnames(mtx) <- c("x", "y")
				mtx_scaled <- scaleSamplesToChromosomes(mtx, chrDim)
				colnames(mtx_scaled) <- c("x", "y")
				
				index <- get.new.collection.index(manifest, datasetName, dataType)
				process= c("scaled", scaleFactor);				processName <- paste(process, collapse="-")
				file= paste(datasetName, dataType, index, processName, sep="_")
				parent <- list(c(datasetName, dataType, dataObj$id))
				
				collection <- data.frame(id=index, date=date,directory=mds_scaled_dir, file=file)
				collection$parent <- list(parent)
				collection$process <- list(process)

				manifest <- add.new.collection(manifest, datasetName, dataType, collection)
				os.data.save(mtx_scaled, mds_scaled_dir, file, format="JSON")
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
			
			process <- c("scaled", scaleFactor, genesetName); processName <- paste(process, collapse="-")
			
			index <- get.new.collection.index(manifest,geneObj$dataset, geneObj$dataType)
			file= paste(geneObj$dataset, geneObj$dataType, index, processName, sep="_")
			parent <- list(c(geneObj$dataset, geneObj$dataType, genePosObj$id),c(genesetObj$dataset, genesetObj$dataType, genesetSymObj$id))
			collection <- data.frame(id=index, date=date,directory=hg19_dir, file=file)
			collection$parent <- list(parent)
			collection$process <- list(process)

			manifest <- add.new.collection(manifest, geneObj$dataset, geneObj$dataType, collection)
			os.data.save(genesetPos, hg19_dir, file, format="JSON")
		}	
		
		return(manifest)
}


#--------------------------------------------------------------#
if("chromosomeScaled" %in% commands) 
  Manifest_hg19<- run.scale.chr.genes(Manifest_hg19,scaleFactor)
if("mdsScaled" %in% commands){
   chrDim <- get.chromosome.dimensions(Manifest_hg19, scaleFactor) 
   Manifest_mds <- run.batch.mds(Manifest_mds,chrDim, scaleFactor)
}
if("geneScaled" %in% commands) 
  Manifest_hg19<- run.batch.genesets(Manifest_hg19, scaleFactor)

os.data.save(Manifest_hg19,"../manifests/", "os.hg19.scaled.manifest", format="JSON")
os.data.save(Manifest_mds,"../manifests/", "os.mds.network.manifest", format="JSON")
