library(org.Hs.eg.db)
library(jsonlite)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

#--------------------------------------------------------------#

mol_dir<- "../molecular_data/"
hg19_dir <- "../molecular_data/hg19"
ucsc_dir <- "../molecular_data/UCSC"
mds_orig_dir <- "../networks/mds/original"
mds_scaled_dir <- "../networks/mds/scaled"

chromosome_file <- "chromosome_lengths_hg19.json"
centromere_file <- "centromere_position_hg19.json"
genepos_file <- "gene_symbol_min_abs_start_hg19.json"
geneset_file <- "genesets_by_symbol.json"

mol_metadata_file <- "molecular_collections_metadata.txt"

#----------------------------------------------------------------------------------------------------
get.json <- function(directory, file)
{
	if(!grepl("/$", directory)) directory <- paste(directory, "/", sep="")

	data.json <- fromJSON(paste(directory, file, sep=""))
	return(data.json)
}

#----------------------------------------------------------------------------------------------------
save.json <- function(dataObj, directory, file)
{
  if(!dir.exists(directory))
    dir.create(file.path(directory), recursive=TRUE)
  
  if(!grepl("/$", directory)) directory <- paste(directory, "/", sep="")
 
  outFile = paste(directory, file, sep="")
  write(toJSON(dataObj, pretty=TRUE, digits=I(8)), file=paste(outFile,".json", sep = "") )
  
} # saveGraph

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
run.batch <- function(){

	chrLengths <- get.json(hg19_dir,chromosome_file)
	pLength <- get.json(hg19_dir,centromere_file)
	genePos <- get.json(hg19_dir, genepos_file)
	genesets <- get.json(hg19_dir, geneset_file)
	chromosomes <- c(seq(1:22), "X", "Y")

	chrSpecs <- getChromosomeOffsets(chromosomes, chrLengths, pLength)
	chrPos <- getChromosomePositions(chromosomes, chrSpecs$chrCoordinates)
	save.json(chrPos, hg19_dir, "chromosome_coordinates_scaled")
	
	genePos_scaled <- scaleGenesToChromosomes(genePos, chrSpecs$chrCoordinates)
	save.json(genePos_scaled, hg19_dir, paste(genepos_file, "scaled", sep="_"))

	mds_Files<- list.files(mds_orig_dir)
#	mdsFiles <- networkFiles[grep("^mds_", networkFiles)]
	
	for(mdsFile in mds_Files){
		mtx <- get.json(mds_orig_dir, mdsFile)
		mtx <- t(as.data.frame(mtx)); 
		colnames(mtx) <- c("x", "y")
		mtx_scaled <- scaleSamplesToChromosomes(mtx, chrSpecs$dim)
		file= paste(gsub(".json", "", mdsFile), "scaled", sep="_")
		save.json(mtx_scaled, mds_scaled_dir, file)
	}
	
	for (genesetName in names(genesets)){	
		genes <- genesets[[genesetName]]
		map_genes <- intersect(genes, names(genePos_scaled))
		genesetPos <- genePos_scaled[map_genes]
		file= paste("network_chrPos_", genesetName, sep="")
		save.json(genesetPos, hg19_dir, file)
	}	
}

#--------------------------------------------------------------#

	run.batch()
