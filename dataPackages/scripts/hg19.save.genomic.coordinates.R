
# Configuration -----------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

source("common.R")
library(org.Hs.eg.db)

cytoband_url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
chromosomes <- c(seq(1:22), "X", "Y")
date <- as.character(Sys.Date())

#----------------------------------------------------------------------------------------------------
getChromosomeLengths <- function(){
	return( org.Hs.egCHRLENGTHS )
}
#----------------------------------------------------------------------------------------------------
getEntrezSymbolMap <- function(){
	symbolMap <- org.Hs.egSYMBOL
	mapped_genes <- mappedkeys(symbolMap)
	
	return(as.list(symbolMap[mapped_genes]) )
	
}#----------------------------------------------------------------------------------------------------
getGenePositions_Entrez <- function(){

	genePos <- org.Hs.egCHRLOC
	# Get the entrez gene identifiers that are mapped to chromosome locations

	mapped_genes <- mappedkeys(genePos)
	return( as.list(genePos[mapped_genes]) )
}
#----------------------------------------------------------------------------------------------------
getGenePositions_Symbol <- function(){
	EntrezGenePos <- getGenePositions_Entrez()
	  # list of gene start positions, with entrez ID as name.  
	  # Multiple positions reported as named integer vector (name=chromosome)
	  # Locations measured as the number of base pairs from the p (5' end of the sense strand) to q (3' end of the sense strand) arms. 
	  # Locations on antisense strand have a leading "-" sign (e. g. -1234567)
	
	EntrezSymbolMap <- getEntrezSymbolMap()
	
	noSymbol<- which(sapply(EntrezSymbolMap, function(gene){is.na(gene)}))
	stopifnot(length(noSymbol) == 0)
		# all EntrezGene IDs should have symbol
	stopifnot(length(setdiff(names(EntrezGenePos), names(EntrezSymbolMap))) ==0)

	names(EntrezGenePos) <- sapply(names(EntrezGenePos), function(entrezID) { EntrezSymbolMap[[entrezID]] })
		# Set EntrezGenePos names to gene symbol using EntrezSymbolMap
	
	return (EntrezGenePos)
}

#----------------------------------------------------------------------------------------------------
saveChromosome_Coordinates <- function(){
	
	chrLengths <- getChromosomeLengths()
	df<-data.frame(t(chrLengths))
	names(df) <- names(chrLengths)

	result = list(dataset="hg19", type="chromosome", process="length", data=df)
	
	save.collection(mongo, dataset="hg19", dataType="chromosome", source="orgHs",
	                result=list(result),parent=NA, process="length",processName="length")
}	
#----------------------------------------------------------------------------------------------------
saveGene_Coordinates <- function(){

	genePos <- getGenePositions_Symbol()
		# list of all start locations for each gene symbol
			
	genePos_min <- lapply(genePos, function(geneLocations) { 
		geneLocations <- geneLocations[names(geneLocations) %in% chromosomes]
		if(length(geneLocations)==0) return("NULL")

		minLoc <- min(abs(geneLocations))
		chrName <- unique(names(geneLocations)[which(abs(geneLocations) == minLoc)])
		c(chrName, as.integer(minLoc))
	})
	# list with chr name & min (absolute value) chr location for each gene symbol
	
	notMapped <- which(genePos_min == "NULL")
	genePos_min[notMapped] <- NULL
		# removes gene positions that map to chromosomes outside our list (ie 1-22, X, Y)

	process = list("position", "min", "abs", "start")
	processName = paste(unlist(process), collapse="-")
	result = list(dataset="hg19", type="genes", process=process, data=genePos_min)
	
	save.collection(mongo, dataset="hg19", dataType="genes", source="orgHs",
	                result=list(result),parent=NA, process=process,processName=processName)
	
}
#----------------------------------------------------------------------------------------------------
saveCentromere_Coordinates <- function(cytoband_url){
	
	temp <- tempfile()
	download.file(cytoband_url,temp)
	cytobands <- read.delim(gzfile(temp), header=F)
	unlink(temp)

	colnames(cytobands) <- c("chr", "start", "end", "cytoband", "other")
	cytobands$chr <- sub("chr", "", cytobands$chr)
	
	centromere <- sapply(chromosomes, function(chrName){
		chr_cyto <- subset(cytobands, chr == chrName)
		chr_cyto_p <- subset(chr_cyto, grepl("^p", cytoband))
		chr_p_end <- max(chr_cyto_p$end)
			# == chr_q_start 
		
		chr_p_end
	})
	
	df<-data.frame(t(centromere))
	names(df) <- names(centromere)
	
	result = list(dataset="hg19", type="centromere", process="position", data=df)
	
	save.collection(mongo, dataset="hg19", dataType="centromere", source="orgHs",
	                result=list(result),parent=NA, process="position",processName="position")
	
}

#----------------------------------------------------------------------------------------------------
## must first initialize server (through shell >mongod)
mongo <- connect.to.mongo()

	saveChromosome_Coordinates()
	saveGene_Coordinates()
	saveCentromere_Coordinates(cytoband_url)	

close.mongo(mongo)



	