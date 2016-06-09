library(org.Hs.eg.db)
library(jsonlite)

#--------------------------------- make plot data -----------------------------#
directory <- "../molecular_data/hg19"
chr_file <- "chromosome_lengths_hg19"
gene_file <- "gene_symbol_min_abs_start_hg19"
cent_file <- "centromere_position_hg19"
cytoband_url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"

chromosomes <- c(seq(1:22), "X", "Y")

#----------------------------------------------------------------------------------------------------
save.json <- function(dataObj, directory, file)
{
  if(!dir.exists(directory))
    dir.create(file.path(directory), recursive=TRUE)
  
  if(!grepl("/$", directory)) directory <- paste(directory, "/", sep="")
 
  outFile = paste(directory, file, sep="")
  write(toJSON(dataObj, pretty=TRUE, digits=I(8)), file=paste(outFile,".json", sep = "") )
  
} # saveGraph


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
saveChromosome_Coordinates <- function(out_file){
	
	chrLengths <- getChromosomeLengths()
	df<-data.frame(t(chrLengths))
	names(df) <- names(chrLengths)
	save.json(df, directory, file=out_file)
}	
#----------------------------------------------------------------------------------------------------
saveGene_Coordinates <- function(out_file){

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

	save.json(genePos_min, directory, out_file)

}
#----------------------------------------------------------------------------------------------------
saveCentromere_Coordinates <- function(cytoband_url, out_file){
	
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
	
	save.json(df, directory, out_file)
}
#----------------------------------------------------------------------------------------------------

	saveChromosome_Coordinates(chr_file)
	saveGene_Coordinates(gene_file)
	saveCentromere_Coordinates(cytoband_url, cent_file)	
	
	