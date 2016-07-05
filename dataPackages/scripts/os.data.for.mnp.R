###
#       
###

# Library Imports ---------------------------------------------------------
library(jsonlite)
source("common.R")

# Configuration -----------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

date <- as.character(Sys.Date())
#commands <- c("categories", "layouts", "edges")
commands <- "chromosome"

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args


os.mds.manifest   <- fromJSON("../manifests/os.mds.scaled.manifest.json")
os.edges.manifest   <- fromJSON("../manifests/os.edges.manifest.json")
os.hg19.scaled.manifest   <- fromJSON("../manifests/os.hg19.scaled.manifest.json")
os.categories.data.filename <-"os.categories.color.data"

geneset_file <- "../data/molecular/hg19/hg19_genesets_1_hgnc.json"
genesets <- fromJSON(geneset_file)

output.dir <- "../data/tools/markerPatient/"

#----------------------------------------------------------------------------------------------------
os.save.ptLayouts <- function(manifest, datasetNames){

	datatypeName= "ptLayout"
	MnP.ptLayouts <- list()
	# 
	## Patient Layout by MDS CNV/SNV OV
	name = "MDS-CNV;SNV-OV"
	
	for(i in 1:length(datasetNames)){
	  	datasetName = datasetNames[i]
  		mdsCollections <- subset(manifest, dataset==datasetName & dataType=="mds")$collections[[1]]
	  	matchColl <- apply(mdsCollections,1, function(coll){ coll$process$calculation == "scaled" & coll$process$geneset=="oncoVogel274"})
  		ptLayoutObj <- mdsCollections[matchColl,]
  		
	    ptLayout.data <- fromJSON(paste(ptLayoutObj$directory, ptLayoutObj$file,".json", sep=""))
	    MnP.ptLayouts[[i]] <- data.frame(dataset=datasetName, datatype=datatypeName, name=name)
	    MnP.ptLayouts[[i]]$data=list(ptLayout.data)
  }
	os.data.save(MnP.ptLayouts, output.dir, "os.MnP.ptLayouts.data", format="JSON")
		
}

#----------------------------------------------------------------------------------------------------
os.save.network.edges <- function(manifest, datasetNames){
  
  datatypeName= "networkEdges"
  MnP.edges <- list()
  # 
  ## Patient Layout by MDS CNV/SNV OV
  
  for(i in 1:length(datasetNames)){
    datasetName = datasetNames[i]
    edgeCollections <- subset(manifest, dataset==datasetName & dataType=="network")$collections[[1]]
    
    # get Edge sets for each Gene set
    for(j in 1:length(genesets)){			
      genesetName <- names(genesets)[j]
      edgeObj <- subset(edgeCollections, 
                            process$geneset ==genesetName &
                            all(c("cnv", "mut01") %in% unlist(process$edgeType)))
      edgeObj <- edgeObj[unlist(lapply(edgeObj$process$edgeType, function(eType) all(eType== c("cnv", "mut01")))),]
    
      edge.data <- list(
          edges = fromJSON(paste(edgeObj$directory, edgeObj$file[[1]][1],".json", sep="")),
          ptDegree = fromJSON(paste(edgeObj$directory, edgeObj$file[[1]][2],".json", sep="")),
          geneDegree = fromJSON(paste(edgeObj$directory, edgeObj$file[[1]][3],".json", sep=""))
      )          
      MnP.edges[[j]] <- data.frame(dataset=datasetName, datatype=datatypeName, name=genesetName)
      MnP.edges[[j]]$data=list(edge.data)
    }
  }
  os.data.save(MnP.edges, output.dir, "os.MnP.edges.data", format="JSON")
}

##----------------------------
os.copy.chromosome.layout <- function(manifest, outputDir= "./"){
  
  datasetName="hg19"
  scaleFactor = 10000
  
  chr_collections <- subset(manifest, dataset==datasetName & dataType=="chromosome")$collections[[1]]
  process <- data.frame(field=c("calculation", "scaleFactor"), value=c("scaled", scaleFactor))
  chr_collections <- subset.collections(chr_collections, process)
 
  cc <- apply(chr_collections,1, function(coll){
    data <- fromJSON(paste(coll$directory, coll$file,".json", sep=""))
    output <- data.frame(dataset=datasetName, datatype="chromosome", 
                         name=paste(coll$process$calculation, coll$process$scaleFactor, sep="-"))
    output$data <- list(data)
    
    os.data.save(output, outputDir, coll$file, format= "JSON")
  });
  
  gene_collections <- subset(manifest, dataset==datasetName & dataType=="genes")$collections[[1]]
  process <- data.frame(field=c("calculation", "scaleFactor"), value=c("scaled", scaleFactor))
  gene_collections <- subset.collections(gene_collections, process)
  
  gc <- apply(gene_collections,1, function(coll){
    data <- fromJSON(paste(coll$directory, coll$file,".json", sep=""))
    geneset <- coll$process$geneset; if(is.null(geneset)) geneset <- ""
    output <- data.frame(dataset=datasetName, datatype="genes", name=geneset)
    output$data <- list(data)
      
    os.data.save(output, outputDir,coll$file,  format= "JSON")
  });
  
}

##----------------------------
if("layouts" %in% commands) 
 os.save.ptLayouts(os.mds.manifest, datasetNames = c("brca"))
if("edges" %in% commands) 
 os.save.network.edges(os.edges.manifest, datasetNames = c("brca"))
if("categories" %in% commands) 
  os.copy.file(inputDir ="../data/categories/", filename=os.categories.data.filename, outputDir="../data/tools/MarkerPatient/")
if("chromosome" %in% commands) 
  os.copy.chromosome.layout(os.hg19.scaled.manifest, outputDir="../data/tools/MarkerPatient/")
