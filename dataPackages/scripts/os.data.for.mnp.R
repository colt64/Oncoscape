###
#       
###

# Library Imports ---------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

source("common.R")

# Configuration -----------------------------------------------------------

date <- as.character(Sys.Date())
commands <- c("chromosome", "categories", "layouts", "edges")
#commands <- c("chromosome")
scaleFactor = 100000

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args


#os.mds.manifest   <- fromJSON("../manifests/os.mds.scaled.manifest.json")
os.ptLayout.manifest   <- fromJSON("../manifests/os.pca.scaled.manifest.json")
#os.ptlayout.manifest <- rbind(os.mds.manifest, os.pca.manifest)

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
	 

	for(i in 1:length(datasetNames)){
	  	datasetName = datasetNames[i]
  		mdsCollections <- subset(manifest, dataset==datasetName & dataType == "mds")$collections[[1]]
  		pcaCollections <- subset(manifest, dataset==datasetName & dataType == "pcaScores")$collections[[1]]
  		ptLayoutSet <- rbind(mdsCollections, pcaCollections)
  		ptLayoutSet <- subset.collections(ptLayoutSet, process=list(calculation="scaled", input=scaleFactor))

  		for(j in 1:nrow(ptLayoutSet)){
  		  ptLayoutObj <- ptLayoutSet[j,]
  		  
  		  parentObj <- get.parent.collection(manifest, ptLayoutObj$parent)
  		  name <- paste(c(unlist(ptLayoutObj$source), parentObj$process[[1]]$calculation, unlist(parentObj$process[[1]]$input), parentObj$process[[1]]$geneset), collapse="-")
  		  
  	    ptLayout.data <- fromJSON(paste(ptLayoutObj$directory, ptLayoutObj$file,".json", sep=""))
	      MnP.ptLayouts <- c(MnP.ptLayouts, list(data.frame(dataset=datasetName, datatype=datatypeName, name=name)))
	      MnP.ptLayouts[[length(MnP.ptLayouts)]]$data=list(ptLayout.data)
	    }
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
    edgeCollections <- subset(manifest, dataset==datasetName & dataType=="edges")$collections[[1]]
      ptCollections <- subset(manifest, dataset==datasetName & dataType=="ptDegree")$collections[[1]]
    geneCollections <- subset(manifest, dataset==datasetName & dataType=="geneDegree")$collections[[1]]
    
    # get Edge sets for each Gene set
    for(j in 1:length(genesets)){			
      genesetName <- names(genesets)[j]
#      edgeObj <- subset.collections(edgeCollections, process <- list(geneset=genesetName, edgeType=c("cnv", "mut01")))
      edgeObj <- subset(edgeCollections, 
                            process$geneset ==genesetName &
                            all(c("cnv", "mut01") %in% unlist(process$edgeType)))
      edgeObj <- edgeObj[unlist(lapply(edgeObj$process$edgeType, function(eType) all(eType== c("cnv", "mut01")))),]
      edgeObj <- edgeObj[unlist(lapply(edgeObj$parent, function(parent) all( parent[,3] == 2))),]
      ### HACK!! currently hardcoding in that edge types come from UCSC data - which is defined 2nd in manifest
      #      parentObj <- get.parent.collection(manifest, edgeObj$parent)
      ptObj   <- subset(ptCollections,   id==edgeObj$id )  
      geneObj <- subset(geneCollections, id==edgeObj$id )  
      
      edge.data <- list(
          edges = fromJSON(paste(edgeObj$directory, edgeObj$file,".json", sep="")),
          ptDegree = fromJSON(paste(ptObj$directory, ptObj$file,".json", sep="")),
          geneDegree = fromJSON(paste(geneObj$directory, geneObj$file,".json", sep=""))
      )          
      MnP.edges[[j]] <- data.frame(dataset=datasetName, datatype=datatypeName, name=genesetName)
      MnP.edges[[j]]$data=list(edge.data)
    }
  }
  os.data.save(MnP.edges, output.dir, "os.MnP.edges.data", format="JSON")
}

##----------------------------
os.copy.chromosome.layout <- function(manifest, scaleFactor, outputDir= "./"){
  
  datasetName="hg19"
  
  chr_collections <- subset(manifest, dataset==datasetName & dataType=="chromosome")$collections[[1]]
  process <- list(calculation="scaled", input= scaleFactor)
  chr_collections <- subset.collections(chr_collections, process)
 
  cc <- apply(chr_collections,1, function(coll){
    data <- fromJSON(paste(coll$directory, coll$file,".json", sep=""))
    output <- data.frame(dataset=datasetName, datatype="chromosome", 
                         name=paste(coll$source, coll$process$calculation, coll$process$scaleFactor, sep="-"))
    output$data <- list(data)
    
    os.data.save(output, outputDir, coll$file, format= "JSON")
  });
  
  gene_collections <- subset(manifest, dataset==datasetName & dataType=="genes")$collections[[1]]
  process <- list(calculation="scaled", input= scaleFactor)
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
 os.save.ptLayouts(os.ptLayout.manifest, datasetNames = c("brca"))
if("edges" %in% commands) 
 os.save.network.edges(os.edges.manifest, datasetNames = c("brca"))
if("categories" %in% commands) 
  os.copy.file(inputDir ="../data/categories/", filename=os.categories.data.filename, outputDir="../data/tools/MarkerPatient/")
if("chromosome" %in% commands) 
  os.copy.chromosome.layout(os.hg19.scaled.manifest,scaleFactor, outputDir="../data/tools/MarkerPatient/")
