library(org.Hs.eg.db)
library(jsonlite)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

commands <- c("mds", "edges")
args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

#--------------------------------- Configuration -----------------------------#

molecular_manifest <- "../manifests/os.ucsc.molecular.manifest.json"

network_output_directory <- "../data/molecular/edges"
mds_output_directory <- "../data/molecular/mds"

geneset_file <- "../data/molecular/hg19/genesets_by_symbol.json"
genesets <- fromJSON(geneset_file)

Manifest <- data.frame()
date <- Sys.Date()

#----------------------------------------------------------------------------------------------------
getGeneSet <- function(geneset_name){
	stopifnot(geneset_name %in% names(genesets))
	return(genesets[[geneset_name]])
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
  write(toJSON(dataObj, pretty=TRUE, digits=I(8)), file=paste(outFile,".json", sep = "") )
  
  # Return DataFrame For Chaining
  return(df)
}

#----------------------------------------------------------------------------------------------------
calcSimilarity <- function(indicatorMatrix) {
	similarity=NULL
	similarity <- apply(indicatorMatrix,2, function(ptCol){
		ptCol %*% indicatorMatrix
	})
	diag(similarity) <- 1
	rownames(similarity) <- colnames(indicatorMatrix)
	colnames(similarity) <- colnames(indicatorMatrix)
	return(similarity)
}

#----------------------------------------------------------------------------------------------------
calculateSampleSimilarityMatrix <- function (mut, cn, samples=NA, genes=NA, copyNumberValues=c(-2, 2), threshold=NA) {

    if(!all(is.na(samples))){
        mut <- mut[,intersect(colnames(mut), samples)]
        cn  <- cn[ ,intersect(colnames(cn) , samples)]
    }

    if(!all(is.na(genes))){
        mut <- mut[intersect(rownames(mut), genes),]
        cn  <- cn[ intersect(rownames(cn) , genes),]
    }

	#remove any genes with NA in mutation
	tmp <- apply(mut, 1, function(x) any(is.na(x)))
	if(length(which(tmp))>0)
		mut <- mut[-which(tmp),]

	stopifnot(all(sort(unique(as.integer(mut))) == c(0,1)))
     
    cn[!cn %in% copyNumberValues] <- 0

	similaritySNV <- calcSimilarity(as.matrix(mut))
	similarityCNV <- calcSimilarity(as.matrix(cn))

	sharedSnvCnv <- intersect(rownames(similaritySNV), rownames(similarityCNV))
	simSNV <- similaritySNV[sharedSnvCnv, sharedSnvCnv]
	simCNV <- similarityCNV[sharedSnvCnv, sharedSnvCnv]

	SNV.CNV <- ((simSNV)/sum(simSNV)) + 
			   ((simCNV)/sum(simCNV))

	D <- as.dist(max(SNV.CNV) - SNV.CNV)
	tbl.pos <- cmdscale(D, k=2) #MDS.SNV.CNV
	colnames(tbl.pos) <- c("x", "y")
	tbl.pos <- as.data.frame(tbl.pos)


#	 ptIDs <- canonicalizePatientIDs(obj@pkg, rownames(tbl.pos))
#	 tbl.pos <- tbl.pos[!duplicated(ptIDs),]
#     rownames(tbl.pos) <- ptIDs[!duplicated(ptIDs)]
     
     return(tbl.pos)
}

#----------------------------------------------------------------------------------------------------
run.batch.patient_similarity <- function(manifest_file, geneset_name=NA, output_directory="./"){

  # Load Input File 
  datasets <- fromJSON(manifest_file)
  dataType <- "mds"
  process <- c("mds", "cnv;mut01")

  gistic.scores <-c(-2,-1,1, 2)
  goi = getGeneSet(geneset_name)
  
  # Loop for each dataset
  datasetNames <- unique(datasets$dataset)
  for (datasetName in datasetNames){
      molTables <- subset(datasets, dataset==datasetName && dataType %in% c("cnv", "mut01"))
    
      cnvTables <- subset(molTables, dataType == "cnv")
      mutTables <- subset(molTables, dataType == "mut01")
    
      if(nrow(cnvTables)==0 & nrow(mutTables) ==0) next;
      cat(dataset, "\n")		  
    
			## ----- Configuration ------
			regex = ".01$"; threshold = NA;
			if(datasetName == "laml"){        regex = "-03$|-09$";
			} else if(datasetName == "luad"){ regex = "TCGA-(17)^-\\d{4}-01$" }

			if(datasetName == "brca" | datasetName == "brain")  threshold = -1e-04

		  PtLayouts[[datasetName]] <- c()
			cnvTables <- subset(collections, molecular_type=="cnv")
			mutTables <- subset(collections, molecular_type=="mutation01")

			if(nrow(cnvTables)==0 | nrow(mutTables) ==0) next;
			
			## ----- MDS on All Combinations of CNV and MUT Tables ------
			for(i in 1:nrow(cnvTables)){
				cnv.json <- fromJSON(paste(cnvTables[i,"directory"], cnvTables[i,"file"], sep="/")) 
				mtx.cnv <- cnv.json$data[[1]]
				rownames(mtx.cnv) <- cnv.json$rows[[1]]; colnames(mtx.cnv) <- cnv.json$cols[[1]]
				cnv.samples <- grep(regex, cnv.json$cols[[1]],  value=TRUE)
				cnv.index <- cnvTables[i,"id"]

				for(j in 1:nrow(mutTables)){
					mut.json <- fromJSON(paste(mutTables[j,"directory"],mutTables[j,"file"], sep="/")) 
					mtx.mut <- mut.json$data[[1]]
					rownames(mtx.mut) <- mut.json$rows[[1]]; colnames(mtx.mut) <- mut.json$cols[[1]]
   				mut.samples <- grep(regex, mut.json$cols[[1]],  value=TRUE)
					mut.index <- mutTables[j,"id"]
   					
   				samples <- unique(cnv.samples, mut.samples)
	   			sample_similarity <- calculateSampleSimilarityMatrix(
	   										mtx.mut, mtx.cnv,
	   										copyNumberValues=gistic.scores,
   											genes = goi, samples=samples)
   													 
					if(!is.na(threshold)){
						outliers <- names(which(sample_similarity[,1]<threshold))
						sample_similarity <- sample_similarity[setdiff(rownames(sample_similarity), outliers), ]
					}
   					
					## ----- Save  ------
	   			index <- get.new.collection.index(datasetName, dataType)
	   			processName <- paste(c(process,geneset_name),collapse=";")
	   			outputFile <- paste(datasetName, dataType, index, processName, sep="_")
	   			parent <- list(c(datasetName, "cnv", cnv.index),c(datasetName, "mut01", mut.index))
	   			collection <- data.frame(id=index, date=date,process=list(c(process,geneset_name)), directory=output_directory, file=outputFile)
	   			collection$parent <- parent
	   			add.new.collection(datasetName, dataType, collection)
	   			
					mds.list<- lapply(rownames(sample_similarity), function(name) c(x=sample_similarity[name,"x"], y=sample_similarity[name, "y"]))
					names(mds.list) <- rownames(sample_similarity)
					os.data.save(mds.list, output_directory, outputFile, format="JSON")

				} # mut files
			} #cnv files
		} # for diseaseName	
}
#---------------------------------------------------------
get.new.collection.index <- function(datasetName, dataTypeName){
  
  if(nrow(Manifest) == 0) return(1)
  
  dataObj <- subset(Manifest, dataset == datasetName && dataType == dataTypeName)
  if(nrow(dataObj) == 0) return(1)
  
  return(nrow(dataObj$collection) +1)
}
#---------------------------------------------------------
add.new.collection <- function(datasetName, dataTypeName, collection){
  
  if(nrow(Manifest) == 0){	
    newCollection <- data.frame(dataset=datasetName, dataType=dataTypeName)
    newCollection$collections <- list(collection)
    Manifest <<- newCollection
  }
  
  dataObj <- subset(Manifest, dataset == datasetName & dataType == dataTypeName)
  if(nrow(dataObj) == 1){
    Manifest[Manifest$dataset==datasetName & Manifest$dataType ==dataTypeName,"collection"] <<- c(dataObj$collection, collection)
  }
  if(nrow(dataObj) == 0){	
    newCollection <- data.frame(dataset=datasetName, dataType=dataTypeName)
    newCollection$collections <- list(collection)
    Manifest <<- rbind(Manifest, newCollection)
  }
  stop(printf("add.new.collection found %d instances of dataset %s and dataType %s", length(dataObj), datasetName, dataTypeName))
  
}
#----------------------------------------------------------------------------------------------------
get.network_edges <- function(mtx,samples, genes, edgeTypes){

  if(all(is.na(samples))) samples <- colnames(mtx)
  if(all(is.na(genes))) genes <- rownames(mtx)
	
  samples <- intersect(samples, colnames(mtx))
  genes <- intersect(genes, rownames(mtx))
  
  mtx <- mtx[genes, samples]
  rows <- rownames(mtx); cols <- colnames(mtx)
  allEdges <- c()
  
  for(edgeName in names(edgeTypes)){
  	matchingIndex <- which(mtx==edgeTypes[[edgeName]], arr.ind=T)
	  edgeMap <- apply(matchingIndex, 1, function(matchPair){
	    c(edgeName, rows[matchPair[1]], cols[matchPair[2]])
	  })
	  allEdges <- rbind(allEdges, t(edgeMap))
  }
  colnames(allEdges) <- c("m", "g", "p")
  
  return(allEdges)
}
#----------------------------------------------------------------------------------------------------
save.edge.files <- function(edgePairs, outputDirectory, datasetName, dataType, index, genesetName){

	## get and save node degrees
	node1_counts <- as.data.frame(table(edgePairs[,2]))  
	colnames(node1_counts) <- NULL

	node2_counts <- as.data.frame(table(edgePairs[,3]))  
	colnames(node2_counts) <- NULL
	
	edgeFile <- paste("edges"     ,datasetName, dataType, index, genesetName, sep="_")
	geneDegreeFile <- paste("geneDegree",datasetName, dataType, index, genesetName, sep="_")
	ptDegreeFile <- paste("ptDegree"  ,datasetName, dataType, index, genesetName, sep="_")

	os.data.save(edgePairs,    outputDirectory, edgeFile , format="JSON")
	os.data.save(node1_counts, outputDirectory, geneDegreeFile , format="JSON")
	os.data.save(node2_counts, outputDirectory, ptDegreeFile, format="JSON")
		
	return( c(edgeFile, geneDegreeFile, ptDegreeFile) )

}
#----------------------------------------------------------------------------------------------------
get.edgePairs <- function(collection, genesetName, ...){				
  
  goi <- genesets[[genesetName]]
 
    File <- collection$file
    dataObj <- fromJSON(paste(collections$directory, File, sep="/")) 
    mtx <- dataObj$data[[1]]
    rownames(mtx) <- dataObj$rows[[1]]
    colnames(mtx) <- dataObj$cols[[1]]
    
    ## get and save edge pairs
    edgePairs <- get.network_edges(mtx, samples=NA, ...)

  return(edgePairs)
}

#----------------------------------------------------------------------------------------------------
run.batch.network_edges <- function(manifest_file, output_directory="./"){

    # Load Input File 
    datasets <- fromJSON(manifest_file)
    dataType <- "network"
    date <- Sys.Date()
    
    newManifest <-data.frame()

    # Loop for each dataset
    datasetNames <- unique(datasets$dataset)
		for (datasetName in datasetNames){
		  molTables <- subset(datasets, dataset==datasetName && dataType %in% c("cnv", "mut01"))

		  cnvTables <- subset(molTables, dataType == "cnv")
		  mutTables <- subset(molTables, dataType == "mut01")

		  if(nrow(cnvTables)==0 & nrow(mutTables) ==0) next;
		  cat(dataset, "\n")		  
		  
			# save Edge sets for each Gene set
			for(genesetName in names(genesets)){			
	
				EdgeList <- list()
			  	 
				# save individual CNV edges
				if(nrow(cnvTables) >0){
				  
				  for(i in 1:nrow(cnvTables)){
					  newEdges <- get.edgePairs(cnvTables[i,], genesetName, edgeTypes=list("-2"="-2", "-1"="-1", "1"="1", "2"="2"))

					  index <- get.new.collection.index(datasetName, dataType)
					  edgeFiles <- save.edge.files(newEdges,output_directory, datasetName, dataType,index, genesetName)
					  parent <- list(c(datasetName, "cnv", index))
					  collection <- data.frame(id=index, process="cnv", date=date,directory=output_directory, file=edgeFiles)
					  collection$parent <- parent
					  add.new.collection(datasetName, dataType, collection)
					  
					  EdgeList$cnv[[as.character(index)]] <- newEdges
				  }
				}
				# save individual MUT edges				
				if(nrow(mutTables)==0)next;

				for(i in 1:nrow(mutTables)){
				  newEdges <- get.edgePairs(mutTables[i,], genesetName, edgeTypes=list("0"="1"))
				  
				  index <- get.new.collection.index(datasetName, dataType)
				  edgeFiles <- save.edge.files(newEdges,output_directory, datasetName, dataType,index, genesetName)
				  parent <- list(c(datasetName, "mut", index))
				  collection <- data.frame(id=index, process="mut01", date=date,directory=output_directory, file=edgeFiles)
				  collection$parent <- parent
				  add.new.collection(datasetName, dataType, collection)
				  
				  EdgeList$mut[[as.character(index)]] <- newEdges
				}

				# save all combinations of MUT & CNV edges
				numCNV = length(EdgeList$cnv)
				numMut = length(EdgeList$mut)
				process <- c("cnv;mut01", genesetName)
				processName <- paste(process, collapse=";")
				
				for(k in 1:numCNV){
					cnvEdges <- EdgeList$cnv[[k]]
					
					for(m in 1:numMut){
						mutEdges <- EdgeList$mut[[m]]
																	
						allEdges <- rbind(cnvEdges, mutEdges)
						#compIndex = (k+m)+((k-1)*numCNV)+m
						index <- get.new.collection.index(datasetName, dataType)
						edgeFiles <- save.edge.files(allEdges, output_directory, datasetName, dataType, index=index, processName)
						parent <- list(c(datasetName, "cnv", names(EdgeList$cnv)[k]),c(datasetName, "mut01", names(EdgeList$mut)[m]) )
						collection <- data.frame(id=index, process=processName, date=date,directory=output_directory, file=edgeFiles)
						collection$parent <- parent
						add.new.collection(datasetName, dataType, collection)
					}
				}
				
				
			} # for genesetName
			
 		} # for diseaseName	

}


#----------------------------------------------------------------------------------------------------

if("mds" %in% commands)
	run.batch.patient_similarity(molecular_manifest,geneset_name="oncoVogel274", outputDirectory = mds_output_directory)
		# calculate patient similarity
		# save json

if("edges" %in% commands) 
	run.batch.network_edges(molecular_manifest, outputDirectory=network_output_directory)
		# map edges for all patients between CNV/Mut and Geneset tables
