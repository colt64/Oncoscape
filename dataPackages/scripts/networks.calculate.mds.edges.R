library(org.Hs.eg.db)
library(jsonlite)

rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

source("common.R")

#commands <- c("mds", "edges")
commands <- c("mds")
#commands <- c("edges")
args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

#--------------------------------- Configuration -----------------------------#


molecular_manifest <- "../manifests/os.import.molecular.ucsc.manifest.json"

network_output_directory <- "../data/molecular/edges/"
mds_output_directory <- "../data/molecular/mds/"

geneset_file <- "../data/molecular/hg19/hg19_genesets_1_hgnc.json"
genesets <- fromJSON(geneset_file)

date <- as.character(Sys.Date())

#----------------------------------------------------------------------------------------------------
getGeneSet <- function(geneset_name){
	stopifnot(geneset_name %in% names(genesets))
	return(genesets[[geneset_name]])
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

     return(tbl.pos)
}

#----------------------------------------------------------------------------------------------------
save.pca<- function(Manifest, tbl, datasetName, dataType, geneset=NA, output_directory=output_directory, ...){

#   mtx[is.na(mtx)] <- 0.0
  	## ----- Configuration ------
	process <- data.frame(calculation="prcomp", geneset= geneset)
	process$input=dataType
	processName <- paste(unlist(process), collapse="-")
	process$center="TRUE"; process$scale="TRUE"

	for(i in 1:nrow(tbl)){
	  	collection <- tbl[i, "collections"][[1]]
		tbl.json <- fromJSON(paste(collection$directory, collection$file,".json", sep="")) 
		mtx <- tbl.json$data[[1]]
		dimnames(mtx) <- list(tbl.json$rows[[1]], tbl.json$cols[[1]])
		tbl.index <- tbl.json$id
   
		if(tbl.json$rowType == "gene")
		  mtx <- t(mtx)
		
		column.sums <- colSums(mtx, na.rm=TRUE)
		removers <- as.integer(which(column.sums == 0))
		if(length(removers) > 0) {
		   printf("removing %d columns", length(removers))
		   mtx <- mtx[, -removers]
		} # if removers

		if(!is.na(geneset)){
			genes <- getGeneSet(geneset)
			mtx <- mtx[, intersect(colnames(mtx), genes)]
		}
   
	   PCs <- tryCatch(
		  prcomp(mtx,center=T,scale=T),
		  error=function(error.message){
			 print(error.message)
			 stop("error with PCA calculation.  See R log");
			 })
   
	   if(all(is.na(PCs)))
		   stop("error with PCA calculation.  See R log");
	
	   parent <- list(c(datasetName, dataType, tbl.index))
	   
	   scores <- PCs$x
     result <- list(rowType="samples", colType="PC", rows=rownames(scores), cols=colnames(scores), data=scores)
	   ## ----- Save  ------
     Manifest <- save.collection(Manifest=Manifest, dataset=datasetName, dataType="pcaScores", result=result,
                     parent=parent, process=process,processName=processName, outputDirectory=outputDirectory)
 
	   loadings <- PCs$rotation
	   result <- list(rowType="genes", colType="PC", rows=rownames(loadings), cols=colnames(loadings), data=loadings)
	   Manifest <- save.collection(Manifest=Manifest, dataset=datasetName, dataType="pcaLoadings", result=result,
	                               parent=parent, process=process,processName=processName, outputDirectory=outputDirectory)
	   
	   	#importance <- summary(PCs)$importance   
	   	
	}# for each collection

	return(Manifest)
}


#----------------------------------------------------------------------------------------------------
save.mds.innerProduct <- function(Manifest,datasetName, tbl1, tbl2, tbl1Type, tbl2Type, geneset, output_directory=output_directory, ...){
    ## ----- MDS on All Combinations of CNV and MUT Tables ------

  	## ----- Configuration ------
	genes = getGeneSet(geneset)
  	dataType <- "mds"

	regex = ".01$"; threshold = NA;
	if(datasetName == "laml"){        regex = "-03$|-09$";
	} else if(datasetName == "luad"){ regex = "TCGA-(17)^-\\d{4}-01$" }

	if(datasetName == "brca" | datasetName == "brain")  threshold = -1e-04
  
  	process <- data.frame(calculation="mds", geneset= geneset)
  	process$input=list( c(tbl1Type, tbl2Type))
  	processName <- paste(unlist(process), collapse="-")
  	process$regex=regex; process$threshold=threshold

	for(i in 1:nrow(tbl1)){
	  collection <- tbl1[i, "collections"][[1]]
	  tbl1.json <- fromJSON(paste(collection$directory, collection$file,".json", sep="")) 
		mtx.tbl1 <- tbl1.json$data[[1]]
		rownames(mtx.tbl1) <- tbl1.json$rows[[1]]; colnames(mtx.tbl1) <- tbl1.json$cols[[1]]
		tbl1.samples <- grep(regex, tbl1.json$cols[[1]],  value=TRUE)
		tbl1.index <- collection$id

		for(j in 1:nrow(tbl2)){
		  collection <- tbl2[i, "collections"][[1]]
		  tbl2.json <- fromJSON(paste(collection$directory, collection$file,".json", sep="")) 
		  mtx.tbl2 <- tbl2.json$data[[1]]
			rownames(mtx.tbl2) <- tbl2.json$rows[[1]]; colnames(mtx.tbl2) <- tbl2.json$cols[[1]]
			tbl2.samples <- grep(regex, tbl2.json$cols[[1]],  value=TRUE)
			tbl2.index <-collection$id
			
			samples <- unique(tbl1.samples, tbl2.samples)
			sample_similarity <- calculateSampleSimilarityMatrix(
									mtx.tbl1, mtx.tbl2,samples=samples, ...)
											 
			if(!is.na(threshold)){
				outliers <- names(which(sample_similarity[,1]<threshold))
				sample_similarity <- sample_similarity[setdiff(rownames(sample_similarity), outliers), ]
			}

			parent <- list(c(datasetName, tbl1Type, tbl1.index),c(datasetName, tbl2Type, tbl2.index))
			mds.list<- lapply(rownames(sample_similarity), function(name) data.frame(x=sample_similarity[name,"x"], y=sample_similarity[name, "y"]))
			names(mds.list) <- rownames(sample_similarity)
			Manifest <- save.collection(Manifest=Manifest, dataset=datasetName, dataType=dataType, result=mds.list,
			                            parent=parent, process=process,processName=processName, outputDirectory=outputDirectory)
			
		} # mut files
	} #cnv files

	return(Manifest)
}


#----------------------------------------------------------------------------------------------------
run.batch.patient_similarity <- function(manifest_file, geneset_name=NA, output_directory="./"){

  Manifest <- data.frame()
  
  datasets <- fromJSON(manifest_file)
  gistic.scores <-c(-2,-1,1, 2)
  
  # Loop for each dataset
  datasetNames <- unique(datasets$dataset)
  for (datasetName in datasetNames){
      molTables <- subset(datasets, dataset==datasetName)
    
      cnvTables <- subset(molTables, dataType == "cnv")
      mutTables <- subset(molTables, dataType == "mut01")
    
      if(nrow(cnvTables)==0 | nrow(mutTables) ==0) next;
      cat(datasetName, "\n")		  
    
	  Manifest <- save.mds.innerProduct(Manifest=Manifest, datasetName=datasetName, tbl1=cnvTables, tbl2=mutTables, tbl1Type="cnv", tbl2Type="mut01", 
	                                    copyNumberValues=gistic.scores, geneset = geneset_name, output_directory=output_directory)
	  Manifest <- save.pca(Manifest=Manifest, tbl=cnvTables,datasetName=datasetName, dataType="cnv", geneset = geneset_name, output_directory=output_directory)

	} # for diseaseName	
  
  return(Manifest)
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
save.edge.files <- function(Manifest, dataset, result,
                            parent, process,processName, outputDirectory="./"){

  Manifest <- save.collection(Manifest=Manifest, dataset=datasetName, dataType="edges", result=result,
                              parent=parent, process=process,processName=processName, outputDirectory=output_directory)
  
  node1_counts <- as.data.frame(table(result[,2]))  
  colnames(node1_counts) <- NULL
  Manifest <- save.collection(Manifest=Manifest, dataset=datasetName, dataType="geneDegree", result=node1_counts,
                              parent=parent, process=process,processName=processName, outputDirectory=output_directory)
  
  node2_counts <- as.data.frame(table(result[,3]))  
  colnames(node2_counts) <- NULL
  Manifest <- save.collection(Manifest=Manifest, dataset=datasetName, dataType="ptDegree", result=node2_counts,
                              parent=parent, process=process,processName=processName, outputDirectory=output_directory)

	return( Manifest )

}
#----------------------------------------------------------------------------------------------------
get.edgePairs <- function(collection, genesetName, ...){				
  
    goi <- genesets[[genesetName]]
 
    File <- collection$file
    dataObj <- fromJSON(paste(collection$directory, File,".json", sep="")) 
    mtx <- dataObj$data[[1]]
    rownames(mtx) <- dataObj$rows[[1]]
    colnames(mtx) <- dataObj$cols[[1]]
    
    ## get and save edge pairs
    edgePairs <- get.network_edges(mtx, samples=NA, genes=goi, ...)

  return(edgePairs)
}

#----------------------------------------------------------------------------------------------------
run.batch.network_edges <- function(manifest_file, output_directory="./"){

    # Load Input File 
    datasets <- fromJSON(manifest_file)
    dataType <- "network"

    Manifest <-data.frame()

    # Loop for each dataset
    datasetNames <- unique(datasets$dataset)
		for (datasetName in datasetNames){
		  molTables <- subset(datasets, dataset==datasetName)

		  cnvTables <- subset(molTables, dataType == "cnv")
		  mutTables <- subset(molTables, dataType == "mut01")

		  if(nrow(cnvTables)==0 & nrow(mutTables) ==0) next;
		  cat(datasetName, "\n")		  
		  
			# save Edge sets for each Gene set
			for(genesetName in names(genesets)){			
	
				EdgeList <- list()
			  	 
				# save individual CNV edges
				if(nrow(cnvTables) >0){
				  
				  for(i in 1:nrow(cnvTables)){
				    collection <- cnvTables[i,"collections"][[1]]
					  newEdges <- get.edgePairs(collection, genesetName, edgeTypes=list("-2"="-2", "-1"="-1", "1"="1", "2"="2"))

					  parent <- list(c(datasetName, "cnv", collection$id))
					  process <- list(edgeType="cnv", geneset= genesetName); processName=paste(process, collapse="-")
					  Manifest <- save.edge.files(Manifest=Manifest, dataset=datasetName, result=result,
					                              parent=parent, process=process,processName=processName, outputDirectory=output_directory)				  
					  
					  EdgeList$cnv[[as.character(index)]] <- newEdges
				  }
				}
				# save individual MUT edges				
				if(nrow(mutTables)==0)next;

				for(i in 1:nrow(mutTables)){
				  collection <- mutTables[i, "collections"][[1]]
				  newEdges <- get.edgePairs(collection, genesetName, edgeTypes=list("0"="1"))
				  
				  parent <- list(c(datasetName, "mut", collection$id))
				  process <- list(edgeType="mut01", geneset= genesetName); processName=paste(process, collapse="-")
				  
          Manifest <- save.edge.files(Manifest=Manifest, dataset=datasetName, result=result,
                                      parent=parent, process=process,processName=processName, outputDirectory=output_directory)				  

				  EdgeList$mut[[as.character(index)]] <- newEdges
				}

				# save all combinations of MUT & CNV edges
				numCNV = length(EdgeList$cnv)
				numMut = length(EdgeList$mut)
				process <- list(edgeType=c("cnv", "mut01"), geneset= genesetName)
				processName <- paste(unlist(process), collapse="-")
				
				for(k in 1:numCNV){
					cnvEdges <- EdgeList$cnv[[k]]
					
					for(m in 1:numMut){
						mutEdges <- EdgeList$mut[[m]]
																	
						allEdges <- rbind(cnvEdges, mutEdges)

						parent <- list(c(datasetName, "network", names(EdgeList$cnv)[k]),c(datasetName, "network", names(EdgeList$mut)[m]) )
						Manifest <- save.edge.files(Manifest=Manifest, dataset=datasetName, result=result,
						                            parent=parent, process=process,processName=processName, outputDirectory=output_directory)				  
						
					}
				}
				
				
			} # for genesetName
			
 		} # for diseaseName	

      return(Manifest)
}


#----------------------------------------------------------------------------------------------------

if("mds" %in% commands){
	Manifest_mds <- run.batch.patient_similarity(molecular_manifest,geneset_name="oncoVogel274", output_directory = mds_output_directory)
		# calculate patient similarity
		# save json
	os.data.save(
	  df = Manifest_mds,
	  directory="../manifests",
	  file= paste("os","mds","manifest", sep="."),
	  format = "JSON") 
}

if("edges" %in% commands){
	Manifest_edges <- run.batch.network_edges(molecular_manifest, output_directory=network_output_directory)
		# map edges for all patients between CNV/Mut and Geneset tables
	os.data.save(
	  df = Manifest_edges,
	  directory="../manifests",
	  file= paste("os","edges","manifest", sep="."),
	  format = "JSON") 
	
}


