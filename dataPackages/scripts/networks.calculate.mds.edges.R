library(org.Hs.eg.db)
library(jsonlite)

rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

source("common.R")

commands <- c("mds", "edges")
#commands <- c("mds")
#commands <- c("edges")
args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

#--------------------------------- Configuration -----------------------------#
source("common.R")


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
calculateSampleSimilarityMatrix <- function (mut, cn, samples=NA, genes=NA) {

  cat("--- sample similarity matrix\n")
  
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

  cat("-calculating pca\n")
  
  ## ----- Configuration ------
	process <- data.frame(calculation="prcomp", geneset= geneset)
	process$input=dataType
	processName <- paste(unlist(process), collapse="-")
	process$center="TRUE"; process$scale="TRUE"
	process <- list(process)

	coll <- tbl$collections[[1]]
	
	
	for(i in 1:nrow(coll)){
	  collection <- coll[i,]
		tbl.json <- fromJSON(paste(collection$directory, collection$file,".json", sep="")) 
		mtx <- tbl.json$data[[1]]
		dimnames(mtx) <- list(tbl.json$rows[[1]], tbl.json$cols[[1]])
		tbl.index <- collection$id
   
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
     Manifest <- save.collection(Manifest=Manifest, dataset=datasetName, dataType="pcaScores",source=collection$source, result=result,
                     parent=parent, process=process,processName=processName, outputDirectory=output_directory)
 
	   loadings <- PCs$rotation
	   result <- list(rowType="genes", colType="PC", rows=rownames(loadings), cols=colnames(loadings), data=loadings)
	   Manifest <- save.collection(Manifest=Manifest, dataset=datasetName, dataType="pcaLoadings", source=collection$source, result=result,
	                               parent=parent, process=process,processName=processName, outputDirectory=output_directory)
	   
	   	#importance <- summary(PCs)$importance   
	   	
	}# for each collection

	return(Manifest)
}


#----------------------------------------------------------------------------------------------------
save.mds.innerProduct <- function(datasetName, tbl1, tbl2, tbl1Type, tbl2Type, geneset, ...){
    ## ----- MDS on All Combinations of CNV and MUT Tables ------

  cat("-calculating mds\n")
  
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
  	process <- list(process)

  	coll1 <- tbl1$collections[[1]]
  	coll2 <- tbl2$collections[[1]]
  	
	for(i in 1:nrow(coll1)){
	  tbl1.json <- fromJSON(paste(coll1[i,"directory"], coll1[i,"file"],".json", sep="")) 
		mtx.tbl1 <- tbl1.json$data[[1]]
		rownames(mtx.tbl1) <- tbl1.json$rows[[1]]; colnames(mtx.tbl1) <- tbl1.json$cols[[1]]
		tbl1.samples <- grep(regex, tbl1.json$cols[[1]],  value=TRUE)
		tbl1.index <- coll1[i, "id"]

		cat("--- tbl1: ", coll1[i, "file"], "\n")
		
		for(j in 1:nrow(coll2)){
		  tbl2.json <- fromJSON(paste(coll2[j , "directory"], coll2[j, "file"],".json", sep="")) 
		  mtx.tbl2 <- tbl2.json$data[[1]]
			rownames(mtx.tbl2) <- tbl2.json$rows[[1]]; colnames(mtx.tbl2) <- tbl2.json$cols[[1]]
			tbl2.samples <- grep(regex, tbl2.json$cols[[1]],  value=TRUE)
			tbl2.index <-coll2[j,"id"]

			cat("--- tbl2: ", coll2[j, "file"], "\n")
			
			samples <- unique(tbl1.samples, tbl2.samples)
			sample_similarity <- calculateSampleSimilarityMatrix(mtx.tbl1, mtx.tbl2,samples=samples, genes=genes)
											 
			if(!is.na(threshold)){
				outliers <- names(which(sample_similarity[,1]<threshold))
				sample_similarity <- sample_similarity[setdiff(rownames(sample_similarity), outliers), ]
			}

			parent <- list(list(c(datasetName, tbl1Type, tbl1.index),c(datasetName, tbl2Type, tbl2.index)))
			mds.list<- lapply(rownames(sample_similarity), function(name) data.frame(x=sample_similarity[name,"x"], y=sample_similarity[name, "y"]))
			names(mds.list) <- rownames(sample_similarity)
			save.collection(mongo, dataset=datasetName, dataType=dataType,source=c(coll1[i,"source"], coll2[j,"source"]), result=mds.list,
			                            parent=parent, process=process,processName=processName)
			
		} # mut files
	} #cnv files

}


#----------------------------------------------------------------------------------------------------
run.batch.patient_similarity <- function(datasets, geneset_name=NA, output_directory="./"){

  gistic.scores <-c(-2,-1,1, 2)
  
  # Loop for each dataset
  datasetNames <- unique(datasets$dataset)
  for (datasetName in datasetNames){
      molTables <- subset(datasets, dataset==datasetName)
    
      cnvTables <- subset(molTables, dataType == "cnv")
      mutTables <- subset(molTables, dataType == "mut01")
      rnaTables <- subset(molTables, dataType == "rna")
      protTables <- subset(molTables, dataType == "protein")
      
      if(nrow(cnvTables)==0 | nrow(mutTables) ==0) next;
      cat(datasetName, "\n")		  
    
	  save.mds.innerProduct(datasetName=datasetName, tbl1=cnvTables, tbl2=mutTables, tbl1Type="cnv", tbl2Type="mut01", 
	                                    copyNumberValues=gistic.scores, geneset = geneset_name)
	  save.pca(tbl=cnvTables,datasetName=datasetName, dataType="cnv", geneset = geneset_name)
	  save.pca(tbl=mutTables,datasetName=datasetName, dataType="mut01", geneset = geneset_name)
	  save.pca(tbl=rnaTables,datasetName=datasetName, dataType="rna", geneset = geneset_name)
	  save.pca(tbl=rnaTables,datasetName=datasetName, dataType="rna", geneset = NA)
	  save.pca(tbl=protTables,datasetName=datasetName, dataType="protein", geneset = geneset_name)
	  save.pca(tbl=protTables,datasetName=datasetName, dataType="protein", geneset = NA)
	  
	} # for diseaseName	
  
  
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
save.edge.files <- function(dataset, result, source, parent, process,processName){

  save.collection(mongo, dataset=dataset, dataType="edges",source=source, result=result,
                              parent=parent, process=process,processName=processName)
  
  node1_counts <- as.data.frame(table(result[,2]))  
  colnames(node1_counts) <- NULL
  save.collection(mongo, dataset=dataset, dataType="geneDegree",source=source, result=node1_counts,
                              parent=parent, process=process,processName=processName)
  
  node2_counts <- as.data.frame(table(result[,3]))  
  colnames(node2_counts) <- NULL
  save.collection(mongo, dataset=dataset, dataType="ptDegree", source=source, result=node2_counts,
                              parent=parent, process=process,processName=processName)
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
run.batch.network_edges <- function(datasets){

  cat("-calculating edges\n")
  
    # Load Input File 
#    datasets <- fromJSON(manifest_file)
    dataType <- "network"

    Manifest <-data.frame()

    # Loop for each dataset
    datasetNames <- unique(datasets$dataset)
		for (datasetName in datasetNames){
		  molTables <- subset(datasets, dataset==datasetName)

		  cnvTables <- subset(molTables, dataType == "cnv")$collections[[1]]
		  mutTables <- subset(molTables, dataType == "mut01")$collections[[1]]

		  if(nrow(cnvTables)==0 & nrow(mutTables) ==0) next;
		  cat(datasetName, "\n")		  
		  
			# save Edge sets for each Gene set
			for(genesetName in names(genesets)){			
	
				EdgeList <- list()
			  	 
				# save individual CNV edges
				if(nrow(cnvTables) >0){
				  
				  for(i in 1:nrow(cnvTables)){
				    collection <- cnvTables[i,]
					  newEdges <- get.edgePairs(collection, genesetName, edgeTypes=list("-2"="-2", "-1"="-1", "1"="1", "2"="2"))

					  parent <- list(c(datasetName, "cnv", collection$id))
					  process <- list(edgeType="cnv", geneset= genesetName); processName=paste(process, collapse="-")
					  save.edge.files(dataset=datasetName, result=newEdges, source=collection$source,
					                              parent=parent, process=list(process),processName=processName)				  
					  
					  EdgeList$cnv[[as.character(collection$id)]] <- newEdges
				  }
				}
				# save individual MUT edges				
				if(nrow(mutTables)==0)next;

				for(i in 1:nrow(mutTables)){
				  collection <- mutTables[i,]
				  newEdges <- get.edgePairs(collection, genesetName, edgeTypes=list("0"="1"))
				  
				  parent <- list(c(datasetName, "mut", collection$id))
				  process <- list(edgeType="mut01", geneset= genesetName); processName=paste(process, collapse="-")
				  
          save.edge.files(dataset=datasetName,source=collection$source, result=newEdges,
                                      parent=parent, process=list(process),processName=processName)				  

				  EdgeList$mut[[as.character(collection$id)]] <- newEdges
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

						parent <- list(c(datasetName, "cnv", names(EdgeList$cnv)[k]),c(datasetName, "mut01", names(EdgeList$mut)[m]) )
						source1 <- subset(datasets, dataset==datasetName & dataType=="cnv")$collections[[1]]
						source1 <- source1[source1$id==as.integer(names(EdgeList$cnv)[k]),"source"]
						source2 <- subset(datasets, dataset==datasetName & dataType=="mut")$collections[[1]]
						source2 <- source2[source2$id==as.integer(names(EdgeList$mut)[m]),"source"]
						save.edge.files(dataset=datasetName,source =list(c(source1,source2)),  result=allEdges,
						                            parent=list(parent), process=list(process),processName=processName)				  
						
					}
				}
				
				
			} # for genesetName
			
 		} # for diseaseName	

      return(Manifest)
}


#----------------------------------------------------------------------------------------------------
## must first initialize server (through shell >mongod)
mongo <- connect.to.mongo()

genesets <-     mongo.find.all(mongo, "oncoscape.hg19_genesets_hgnc_import", 
                                       query=list())

molecular_manifest <- mongo.find.all(mongo, "oncoscape.manifest", 
                                    query='{"dataType":{"$in":["cnv","mut01", "mut", "rna", "protein", "methylation"]}}')

if("mds" %in% commands){
	run.batch.patient_similarity(molecular_manifest,geneset_name="oncoVogel274")
		# calculate patient similarity
}

if("edges" %in% commands){
	run.batch.network_edges(molecular_manifest)
		# map edges for all patients between CNV/Mut and Geneset tables

}


