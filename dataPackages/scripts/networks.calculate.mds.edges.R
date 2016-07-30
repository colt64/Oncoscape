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
  match_name <- which(sapply(genesets, function(set){set$name ==geneset_name}))
  if(length(match_name) == 0)
    return(NA)
  
	return(genesets[[match_name]]$genes)
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
save.pca<- function(collection, geneset=NA, scaleFactor=NA){

  cat("-calculating pca\n")
  
  ## ----- Configuration ------
	process <- data.frame(calculation="prcomp", geneset= geneset)
	process$input=collection$dataType
	processName <- paste(unlist(process), collapse="-")
	process$center="TRUE"; process$scale="TRUE"
	process <- list(process)

	prev.run <- collection.exists(mongo, collection$dataset, dataType="pcaScores",
	                              source=collection$source,processName=processName)
	if(prev.run){
	  print("Skipping.")
	  return();
	}
	
	coll <- mongo.find.all(mongo, paste("oncoscape",collection$collection, sep="."))
	
	mtx <- convert.to.mtx(coll);
	rm(coll);
	
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
		  prcomp(na.omit(mtx),center=T,scale=T),
		  error=function(error.message){
			 print(error.message)
			 stop("error with PCA calculation.  See R log");
			 })
   
	   if(all(is.na(PCs)))
		   stop("error with PCA calculation.  See R log");
	
	   parent <- collection$`_id`
	   
	   scores <- PCs$x
	   if(!is.na(scaleFactor)){
	     chrDim <- get.chromosome.dimensions(scaleFactor) 
	     scores.list <- scaleSamplesToChromosomes(scores, chrDim)
	   }else{
	      colnames(scores) <- NULL
	      scores.list <- lapply(rownames(scores), function(name){ scores[name,1:3]})
	   }
	   
	   names(scores.list) <- rownames(scores)
	   importance <- summary(PCs)$importance   
	   propVar <- importance[2,] *100
	   names(propVar) <- NULL
	   result <- list(disease=collection$dataset, geneset=geneset,scale=scaleFactor, pc1=propVar[1], pc2=propVar[2] ,pc3=propVar[3],data=scores.list)
	   ## ----- Save  ------
     save.collection(mongo, dataset=collection$dataset, dataType="pcaScores",source=collection$source, result=list(result),
                     parent=parent, process=process,processName=processName)

#	   loadings <- PCs$rotation
#	   result <- list(rowType="genes", colType="PC", rows=rownames(loadings), cols=colnames(loadings), data=loadings)
#	   Manifest <- save.collection(mongo, dataset=collection$dataset, dataType="pcaLoadings", source=collection$source, result=result,
#	                               parent=parent, process=process,processName=processName)

}


#----------------------------------------------------------------------------------------------------
save.mds.innerProduct <- function(tbl1, tbl2, geneset=NA, scaleFactor=NA, ...){
    ## ----- MDS on All Combinations of CNV and MUT Tables ------

  cat("-calculating mds\n")
  
  ## ----- Configuration ------
  dataType <- "mds"
  datasetName <- tbl1$dataset
  process <- data.frame(calculation="mds", geneset= geneset)
  process$input=list( c(tbl1$dataType, tbl2$dataType))
  processName <- paste(unlist(process), collapse="-")
  
  prev.run <- collection.exists(mongo, dataset=datasetName, dataType=dataType,source=c(tbl1$source, tbl2$source),processName=processName)
  if(prev.run){
    print("Skipping.")
    return()
  }
  
	regex = ".01$"; threshold = NA;
	if(datasetName == "laml"){        regex = "-03$|-09$";
	} else if(datasetName == "luad"){ regex = "TCGA-(17)^-\\d{4}-01$" }
	process$regex=regex; process$threshold=threshold
	process <- list(process)
	
	if(datasetName == "brca" | datasetName == "brain")  threshold = -1e-04
  
 
  	coll1 <- mongo.find.all(mongo, paste("oncoscape",tbl1$collection, sep="."))
  	coll2 <- mongo.find.all(mongo, paste("oncoscape",tbl1$collection, sep="."))
  	
		mtx.tbl1 <- convert.to.mtx(coll1);
		mtx.tbl2 <- convert.to.mtx(coll2);

		rm(coll1); rm(coll2);
		
		tbl1.samples <- grep(regex, rownames(mtx.tbl1),  value=TRUE)
		tbl2.samples <- grep(regex, rownames(mtx.tbl2),  value=TRUE)
	
		if(is.na(geneset)){
		       genes <- intersect(colnames(mtx.tbl1), colnames(mtx.tbl2))
		}else{ genes = getGeneSet(geneset) }
		
			
			samples <- unique(tbl1.samples, tbl2.samples)
			sample_similarity <- calculateSampleSimilarityMatrix(t(mtx.tbl1), t(mtx.tbl2),samples=samples, genes=genes)
											 #expects rows as genes and cols as samples
			
			if(!is.na(threshold)){
				outliers <- names(which(sample_similarity[,1]<threshold))
				sample_similarity <- sample_similarity[setdiff(rownames(sample_similarity), outliers), ]
			}

			parent <- list(tbl1$`_id`, tbl2$`_id`)
			
			if(!is.na(scaleFactor)){
			  chrDim <- get.chromosome.dimensions(scaleFactor) 
			  mds.list <- scaleSamplesToChromosomes(sample_similarity, chrDim)
			}else{
			  mds.list<- lapply(rownames(sample_similarity), function(name) data.frame(x=sample_similarity[name,"x"], y=sample_similarity[name, "y"]))
			  names(mds.list) <- rownames(sample_similarity)
			}
			
			result <- list(type="cluster", dataset=tbl1$dataset, name=processName, scale=scaleFactor, data=mds.list)
			save.collection(mongo, dataset=datasetName, dataType=dataType,source=c(tbl1$source, tbl2$source), result=list(result),
			                            parent=parent, process=process,processName=processName)
			
}


#----------------------------------------------------------------------------------------------------
run.batch.patient_similarity <- function(datasets){

  gistic.scores <-c(-2,-1,1, 2)
  
  # Loop for each dataset
  for (collection in datasets){

    ## MDS
    if(collection$dataType =="cnv"){
      mut01_colls <- mongo.find.all(mongo, "oncoscape.manifest", 
                     query=list(dataset=collection$dataset, dataType="mut01"))
      for(mut01_coll in mut01_colls){
        save.mds.innerProduct(collection, mut01_coll, copyNumberValues=gistic.scores, geneset = NA)
        for(geneset in genesets){
          save.mds.innerProduct(collection, mut01_coll, copyNumberValues=gistic.scores, geneset = geneset$name)
        }
        
      }
    }
    else if(collection$dataType =="mut01"){
      cnv_colls <- mongo.find.all(mongo, "oncoscape.manifest", 
                                    query=list(dataset=collection$dataset, dataType="cnv"))
      for(cnv_coll in cnv_colls){
        save.mds.innerProduct(cnv_coll, collection, copyNumberValues=gistic.scores, geneset = NA)
        for(geneset in genesets){
          save.mds.innerProduct(cnv_coll, collection, copyNumberValues=gistic.scores, geneset = geneset$name)
        }
        
      }
      
    }
    
    ## PCA
      save.pca(collection, geneset = NA)
      for(geneset in genesets){
        save.pca(collection, geneset = geneset$name)
      }
	      

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

}


#----------------------------------------------------------------------------------------------------
## must first initialize server (through shell >mongod)
mongo <- connect.to.mongo()

genesets <-     mongo.find.all(mongo, "oncoscape.hg19_genesets_hgnc_import", query=list())

molecular_manifest <- mongo.find.all(mongo, "oncoscape.manifest", 
                                    query='{"dataType":{"$in":["cnv","mut01", "rna", "protein", "methylation"]}}')

if("mds" %in% commands){
	run.batch.patient_similarity(molecular_manifest)
		# calculate patient similarity
}

if("edges" %in% commands){
	run.batch.network_edges(molecular_manifest)
		# map edges for all patients between CNV/Mut and Geneset tables

}


