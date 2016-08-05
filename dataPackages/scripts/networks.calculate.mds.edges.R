library(org.Hs.eg.db)
library(jsonlite)

rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
source("common.R")

commands <- c("mds", "edges")
#commands <- c("mds")
#commands <- c("edges")
<<<<<<< HEAD
=======
#commands <- c("mds", "edges")
commands <- c("edges")
>>>>>>> develop
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

#--------------------------------- Configuration -----------------------------#
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
source("common.R")


date <- as.character(Sys.Date())

#----------------------------------------------------------------------------------------------------
getGeneSet <- function(geneset_name){
  match_name <- which(sapply(genesets, function(set){set$name ==geneset_name}))
  if(length(match_name) == 0)
    return(NA)
  
	return(genesets[[match_name]]$genes)
}

<<<<<<< HEAD
=======

molecular_manifest <- "../manifests/os.import.2016-06-24.manifest.json"

network_output_directory <- "../data/molecular/edges/"
mds_output_directory <- "../data/molecular/mds/"

geneset_file <- "../data/molecular/hg19/hg19_genesets_1_hgnc.json"
genesets <- fromJSON(geneset_file)

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
  write(toJSON(df, pretty=TRUE, digits=I(8)), file=paste(outFile,".json", sep = "") )
  
}
>>>>>>> develop
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de

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
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
calculateSampleSimilarityMatrix <- function (mut, cn, samples=NA, genes=NA) {

  cat("--- sample similarity matrix\n")
  
<<<<<<< HEAD
=======
calculateSampleSimilarityMatrix <- function (mut, cn, samples=NA, genes=NA, copyNumberValues=c(-2, 2), threshold=NA) {

>>>>>>> develop
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
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
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
	tmp <- apply(cn, 1, function(x) any(is.na(x)))
	if(length(which(tmp))>0)
	  cn <- cn[-which(tmp),]
	
<<<<<<< HEAD
=======

	stopifnot(all(sort(unique(as.integer(mut))) == c(0,1)))
     
    cn[!cn %in% copyNumberValues] <- 0

>>>>>>> develop
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
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

<<<<<<< HEAD
<<<<<<< HEAD
=======

#	 ptIDs <- canonicalizePatientIDs(obj@pkg, rownames(tbl.pos))
#	 tbl.pos <- tbl.pos[!duplicated(ptIDs),]
#     rownames(tbl.pos) <- ptIDs[!duplicated(ptIDs)]
     
>>>>>>> develop
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
     return(tbl.pos)
}

#----------------------------------------------------------------------------------------------------
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
save.pca<- function(collection, geneset=NA, scaleFactor=NA){

  cat("-calculating pca\n")
  
  ## ----- Configuration ------
  genesetName <- geneset
  if(is.na(genesetName)) genesetName = "All Genes"
	process <- data.frame(calculation="prcomp", geneset= genesetName)
	process$input=collection$dataType
	outputName <- paste(unlist(process), collapse="-")
	
	process$scale=scaleFactor
	processName <- paste(unlist(process), collapse="-")
	process$center="TRUE"; process$scaled="TRUE"
	process <- list(process)

	prev.run <- collection.exists(mongo, collection$dataset, dataType="pcaScores",
	                              source=collection$source,processName=outputName)
	if(prev.run){
	  print("Skipping.")
	  return();
	}
	
	coll <- mongo.find.all(mongo, paste("oncoscape",collection$collection, sep="."))
	
	mtx <- convert.to.mtx(coll, format="as.numeric");
	rm(coll);

	if(!is.na(geneset)){
	  genes <- getGeneSet(geneset)
	  mtx <- mtx[, intersect(colnames(mtx), genes)]
	}
	
	column.sums <- colSums(mtx, na.rm=TRUE)
	removers <- as.integer(which(column.sums == 0))
	removers <- c(removers, which(apply(mtx, 2, var)== 0))
	if(length(removers) > 0) {
		   printf("removing %d columns", length(removers))
		   mtx <- mtx[, -removers]
	} # if removers

	if(any(dim(mtx)==0)){
	  print("WARNING: mtx is singular.  PCA not computed")
	  return();
	}
	  
	
	   PCs <- tryCatch(
		  prcomp(na.omit(mtx),center=T,scale=T),
		  error=function(error.message){
			 print("ERROR: PRCOMP!")
		    print(error.message)
			 return(NA);
			 })
   
	   
	   if(all(is.na(PCs)))
		   return();
	
	   parent <- collection$`_id`
	   
	   scores <- PCs$x
	   colnames(scores) <- NULL
	   importance <- summary(PCs)$importance   
	   propVar <- round(importance[2,] *100, digits=2)
	   names(propVar) <- NULL

	   
	   ## ----- Save Raw ------
	   scores.list <- lapply(rownames(scores), function(name){ scores[name,1:3]})
	   names(scores.list) <- rownames(scores)
	   process$scale = NA
	   result <- list(disease=collection$dataset,type=collection$dataType, geneset=genesetName,scale=NA, pc1=propVar[1], pc2=propVar[2] ,pc3=propVar[3],data=scores.list)
     save.collection(mongo, dataset=collection$dataset, dataType="pcaScores",source=collection$source, result=list(result),
                     parent=parent, process=process,processName=outputName)

     ## ----- Save Scaled  ------
     if(!is.na(scaleFactor)){
	     chrDim <- get.chromosome.dimensions(scaleFactor) 
	     pc3 <- scores[,1:3]; colnames(pc3) <- c("x", "y", "z")
	     scores.list <- scaleSamplesToChromosomes(pc3, chrDim)
	     names(scores.list) <- rownames(scores)
	     process$scale = scaleFactor
	     result <- list(disease=collection$dataset,type=collection$dataType, geneset=genesetName,scale=scaleFactor, pc1=propVar[1], pc2=propVar[2] ,pc3=propVar[3],data=scores.list)
	     save.collection(mongo, dataset=collection$dataset, dataType="pcaScores",source=collection$source, result=list(result),
	                     parent=parent, process=process,processName=processName)
	     
     }
#	   loadings <- PCs$rotation
#	   result <- list(rowType="genes", colType="PC", rows=rownames(loadings), cols=colnames(loadings), data=loadings)
#	   Manifest <- save.collection(mongo, dataset=collection$dataset, dataType="pcaLoadings", source=collection$source, result=result,
#	                               parent=parent, process=process,processName=processName)

}


#----------------------------------------------------------------------------------------------------
save.mds.innerProduct <- function(tbl1, tbl2, geneset=NA, scaleFactor=NA, ...){
    ## ----- MDS on All Combinations of CNV and MUT Tables ------

  if(tbl1$source != tbl2$source){
    print("currently not computing mds based on different sources")
    return()
  }
  
  cat("-calculating mds\n")
  
  ## ----- Configuration ------
  dataType <- "mds"
  genesetName <- geneset
  if(is.na(genesetName)) genesetName = "All Genes"
  datasetName <- tbl1$dataset
  process <- list(calculation="mds", geneset= genesetName)
  process$input=list( tbl1$dataType, tbl2$dataType)
  outputName <- paste(unlist(process), collapse="-")

  process$scale=scaleFactor
  processName <- paste(unlist(process), collapse="-")
  
  prev.run <- collection.exists(mongo, dataset=datasetName, dataType=dataType,source=c(tbl1$source, tbl2$source),processName=processName)
  if(prev.run){
    print("Skipping.")
    return()
  }
  
	regex = "-01$"; threshold = NA;
	if(datasetName == "laml"){        regex = "-03$|-09$";
	} else if(datasetName == "luad"){ regex = "TCGA-(17)^-\\d{4}-01$" }
	process$regex=regex; process$threshold=threshold

	if(datasetName == "brca" | datasetName == "brain")  threshold = -1e-04
  
 
  	coll1 <- mongo.find.all(mongo, paste("oncoscape",tbl1$collection, sep="."))
  	coll2 <- mongo.find.all(mongo, paste("oncoscape",tbl2$collection, sep="."))
  	
		mtx.tbl1 <- convert.to.mtx(coll1, format="as.numeric");
		mtx.tbl2 <- convert.to.mtx(coll2, format="as.numeric");

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
			
		  mds.list<- lapply(rownames(sample_similarity), function(name) data.frame(x=sample_similarity[name,"x"], y=sample_similarity[name, "y"]))
		  names(mds.list) <- rownames(sample_similarity)

		  process$scale = NA
		  process <- list(process)
			result <- list(type="cluster", dataset=tbl1$dataset, name=outputName, scale=NA, data=mds.list)
			save.collection(mongo, dataset=datasetName, dataType=dataType,source=c(tbl1$source, tbl2$source), result=list(result),
			                            parent=parent, process=process,processName=outputName)

			if(!is.na(scaleFactor)){
			    process[[1]]$scale = scaleFactor
			    chrDim <- get.chromosome.dimensions(scaleFactor) 
			  mds.list <- scaleSamplesToChromosomes(sample_similarity, chrDim, dim.names=c("x", "y"))
			  result <- list(type="cluster", dataset=tbl1$dataset, name=outputName, scale=scaleFactor, data=mds.list)
			  save.collection(mongo, dataset=datasetName, dataType=dataType,source=c(tbl1$source, tbl2$source), result=list(result),
			                  parent=parent, process=process,processName=processName)
			}			
}


#----------------------------------------------------------------------------------------------------
run.batch.patient_similarity <- function(datasets, scaleFactor=NA){

  gistic.scores <-c(-2,-1,1, 2)
  
  # Loop for each dataset
  for (collection in datasets){
    ## MDS
    if(collection$dataType =="cnv"){
      mut01_colls <- mongo.find.all(mongo, "oncoscape.manifest", 
                     query=list(dataset=collection$dataset, dataType="mut01"))
      for(mut01_coll in mut01_colls){
        save.mds.innerProduct(collection, mut01_coll, copyNumberValues=gistic.scores, geneset = NA, scaleFactor=scaleFactor)
        for(geneset in genesets){
          save.mds.innerProduct(collection, mut01_coll, copyNumberValues=gistic.scores, geneset = geneset$name, scaleFactor=scaleFactor)
        }
        
      }
    }
    else if(collection$dataType =="mut01"){
      cnv_colls <- mongo.find.all(mongo, "oncoscape.manifest", 
                                    query=list(dataset=collection$dataset, dataType="cnv"))
      for(cnv_coll in cnv_colls){
        save.mds.innerProduct(cnv_coll, collection, copyNumberValues=gistic.scores, geneset = NA, scaleFactor=scaleFactor)
        for(geneset in genesets){
          save.mds.innerProduct(cnv_coll, collection, copyNumberValues=gistic.scores, geneset = geneset$name, scaleFactor=scaleFactor)
        }
        
      }
      
    }
    
    ## PCA
      save.pca(collection, geneset = NA, scaleFactor=scaleFactor)
      for(geneset in genesets){
        save.pca(collection, geneset = geneset$name, scaleFactor=scaleFactor)
      }
	      

	} # for diseaseName	
  
<<<<<<< HEAD
=======
run.batch.patient_similarity <- function(manifest_file, geneset_name=NA, output_directory="./"){

  Manifest <- data.frame()
  
  # Load Input File 
  datasets <- fromJSON(manifest_file)
  dataType <- "mds"
  process <- list("calculation"="mds", "input"=c( "cnv", "mut01"), "geneset"= geneset_name)
  processName <- paste(unlist(process), collapse="-")

  gistic.scores <-c(-2,-1,1, 2)
  goi = getGeneSet(geneset_name)
  
  # Loop for each dataset
  datasetNames <- unique(datasets$dataset)
  for (datasetName in datasetNames){
      molTables <- subset(datasets, dataset==datasetName)
    
      cnvTables <- subset(molTables, dataType == "cnv")
      mutTables <- subset(molTables, dataType == "mut01")
    
      if(nrow(cnvTables)==0 | nrow(mutTables) ==0) next;
      cat(datasetName, "\n")		  
    
			## ----- Configuration ------
			regex = ".01$"; threshold = NA;
			if(datasetName == "laml"){        regex = "-03$|-09$";
			} else if(datasetName == "luad"){ regex = "TCGA-(17)^-\\d{4}-01$" }

			if(datasetName == "brca" | datasetName == "brain")  threshold = -1e-04
      
      process$regex=regex; process$threshold=threshold

			## ----- MDS on All Combinations of CNV and MUT Tables ------
			for(i in 1:nrow(cnvTables)){
			  collection <- cnvTables[i, "collections"][[1]]
				cnv.json <- fromJSON(paste(collection$directory, collection$file,".json", sep="")) 
				mtx.cnv <- cnv.json$data[[1]]
				rownames(mtx.cnv) <- cnv.json$rows[[1]]; colnames(mtx.cnv) <- cnv.json$cols[[1]]
				cnv.samples <- grep(regex, cnv.json$cols[[1]],  value=TRUE)
				cnv.index <- collection$id

				for(j in 1:nrow(mutTables)){
				  collection <- mutTables[i, "collections"][[1]]
				  mut.json <- fromJSON(paste(collection$directory, collection$file,".json", sep="")) 
				  mtx.mut <- mut.json$data[[1]]
					rownames(mtx.mut) <- mut.json$rows[[1]]; colnames(mtx.mut) <- mut.json$cols[[1]]
   				mut.samples <- grep(regex, mut.json$cols[[1]],  value=TRUE)
					mut.index <-collection$id
   					
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
	   			index <- get.new.collection.index(Manifest, datasetName, dataType)
	   			outputFile <- paste(datasetName, dataType, index, processName, sep="_")
	   			parent <- list(c(datasetName, "cnv", cnv.index),c(datasetName, "mut01", mut.index))
	   			collection <- data.frame(id=index, date=date, directory=output_directory, file=outputFile)
	   			collection$process <- list(process)
	   			collection$parent <- list(parent)
	   			Manifest <- add.new.collection(Manifest, datasetName, dataType, collection)
	   			
					mds.list<- lapply(rownames(sample_similarity), function(name) data.frame(x=sample_similarity[name,"x"], y=sample_similarity[name, "y"]))
					names(mds.list) <- rownames(sample_similarity)
					os.data.save(mds.list, output_directory, outputFile, format="JSON")

				} # mut files
			} #cnv files
		} # for diseaseName	
  
  return(Manifest)
}
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
>>>>>>> develop
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
  
}
#----------------------------------------------------------------------------------------------------
get.network_edges <- function(mtx,samples, genes, edgeTypes){

  if(all(is.na(samples))) samples <- colnames(mtx)
  if(all(is.na(genes))) genes <- rownames(mtx)
	
  samples <- intersect(samples, colnames(mtx))
  genes <- intersect(genes, rownames(mtx))
  
  mtx <- mtx[genes, samples]
  rows <- rownames(mtx); cols <- colnames(mtx)
<<<<<<< HEAD
<<<<<<< HEAD
  allEdges <- list()
=======
  allEdges <- c()
>>>>>>> develop
=======
  allEdges <- list()
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
  
  for(edgeName in names(edgeTypes)){
  	matchingIndex <- which(mtx==edgeTypes[[edgeName]], arr.ind=T)
	  edgeMap <- apply(matchingIndex, 1, function(matchPair){
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
	    list(m=edgeName, g=rows[matchPair[1]], p=cols[matchPair[2]])
	  })
	  allEdges <- c(allEdges, edgeMap)
  }
  #colnames(allEdges) <- c("m", "g", "p")
  
#  allEdges <- apply(allEdges, 1, function(row){row})
  return(allEdges)
}
#----------------------------------------------------------------------------------------------------
save.edge.files <- function(dataset, result, source, parent, process,processName){

  save.collection(mongo, dataset=dataset, dataType="edges",source=source, result=result,
                              parent=parent, process=process,processName=processName)
  
  temp <- as.list(table(sapply(result,function(edge) edge$p)))
  node1_counts <- lapply(names(temp), function(el) temp[el])
  save.collection(mongo, dataset=dataset, dataType="ptDegree",source=source, result=node1_counts,
                              parent=parent, process=process,processName=processName)
  
  temp <- as.list(table(sapply(result,function(edge) edge$g)))
  node2_counts <- lapply(names(temp), function(el) temp[el])
  save.collection(mongo, dataset=dataset, dataType="geneDegree", source=source, result=node2_counts,
                              parent=parent, process=process,processName=processName)
<<<<<<< HEAD
=======
	    c(edgeName, rows[matchPair[1]], cols[matchPair[2]])
	  })
	  allEdges <- rbind(allEdges, t(edgeMap))
  }
  colnames(allEdges) <- c("m", "g", "p")
  
  return(allEdges)
}
#----------------------------------------------------------------------------------------------------
save.edge.files <- function(edgePairs, outputDirectory, datasetName, dataType, index, process){

	## get and save node degrees
	node1_counts <- as.data.frame(table(edgePairs[,2]))  
	colnames(node1_counts) <- NULL

	node2_counts <- as.data.frame(table(edgePairs[,3]))  
	colnames(node2_counts) <- NULL
	
	edgeFile <- paste("edges"     ,datasetName, dataType, index, process, sep="_")
	geneDegreeFile <- paste("geneDegree",datasetName, dataType, index, process, sep="_")
	ptDegreeFile <- paste("ptDegree"  ,datasetName, dataType, index, process, sep="_")

	os.data.save(edgePairs,    outputDirectory, edgeFile , format="JSON")
	os.data.save(node1_counts, outputDirectory, geneDegreeFile , format="JSON")
	os.data.save(node2_counts, outputDirectory, ptDegreeFile, format="JSON")
		
	return( c(edgeFile, geneDegreeFile, ptDegreeFile) )

>>>>>>> develop
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
}
#----------------------------------------------------------------------------------------------------
get.edgePairs <- function(collection, genesetName, ...){				
  
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
    goi <- getGeneSet(genesetName)
 
    mtx <- convert.to.mtx(collection)
    
    ## get and save edge pairs
    edgePairs <- get.network_edges(t(mtx), samples=NA, genes=goi, ...)
<<<<<<< HEAD
=======
    goi <- genesets[[genesetName]]
 
    File <- collection$file
    dataObj <- fromJSON(paste(collection$directory, File,".json", sep="")) 
    mtx <- dataObj$data[[1]]
    rownames(mtx) <- dataObj$rows[[1]]
    colnames(mtx) <- dataObj$cols[[1]]
    
    ## get and save edge pairs
    edgePairs <- get.network_edges(mtx, samples=NA, genes=goi, ...)
>>>>>>> develop
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de

  return(edgePairs)
}

#----------------------------------------------------------------------------------------------------
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
run.batch.network_edges <- function(datasets){

  cat("-calculating edges\n")

    dataType <- "network"

    origin <- lapply(datasets , function(record){
      c(dataset = record$dataset, source = record$source)
    })
    origin <- unique(origin)
    # get unique dataset & source types
    
    # Loop for each dataset/source type, get mut &/or cnv edges
    for (collection in origin){
            
		  mut01_colls <- sapply(datasets, function(record){record$dataset== collection[["dataset"]] &
		                                                   record$source == collection[["source"]] &
		                                                   record$dataType=="mut01"})
		  mut01_colls <- datasets[mut01_colls]
      cnv_colls <- sapply(datasets, function(record){record$dataset== collection[["dataset"]] &
                                                       record$source == collection[["source"]] &
                                                       record$dataType=="cnv"})
      cnv_colls <- datasets[cnv_colls]
      
		  if(length(mut01_colls)==0 & length(cnv_colls) ==0) next;
		  cat(collection[["dataset"]], "\n")		  
		  
      for(geneset in genesets){
        EdgeList_mut <- EdgeList_cnv <- list()
        parent <- list()
        process <- list(geneset=geneset$name)
        
        prev.run <- collection.exists(mongo, dataset=collection[["dataset"]], dataType="edges",source=collection[["source"]],
                                      processName=paste(geneset$name, "mut01-cnv", sep="-"))
        if(prev.run){
          print("Skipping.")
          next()
        }
        
        
        if(length(mut01_colls)>0){
          mut01_coll <- mut01_colls[[1]]
          coll <- mongo.find.all(mongo, paste("oncoscape",mut01_coll$collection, sep="."))
          EdgeList_mut <- get.edgePairs(coll, geneset$name, edgeTypes=list("0"="1"))
          parent <- list(mut01_coll$`_id`)
          process$edgeType <- "mut01"
        }
        if(length(cnv_colls)>0){
          cnv_coll <- cnv_colls[[1]]
          coll <- mongo.find.all(mongo, paste("oncoscape",cnv_coll$collection, sep="."))
          EdgeList_cnv <- get.edgePairs(coll, geneset$name, edgeTypes=list("-2"="-2", "-1"="-1", "1"="1", "2"="2"))
          parent <- c(parent, cnv_coll$`_id`)
          process$edgeType <- c(process$edgeType, "cnv")
        }
		    processName=paste(unlist(process), collapse="-")
		    newEdges <- c(EdgeList_mut, EdgeList_cnv)
		        
		    save.edge.files(dataset=collection[["dataset"]], result=newEdges, source=collection[["source"]],
		                    parent=parent, process=process,processName=processName)				  
            
		}# for genesetName
  } #collection dataset/source type

<<<<<<< HEAD
=======
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

					  index <- get.new.collection.index(Manifest, datasetName, dataType)
					  parent <- list(c(datasetName, "cnv", collection$id))
					  process <- list(edgeType="cnv", geneset= genesetName); processName=paste(process, collapse="-")
					  edgeFiles <- save.edge.files(newEdges,output_directory, datasetName, dataType,index, processName)
					  newCollection <- data.frame(id=index,date=date,directory=output_directory)
					  newCollection$parent <- parent
					  newCollection$process <- list(process)
					  newCollection$file <- list(edgeFiles)
					  Manifest <- add.new.collection(Manifest, datasetName, dataType, newCollection)
					  
					  EdgeList$cnv[[as.character(index)]] <- newEdges
				  }
				}
				# save individual MUT edges				
				if(nrow(mutTables)==0)next;

				for(i in 1:nrow(mutTables)){
				  collection <- mutTables[i, "collections"][[1]]
				  newEdges <- get.edgePairs(collection, genesetName, edgeTypes=list("0"="1"))
				  
				  index <- get.new.collection.index(Manifest, datasetName, dataType)
				  parent <- list(c(datasetName, "mut", collection$id))
				  process <- list(edgeType="mut01", geneset= genesetName); processName=paste(process, collapse="-")
				  edgeFiles <- save.edge.files(newEdges,output_directory, datasetName, dataType,index, processName)
				  newCollection <- data.frame(id=index, date=date,directory=output_directory)
				  newCollection$parent <- parent
				  newCollection$process <- list(process)
				  newCollection$file <- list(edgeFiles)
				  Manifest <- add.new.collection(Manifest, datasetName, dataType, newCollection)
				  
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
						#compIndex = (k+m)+((k-1)*numCNV)+m
						index <- get.new.collection.index(Manifest, datasetName, dataType)
						edgeFiles <- save.edge.files(allEdges, output_directory, datasetName, dataType, index=index, processName)
						parent <- list(c(datasetName, "network", names(EdgeList$cnv)[k]),c(datasetName, "network", names(EdgeList$mut)[m]) )
						collection <- data.frame(id=index, date=date,directory=output_directory)
						collection$parent <- list(parent)
						collection$process <- list(process)
						collection$file <- list(edgeFiles)
						Manifest <- add.new.collection(Manifest, datasetName, dataType, collection)
					}
				}
				
				
			} # for genesetName
			
 		} # for diseaseName	

      return(Manifest)
>>>>>>> develop
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
}


#----------------------------------------------------------------------------------------------------
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
## must first initialize server (through shell >mongod)
mongo <- connect.to.mongo()

genesets <-     mongo.find.all(mongo, "oncoscape.hg19_genesets_hgnc_import", query=list())

molecular_manifest <- mongo.find.all(mongo, "oncoscape.manifest", 
                                    query='{"dataType":{"$in":["cnv","mut01", "rna", "protein", "methylation"]}}')
run.batch.patient_similarity(molecular_manifest, scaleFactor=10000)
		# calculate patient similarity

#molecular_manifest <- mongo.find.all(mongo, "oncoscape.manifest", 
 #                                    query='{"dataType":{"$in":["cnv","mut01"]}}')
#run.batch.network_edges(molecular_manifest)
		# map edges for all patients between CNV/Mut and Geneset tables


<<<<<<< HEAD
close.mongo(mongo)
=======

if("mds" %in% commands){
	Manifest_mds <- run.batch.patient_similarity(molecular_manifest,geneset_name="oncoVogel274", output_directory = mds_output_directory)
		# calculate patient similarity
		# save json
	os.data.save(
	  df = Manifest_mds,
	  directory="../manifests",
	  file= paste("os","mds", date,"manifest", sep="."),
	  format = "JSON") 
}

if("edges" %in% commands){
	Manifest_edges <- run.batch.network_edges(molecular_manifest, output_directory=network_output_directory)
		# map edges for all patients between CNV/Mut and Geneset tables
	os.data.save(
	  df = Manifest_edges,
	  directory="../manifests",
	  file= paste("os","edges", date,"manifest", sep="."),
	  format = "JSON") 
	
}


>>>>>>> develop
=======
close.mongo(mongo)
>>>>>>> 28fa7a153422d214caffae02005075ac91f165de
