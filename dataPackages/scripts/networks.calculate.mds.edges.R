library(org.Hs.eg.db)
library(jsonlite)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

commands >- c("mds", "edges")
args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

#--------------------------------- Configuration -----------------------------#

molecular_manifest <- "../manifests/os.ucsc.molecular.manifest.json"

network_output_directory <- "../data/molecular/edges"
mds_output_directory <- "../data/molecular/mds"

geneset_file <- "../molecular_data/hg19/genesets_by_symbol.json"
genesets <- fromJSON(geneset_file)

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
run.batch.patient_similarity <- function(input_directory, manifest_file, geneset_name=NA, output_directory="./"){


    # Load Input File 
    dataFiles <- fromJSON(inputFile)
 
 	PtLayouts <- list()
	gistic.scores <-c(-2,-1,1, 2)
  	goi = getGeneSet(geneset_name)

     
		# Loop for each disease/source combo
		for (i in 1:nrow(dataFiles))
		{
 			# Necessary Fields
			sourceObj <- dataFiles[i,]
			stopifnot(all(c("disease", "source", "collections") %in% names(sourceObj)))

			collections <- sourceObj$collections[[1]]
			diseaseName <- sourceObj$disease
			inputDirectory <- sourceObj$inputDirectory
			if(!grepl("/$", inputDirectory)) inputDirectory <- paste(inputDirectory, "/", sep="")
			cat(diseaseName, sourceObj$source,"\n")

#			outputDirectory <- sourceObj$outputDirectory

			## ----- Configuration ------
			regex = ".01$"; threshold = NA;
			if(diseaseName == "laml"){        regex = "-03$|-09$";
			} else if(diseaseName == "luad"){ regex = "TCGA-(17)^-\\d{4}-01$" }

			if(diseaseName == "brca" | diseaseName == "brain")  threshold = -1e-04

		  	PtLayouts[[diseaseName]] <- c()
			cnvTables <- subset(collections, molecular_type=="cnv")
			mutTables <- subset(collections, molecular_type=="mutation01")

			if(nrow(cnvTables)==0 | nrow(mutTables) ==0) next;
			
			## ----- MDS on All Combinations of CNV and MUT Tables ------
			for(i in 1:nrow(cnvTables)){
				cnv.json <- fromJSON(paste(inputDirectory, cnvTables[i,"inputFile"], sep="/")) 
				mtx.cnv <- cnv.json$data[[1]]
				rownames(mtx.cnv) <- cnv.json$rows[[1]]; colnames(mtx.cnv) <- cnv.json$cols[[1]]
				cnv.samples <- grep(regex, cnv.json$cols[[1]],  value=TRUE)
				cnv.index <- cnvTables[i,"index"]

				for(j in 1:nrow(mutTables)){
					mut.json <- fromJSON(paste(inputDirectory,mutTables[j,"inputFile"], sep="/")) 
					mtx.mut <- mut.json$data[[1]]
					rownames(mtx.mut) <- mut.json$rows[[1]]; colnames(mtx.mut) <- mut.json$cols[[1]]
   					mut.samples <- grep(regex, mut.json$cols[[1]],  value=TRUE)
					mut.index <- mutTables[j,"index"]
   					
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
					outputFile <- paste("mds",diseaseName, cnvTables$dataType, cnv.index, mutTables$dataType, mut.index, geneset_name, sep="_")
					mds.list<- lapply(rownames(sample_similarity), function(name) c(x=sample_similarity[name,"x"], y=sample_similarity[name, "y"]))
					names(mds.list) <- rownames(sample_similarity)
					os.data.save(mds.list, outputDirectory, outputFile, format="JSON")

					PtLayouts[[diseaseName]] <- c(PtLayouts[[diseaseName]], outputFile)
					
				} # mut files
			} #cnv files
		} # for diseaseName	

		os.data.save(PtLayouts, output_directory, "mds_collections", format="JSON")
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
save.edge.files(edgePairs, outputDirectory, diseaseName, dataType, index, genesetName)

	## get and save node degrees
	node1_counts <- as.data.frame(table(edgePairs[,2]))  
	colnames(node1_counts) <- NULL

	node2_counts <- as.data.frame(table(edgePairs[,3]))  
	colnames(node2_counts) <- NULL
	
	edgeFile <- paste("edges"     ,diseaseName, dataType, index, genesetName, sep="_")
	geneDegreeFile <- paste("geneDegree",diseaseName, dataType, index, genesetName, sep="_")
	ptDegreeFile <- paste("ptDegree"  ,diseaseName, dataType, index, genesetName, sep="_")

	os.data.save(edgePairs,    outputDirectory, edgeFile , format="JSON")
	os.data.save(node1_counts, outputDirectory, geneDegreeFile , format="JSON")
	os.data.save(node2_counts, outputDirectory, ptDegreeFile, format="JSON")
		
	return( c(edgeFile, geneDegreeFile, ptDegreeFile) )

}
#----------------------------------------------------------------------------------------------------
get.edges.batch <- function(sourceObj, collections, genesetName, startIndex=1, ...){				

		inputDirectory  <- sourceObj$inputDirectory
		outputDirectory <- sourceObj$outputDirectory
		goi <- genesets[[genesetName]]

		edges <- list()
		newCollections <- list()
		
		for(i in 1:nrow(collections)){
		  	File <- collections[i,"inputFile"]
			dataObj <- fromJSON(paste(inputDirectory, File, sep="/")) 
			mtx <- dataObj$data[[1]]
			rownames(mtx) <- dataObj$rows[[1]]
			colnames(mtx) <- dataObj$cols[[1]]
			index <- startIndex+i-1
	
			## get and save edge pairs
			edgePairs <- get.network_edges(mtx, samples=NA, ...)
			edgeFiles <- save.edge.files(edgePairs, outputDirectory, sourceObj$diseaseName, collections$dataType, index, genesetName)
			newCollections <- c(newCollections, data.frame(_id=index, process=process, date=date, inputFile=cFile, outputFile=edgeFiles))

			edges[[as.character(index)]] <-  edgePairs
		}
		
		return(list(edges=edges, collections=newCollections))
}

#----------------------------------------------------------------------------------------------------
run.batch.network_edges <- function(manifest_file, output_directory="./"){

    # Load Input File 
    dataFiles <- fromJSON(manifest_file)
    dataType <- "network"
    date <- System.Date(format="%Y-%m-%d")
    
    newManifest <- list()
     
		# Loop for each disease/source combo
		for (i in 1:nrow(dataFiles))
		{
			sourceObj <- dataFiles[i,]
			stopifnot(all(c("disease", "source","dataType", "collections") %in% names(sourceObj)))

			collections <- sourceObj$collections[[1]]
			diseaseName <- sourceObj$disease
			source      <- sourceObj$source
			
			sourceTokens <- c(diseaseName, source)
			cat(diseaseName, source,"\n")
		  
			cnvTables <- subset(collections, molecular_type=="cnv")
			mutTables <- subset(collections, molecular_type=="mutation01")
			
			if(nrow(cnvTables)==0 & nrow(mutTables) ==0) next;

			newCollections <- list()

			# save Edge sets for each Gene set
			for(index in 1:length(genesets)){			
	
				genesetName <- names(genesets)[index]

				EdgeList <- list()
			  	process <- c("cnv;mut01", genesetName)
			  	processName <- paste(process, collapse=";")
			  	 
				# save individual CNV edges
				if(nrow(cnvTables) >0){
					newEdges <- get.edges.batch(sourceObj, cnvTables, genesetName, edgeTypes=list("-2"="-2", "-1"="-1", "1"="1", "2"="2"))
					EdgeList$cnv <- newEdges$edges
					newCollections <- c(newCollections, newEdges$collections)
				}
				# save individual MUT edges				
				if(nrow(mutTables)==0)next;
				newEdges <- get.edges.batch(sourceObj, mutTables, genesetName, edgeTypes=list("0"="1"), startIndex=nrow(cnvTables)+1)
				EdgeList$mut <- newEdges$edges
				newCollections <- c(newCollections, newEdges$collections)

				# save all combinations of MUT & CNV edges
				numCNV = length(EdgeList$cnv)
				numMut = length(EdgeList$mut)
				for(k in 1:numCNV){
					cnvEdges <- EdgeList$cnv[[k]]
					
					for(m in 1:numMut){
						mutEdges <- EdgeList$mut[[m]]
																	
						allEdges <- rbind(cnvEdges, mutEdges)
						compIndex = (k+m)+((k-1)*numCNV)+m
						edgeFiles <- save.edge.files(allEdges, outputDirectory, diseaseName, dataType, index=compIndex, processName)
						file1 <- newCollections[[k]]$inputFile; file2 <- newCollections[[m]]$inputFile;
						newCollections <- c(newCollections, data.frame(_id=compIndex, process=process, date=date, inputFile=c(file1,file2), outputFile=edgeFiles))
					}
				}
				
				
			} # for genesetName
			
			newManifest <- c(newManifest, 
				data.frame(disease=diseaseName, source=sourceObj$source, dataType=dataType, inputDirectory=sourceObj$outputDirectory, 
							collections = newCollections) )

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
