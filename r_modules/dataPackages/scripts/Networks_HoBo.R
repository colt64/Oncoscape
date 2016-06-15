library(org.Hs.eg.db)
library(jsonlite)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

#--------------------------------- make plot data -----------------------------#

geneset_file <- "../molecular_data/hg19/genesets_by_symbol.json"
genesets <- fromJSON(geneset_file)
#----------------------------------------------------------------------------------------------------
getGeneSet <- function(geneset_name){
	stopifnot(geneset_name %in% names(genesets))
	return(genesets[[geneset_name]])
}

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
run.batch.patient_similarity <- function(input_directory, collection_metadata_file, update_collection_table=F,geneset_name=NA, output_directory="./"){

		if(update_collection_table)
			write_collections_table(input_directory, outputFile=collection_metadata_file)

		# table listing source, disease, process, date, molecular_type for each table (rowname=file)
		mol_metatable <- read.table(collection_metadata_file, sep="\t", header=T)

		PtLayouts <- list()
		gistic.scores <-c(-2,-1,1, 2)
  		goi = getGeneSet(geneset_name)

		diseases <- unique(mol_metatable$disease)
		for(diseaseName in diseases){
		  cat(diseaseName)
		  	PtLayouts[[diseaseName]] <- c()
			cnvTables <- subset(mol_metatable, disease==diseaseName & molecular_type=="cnv")
			mutTables <- subset(mol_metatable, disease==diseaseName & molecular_type=="mutation_01")
			
			if(nrow(cnvTables)==0 | nrow(mutTables) ==0) next;
			
			regex = ".01$"; threshold = NA;
			if(diseaseName == "laml"){       regex = "-03$|-09$";
			} else if(diseaseName == "luad"){ regex = "TCGA-(17)^-\\d{4}-01$" }

			if(diseaseName == "brca" | diseaseName == "brain")  threshold = -1e-04
			
			for(i in 1:nrow(cnvTables)){
			  cnvFile <- cnvTables[i,"file"]
				cnv.json <- fromJSON(paste(input_directory, cnvFile, sep="/")) 
				mtx.cnv <- cnv.json$data[[1]]
				rownames(mtx.cnv) <- cnv.json$rows[[1]]
				colnames(mtx.cnv) <- cnv.json$cols[[1]]
				cnv.samples <- grep(regex, cnv.json$cols[[1]],  value=TRUE)
				cnv.index <- cnvTables[i,"index"]

				for(j in 1:nrow(mutTables)){
				  mutFile <- mutTables[j,"file"]
					mut.json <- fromJSON(paste(input_directory, mutFile, sep="/")) 
					mtx.mut <- mut.json$data[[1]]
					rownames(mtx.mut) <- mut.json$rows[[1]]
					colnames(mtx.mut) <- mut.json$cols[[1]]
					mut.index <- mutTables[j,"index"]
   	
   				mut.samples <- grep(regex, mut.json$cols[[1]],  value=TRUE)
   					
   				samples <- unique(cnv.samples, mut.samples)
   
   				sample_similarity <- calculateSampleSimilarityMatrix(mtx.mut, mtx.cnv, copyNumberValues=gistic.scores,
   													 genes = goi, samples=samples)
   													 
					if(!is.na(threshold)){
						outliers <- names(which(sample_similarity[,1]<threshold))
						sample_similarity <- sample_similarity[setdiff(rownames(sample_similarity), outliers), ]
					}
   					
   					outputFile <- paste("mds",diseaseName, cnvTables$molecular_type, cnv.index, mutTables$molecular_type, mut.index, geneset_name, sep="_")
   					tempList<- lapply(rownames(sample_similarity), function(name) c(x=sample_similarity[name,"x"], y=sample_similarity[name, "y"]))
   					names(tempList) <- rownames(sample_similarity)
   					save.json(tempList, output_directory, outputFile)

					PtLayouts[[diseaseName]] <- c(PtLayouts[[diseaseName]], outputFile)
				} # mut files
			} #cnv files
		} # for diseaseName	

		save.json(PtLayouts, output_directory, "mds_collections")
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
  
  return(allEdges)
}
#----------------------------------------------------------------------------------------------------
run.batch.network_edges <- function(input_directory, collection_metadata_file, genesets, output_directory="./"){

		mol_metatable <- read.table(collection_metadata_file, sep="\t", header=T)
		

		diseases <- unique(mol_metatable$disease)
		for(diseaseName in diseases){
		  cat(diseaseName)
			cnvTables <- subset(mol_metatable, disease==diseaseName & molecular_type=="cnv")
			mutTables <- subset(mol_metatable, disease==diseaseName & molecular_type=="mutation_01")
			
			
			if(nrow(cnvTables)==0 & nrow(mutTables) ==0) next;
			for(genesetName in names(genesets)){			

				goi <- genesets[[genesetName]]
				edgeTypes <- list("-2"="-2", "-1"="-1", "1"="1", "2"="2")
				
				if(nrow(cnvTables) >0){
				for(i in 1:nrow(cnvTables)){
				  cnvFile <- cnvTables[i,"file"]
					cnv.json <- fromJSON(paste(input_directory, cnvFile, sep="/")) 
					mtx.cnv <- cnv.json$data[[1]]
					rownames(mtx.cnv) <- cnv.json$rows[[1]]
					colnames(mtx.cnv) <- cnv.json$cols[[1]]
					cnv.index <- cnvTables[i,"index"]
					
					edgePairs <- get.network_edges(mtx.cnv, samples=NA, genes = goi, edgeTypes)
					outputFile <- paste("edges",diseaseName, cnvTables$molecular_type, cnv.index, genesetName, sep="_")
					save.json(edgePairs, output_directory, outputFile)
				}}
				
				if(nrow(mutTables)==0)next;
				for(j in 1:nrow(mutTables)){
					mutFile <- mutTables[j,"file"]
					mut.json <- fromJSON(paste(input_directory, mutFile, sep="/")) 
					mtx.mut <- mut.json$data[[1]]
					rownames(mtx.mut) <- mut.json$rows[[1]]
					colnames(mtx.mut) <- mut.json$cols[[1]]
					mut.index <- mutTables[j,"index"]
					
					edgePairs <- get.network_edges(mtx.mut, samples=NA, genes = goi, edgeTypes=list("0"="1"))
					outputFile <- paste("edges",diseaseName, mutTables$molecular_type, mut.index, genesetName, sep="_")
					save.json(edgePairs, output_directory, outputFile)
					
				}
				

				
			} # for genesetName
 		} # for diseaseName	

}


#----------------------------------------------------------------------------------------------------
run.batch.network_collections <- function(input_directory, collection_metadata_file, genesets, output_directory="./"){

		mol_metatable <- read.table(collection_metadata_file, sep="\t", header=T)
		EdgeSets <- list()
			EdgeSets[[diseaseName]]
				# merge cnv/mut
				# get.node.degree(poi, goi, edgePairs)
				#


#----------------------------------------------------------------------------------------------------
# reads all files in a directory: assumes files follow json schema {disease, molecular_type, source, process, date}
# saves table with rownames = file names
write_collections_table <- function(directory,outputFile){
		
		jsonFiles<- list.files(directory)
		mol_metatable<- sapply(1:length(jsonFiles), function(i){
			jFile <- fromJSON(paste(directory,jsonFiles[i], sep="/"))
			c(index=i,disease=jFile$disease, molecular_type=jFile$molecular_type, source=jFile$source,
			  process=jFile$process,date=jFile$date, file=jsonFiles[i])
		})
		mol_metatable <- t(mol_metatable)
		
		write.table(mol_metatable, file=outputFile, sep="\t", quote=F, col.names=T, row.names=F )
}


#----------------------------------------------------------------------------------------------------
directory <- "../molecular_data/UCSC"
collection_metadata_file <- "molecular_collections_metadata.txt"
output_directory <- "../networks"

	
#run.batch.patient_similarity(directory, collection_metadata_file, update_collection_table=F,geneset_name="oncoVogel274", output_directory=output_directory)
		# calculate patient similarity
		# save json

run.batch.network_edges(directory, collection_metadata_file, genesets, output_directory=output_directory)
		# map edges for all patients between CNV/Mut and Geneset tables
