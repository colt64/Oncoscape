# Library Imports ---------------------------------------------------------
library(RUnit)
library(R.utils)
library(stringr)
library(plyr)
library(jsonlite)
library(rmongodb)


os.dataset.enumerations     <- fromJSON("../manifests/os.dataset.enumerations.json" )
date <- as.character(Sys.Date())
chromosomes <- c(seq(1:22), "X", "Y")

#---------------------------------------------------------
connect.to.mongo <- function(host= "127.0.0.1", name = "", username = "", password = "", db = "admin"){
	mongo <- mongo.create(host = host, name = name, username = username,
  							password = password, db = db, timeout = 0L)
	
	stopifnot(mongo.is.connected(mongo))
	return(mongo)
}

#---------------------------------------------------------
close.mongo <- function(mongo){

	if(mongo.is.connected(mongo) == TRUE) {
	  # close connection
	  mongo.destroy(mongo)
	}
}

#---------------------------------------------------------
get.mongo.table <- function(mongo, db, table){

	if(mongo.is.connected(mongo)){
		collections <- mongo.get.database.collections(mongo, db)
		pop <- mongo.distinct(mongo, coll, "pop")
		pops1 <- mongo.find.all(mongo, coll, query = list('pop' = list('$lte' = 2), 'pop' = list('$gte' = 1)))
	}
}

#---------------------------------------------------------
mongo.collection.as.matrix <- function(collection, format=""){

  cursor <- mongo.find(mongo, paste("oncoscape",collection, sep="."), query=list(), fields=list())
  count <- mong.count(mongo, paste("oncoscape",collection, sep="."))
  result_lst <- vector('list', count)
  i <- 1
  while (mongo.cursor.next(cursor)) {
    result_lst[[i]] <- mongo.bson.to.list(mongo.cursor.value(cursor))
    val <-geneRow$patients; 
    null.val <- which(unlist(lapply(val, is.null)))
    if(length(null.val)>0) val[null.val] <- NA
    val <- unlist(val);
    if(format == "as.numeric") val <- as.numeric(val)
    
    i <- i + 1
  }
  result_dt <- data.table::rbindlist(result_lst)

  colnames(mtx) <- sapply(data.list, function(geneRow){ geneRow$gene})
  rownames(mtx) <- names(data.list[[1]]$patients)
  return(mtx)  
}

#---------------------------------------------------------
convert.to.mtx <- function(data.list, format=""){
  mtx <- sapply(data.list, function(geneRow){ 
    val <-geneRow$patients; 
    null.val <- which(unlist(lapply(val, is.null)))
    if(length(null.val)>0) val[null.val] <- NA
    val <- unlist(val);
    if(format == "as.numeric") val <- as.numeric(val)
#    if(all(is.null(val))){ val <- rep(NA, length(geneRow$patients))} 
    val})
  colnames(mtx) <- sapply(data.list, function(geneRow){ geneRow$gene})
  rownames(mtx) <- names(data.list[[1]]$patients)
  return(mtx)  
}

#---------------------------------------------------------
mapProcess <- function(process){
  
	processFound <-	sapply(os.dataset.enumerations$dataType, function(typeMap){ process %in% unlist(typeMap) })
	numMatches <- length(which(processFound))
	if(numMatches==1)
		return (names(os.dataset.enumerations$dataType)[which(processFound)])

	stop(printf("mapProcess found %d matches for process %s", numMatches, process))
	return(NA)
}
#---------------------------------------------------------
### For any mutation file, create and save an indicator mut01 file
save.mut01.from.mut <- function(mongo, result, dataset, dataType,source, parentID){
  
  mut.list <- result
  
  data.list <- lapply(result, function(geneSet){
    patients <- lapply(geneSet$patients, function(pt){mut <- nchar(pt); mut01 <- ifelse(mut > 0, 1, 0 ); mut01 })
    list(gene=geneSet$gene,min=min(unlist(patients)), max=max(unlist(patients)), patients = patients)
  })    
  
  parent <- parentID
  
  save.collection(mongo, dataset=dataset, dataType="mut01",source=source, result=data.list,
                              parent=parent, process=process,processName=process)

}

#---------------------------------------------------------
collection.exists <- function(mongo, dataset, dataType,source,processName){

  source <- unique(source)
  if(length(source)>1) source <- list(source)
  sourceName <- paste(unlist(source), collapse="-")

  collection.uniqueName <- paste(dataset, dataType, sourceName, processName, sep="_")
  collection.ns <- paste("oncoscape", collection.uniqueName, sep=".")
  if(mongo.count(mongo, collection.ns) != 0){
    print(paste(collection.uniqueName, " already exists.", sep=""))
    return(TRUE)
  }  
  return(FALSE)

}
#---------------------------------------------------------
save.collection <- function(mongo, dataset, dataType,source,result, parent, 
                            process,processName){
  
  cat("-save collection\n")
  
  source <- unique(source)
  if(length(source)>1) source <- list(source)
  sourceName <- paste(unlist(source), collapse="-")
  
  collection.uniqueName <- paste(dataset, dataType, sourceName, processName, sep="_")
  collection.ns <- paste("oncoscape", collection.uniqueName, sep=".")
  if(mongo.count(mongo, collection.ns) != 0){
    print(paste(collection.uniqueName, " already exists. Skipping.", sep=""))
    return()
  }  
  
  newCollection <- list(dataset=dataset, dataType=dataType, date=date) 
  newCollection$collection <- collection.uniqueName
  newCollection$source <- source
  newCollection$process <- process
  newCollection$parent <- parent
  
  ## add to manifest file
  mongo.insert(mongo, "oncoscape.manifest", newCollection)
  
  pass <- lapply(result, function(el){mongo.insert(mongo, collection.ns, as.list(el))})
  if(!all(unlist(pass))){
    print("ERROR: result not inserted into mongodb.")
    return()
  }
  
  newID <-  mongo.find.one(mongo, "oncoscape.manifest", 
                           query=newCollection, fields=list("_id"))
  
  ## add to lookup table  
  lookup.ns <-  "oncoscape.lookup_oncoscape_datasources"
  query <- list("disease"=dataset)
  datasource <- mongo.find.one(mongo, lookup.ns, query)
  
  if(length(datasource)==0){
    data.list <- list(disease = dataset, source = "TCGA",beta = TRUE)
    data.list$name = dataset
    data.list$img = paste(dataset, "png", sep=".")
  }else{
    data.list <- mongo.bson.to.list(datasource)
  }
  
  if(dataType %in% c("cnv","mut01", "mut", "rna", "protein", "methylation")){
    # update molecular
    
    add.collection <- list(data.frame(source=source, type=dataType, collection=collection.uniqueName))
    if("molecular" %in% names(data.list)){
      data.list$molecular <- c(data.list$molecular, add.collection)
    }else{data.list$molecular <- add.collection}
    
  }else if(dataType %in% c("mds", "pcaScores")){
    #update calculated
    add.collection <- list(data.frame(source=source, type=dataType, collection=collection.uniqueName))
    if("calculated" %in% names(data.list)){
      data.list$calculated	<-c(data.list$calculated, add.collection)
    } else {data.list$calculated <- add.collection }
    
  }else if(dataType %in% c("edges")){
    #update edges
    add.collection <- list(data.frame(name=process$geneset,edges=collection.uniqueName, 
                           patientWeights=paste(dataset, "ptDegree", source, processName, sep="_"), 
                           genesWeights=paste(dataset, "geneDegree", source, processName, sep="_")))
    if("edges" %in% names(data.list)){
      data.list$edges	<- c(data.list$edges, add.collection)
    } else {data.list$edges <- add.collection }
    
  }else if(dataType %in% c("ptDegree", "geneDegree")){
    print(paste(dataType, "lookup info processed with edge creation", sep=" "))
    return()
  }else if(dataType %in% c("patient", "drug", "radiation", "followUp-v1p0","followUp-v1p5", "followUp-v2p1", "followUp-v4p0", "newTumor", "newTumor-followUp-v4p0", "otherMalignancy-v4p0")){
    #update patient
    add.collection <- list()
    add.collection[dataType] <- collection.uniqueName
    if("clinical" %in% names(data.list)){
      data.list$clinical	<- c(data.list$clinical, add.collection)
    } else {data.list$clinical <- add.collection }
    
  }else if(dataType %in% c("chromosome", "centromere", "genes")){
    #update patient
    add.collection <- list()
    add.collection[dataType] <- collection.uniqueName
    if("location" %in% names(data.list)){
      data.list$location	<- c(data.list$location, add.collection)
    } else {data.list$location <- add.collection }
    
  }else if(dataType %in% c("genesets")){
    add.collection <- list(data.frame(source=source, type=dataType, collection=collection.uniqueName))
    if("category" %in% names(data.list)){
      data.list$category <- c(data.list$category, add.collection)
    }else{data.list$category <- add.collection}
    
  }else{
    print(paste("WARNING: data type not recognized:", dataType, sep=" "))
    return()
  }
  
  ## insert lookup into mongo collection
  mongo.update(mongo, lookup.ns, query, data.list, mongo.update.upsert)
  
  if(dataType == "mut")
    save.mut01.from.mut(mongo, result, dataset, dataType,source, parentID=newID)
  
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
##----------------------------
## get parent collection given the parentID: [dataset, dataType, collection$id]
get.parent.collection <- function(manifest, parentID){
  
  
  if(is.list(parentID)){
    if(length(parentID)>1){
      collection <- lapply(parentID, function(parent){
        get.parent.collection(manifest, parent)
        #        subset(manifest, dataset==parent[1] & dataType==parent[2] & collections$id==parent[3])
      })
      return(collection)
    }
    if(length(parentID)==1)
      parentID <- unlist(parentID)
  }
  if(is.data.frame(parentID)){
    collection <- apply(parentID, 1, function(parent){ get.parent.collection(manifest, parent)})
  }
  if(is.null(parentID))return(NA)
  Obj <- subset(manifest, dataset==parentID[1] & dataType==parentID[2])$collections[[1]]
  collection <- subset(Obj, id==parentID[3])
  return(collection)
}
##----------------------------
## get all collections with a specific key:value pair within the process
subset.collections <- function(collections, process=NA){

    coll.subset <- collections
    if(!all(is.na(process))){
      matchColl <- sapply(coll.subset$process, function(coll){ 
        all(names(process) %in% names(coll)) &&
        all(unlist(sapply(names(process),function(param){
          length(intersect(unlist(coll[[param]]),process[[param]])) == length(process[[param]]) }) ) )
      })
      matchColl[is.na(matchColl)] <- FALSE
      coll.subset <- coll.subset[matchColl,]
    }
      
    return(coll.subset)
}
#---------------------------------------------------------
# Aggregate unmapped column names and classes into a single list  
appendList <- function (x, val) 
{
    if(!is.list(x) && !is.list(val)) return(x)
    xnames <- names(x)
    for (v in names(val)) {
        x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
            appendList(x[[v]], val[[v]])
        else unique(c(x[[v]], val[[v]]))
    }
    x
}

##----------------------------
os.copy.file <- function(inputDir, filename, outputDir= "./"){

    file.copy(paste(inputDir, filename, sep=""), paste(outputDir, filename, sep=""))
  
   # data <- fromJSON(paste(inputDir, filename, sep=""))
  #  os.data.save(data, outputDir, filename, format= "JSON")
}
  

#--------------------------------------------------------------#
get.chromosome.dimensions <- function(scaleFactor=100000){
  
  chrPosScaledObj <- mongo.find.all(mongo, "oncoscape.manifest", list(dataset="hg19",dataType="chromosome", process=list(scale=scaleFactor)))[[1]]
  
  chrCoord <- mongo.find.all(mongo, paste("oncoscape",chrPosScaledObj$collection, sep="."))[[1]][["data"]]
  chrPos_xy <-t(sapply(chromosomes, function(chr){ return(c(chrCoord[[chr]]$x, chrCoord[[chr]]$q))}))
  chrDim <- c(max(chrPos_xy[,1]), max(chrPos_xy[,2]))
  
  return(chrDim)
}
#--------------------------------------------------------------#
scaleSamplesToChromosomes <- function(mtx, chrDim, dim.names=c("x", "y", "z")){
  
  mtx <- apply(mtx, 2, function(col){ -1* min(col) + col})
  # offset mtx so min val is 0,0
  mtx.max <- apply(mtx, 2, max)
  
  r2Chr <- sum(chrDim*chrDim)
  r2Mtx <- sum(mtx.max*mtx.max)	
  scale <- sqrt(r2Chr/r2Mtx)
  # make diagonal of drawing regions equal
  
  mtx <- mtx * scale
  mtx <- round(mtx)
  
  list.coord <- lapply(rownames(mtx), function(name){
    vals <- data.frame(t(mtx[name,dim.names]))
  })
  names(list.coord) <- rownames(mtx)
  return(list.coord)
}
#--------------------------------------------------------------#
scaleGenesToChromosomes <- function(genePos, chrCoordinates, scaleFactor=1000){
  
  genePos_xy <- lapply(genePos, function(gene){
    x <- chrCoordinates[gene[1], "xOffset"]
    y <- chrCoordinates[gene[1], "yOffset"] + as.numeric(gene[2])/scaleFactor
    data.frame(x=round(x),y=round(y))
  })
  
  return(genePos_xy)	
  
}

#--------------------------------------------------------------#
save.batch.genesets.scaled.pos <- function(scaleFactor=100000){
  
  geneObj<- mongo.find.all(mongo, "oncoscape.manifest", list(dataset="hg19", dataType="genes", process=list(scale=scaleFactor)))[[1]]
  genePos_scaled <- mongo.find.all(mongo, paste("oncoscape",geneObj$collection, sep="."))[[1]]
  
  genesetObj <-  mongo.find.all(mongo, "oncoscape.manifest", list(dataset="hg19",dataType="genesets"))[[1]]
  genesets <- mongo.find.all(mongo, paste("oncoscape",genesetObj$collection, sep="."))
  
  process <- list(calculation="scaled", scaleFactor=scaleFactor); 
  processName <- paste(process, collapse="-")
  parent <- list(geneObj$`_id`,genesetObj$`_id`)
  
  result <- lapply(genesets, function(geneSet){	
    genes <- geneSet$genes
    map_genes <- intersect(genes, names(genePos_scaled$data))
    genesetPos <- genePos_scaled$data[map_genes]
    list(type="geneset", name=geneSet$name, scale=scaleFactor, data=genesetPos)
  }	)
  
  save.collection(mongo, dataset=geneObj$dataset, dataType="genesets",source=geneObj$source, result=result,
                  parent=parent, process=process,processName=processName)
}
