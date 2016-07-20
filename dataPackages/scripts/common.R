# Library Imports ---------------------------------------------------------
library(RUnit)
library(R.utils)
library(stringr)
library(plyr)
library(jsonlite)
library(rmongodb)


os.dataset.enumerations     <- fromJSON("../manifests/os.dataset.enumerations.json" )

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
	  mongo.drop(mongo, icoll)
	  #mongo.drop.database(mongo, db)
	  res <- mongo.get.database.collections(mongo, db)
	  print(res)

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
write.to.mongo <- function(mongo, db, dataObj){

	bson.data <- mongo.bson.from.JSON(toJSON(dataObj))

#a <- mongo.bson.from.JSON( '{"ident":"a", "name":"Markus", "age":33}' )
#b <- mongo.bson.from.JSON( '{"ident":"b", "name":"MongoSoup", "age":1}' )
#c <- mongo.bson.from.JSON( '{"ident":"c", "name":"UseR", "age":18}' )

#	icoll <- paste(db, "test", sep=".")
	mongo.insert.batch(mongo, db, bson.data )

#mongo.update(mongo, icoll, list('ident' = 'b'), list('$inc' = list('age' = 3)))
# mongo.index.create(mongo, icoll, list('ident' = 1))
 
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
save.mut01.from.mut <- function(Manifest, result, dataset, dataType,source, outputDirectory){
  
  resultObj <- data.frame(dataset = dataset, dataType = "mut01",
                          rowType= result$rowType, colType = result$colType)
  resultObj$rows <-  list(result$rows)
  resultObj$cols <-  list(result$cols);
  
  mtx.01 <- result$data[[1]]
  mtx.01[is.na(mtx.01)] <- 0
  mtx.01[mtx.01 == ""] <- 0
  mtx.01[nchar(mtx.01) >1] <- 1
  mtx.01 <- apply(mtx.01, 2, as.integer)
  
  resultObj$data <- list(mtx.01)
  
  parent <- list(c(dataset, dataType, result$id))
  
  Manifest <- save.collection(Manifest=Manifest, dataset=dataset, dataType="mut01",source=source, result=resultObj,
                              parent=parent, process=process,processName=process, outputDirectory=outputDirectory)
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
add.new.collection <- function(mongo, datasetName, dataTypeName, collection){
  
  db <- mongo.get.database.collections(mongo, "oncoscape")
  
  if(length(db) == 0){	
    newCollection <- data.frame(dataset=datasetName, dataType=dataTypeName)
    newCollection$collections <- list(collection)
    mongo.insert(mongo, "Manifest", newCollection)
  }
  
  Manifest <- mongo.find.all(mongo, "Manifest", query = list('dataset' = datasetName, 'dataType' = dataTypeName))
  mongo.insert.batch(mongo, Manifest, newCollection)
  
  if(nrow(Manifest) == 0){	
    newCollection <- data.frame(dataset=datasetName, dataType=dataTypeName)
    newCollection$collections <- list(collection)
    Manifest <- newCollection
    return()
  }
  
  dataObj <- subset(Manifest, dataset == datasetName & dataType == dataTypeName)
  if(nrow(dataObj) == 1){
    newCollection <- list(rbind(dataObj$collections[[1]],collection))
    Manifest[Manifest$dataset==datasetName & Manifest$dataType ==dataTypeName,"collections"] <- list(newCollection)
    return()
  }
  if(nrow(dataObj) == 0){	
    newCollection <- data.frame(dataset=datasetName, dataType=dataTypeName)
    newCollection$collections <- list(collection)
    Manifest <- rbind(Manifest, newCollection)
    return()
  }
  stop(printf("add.new.collection found %d instances of dataset %s and dataType %s", length(dataObj), datasetName, dataTypeName))
  
}
#---------------------------------------------------------
save.collection <- function(mongo, dataset, dataType,source,result, parent, 
							process,processName){

	cat("-save collection\n")

#  index <- get.new.collection.index(Manifest, dataset, dataType)
#  result$id <- index
#  outputFile <- paste(dataset, dataType, index, processName , sep="_")
  
  source <- unique(source)
  if(length(source)>1) source <- list(source)
 
  newCollection <- data.frame(date=date) 
#  newCollection <- data.frame(id=index,
#                              date = date,
 #                             directory= outputDirectory,
#                              file=outputFile)
  newCollection$source <- source
  newCollection$process <- process
  newCollection$parent <- parent
  add.new.collection(mongo, dataset, dataType, newCollection)
  
  # Save Data Frame
#  os.data.save( mongo,
#    df = result,
#    file= outputFile,
#    ) 
    
  if(dataType == "mut")
  	save.mut01.from.mut(result, dataset, dataType,source)
    
	return()
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
  