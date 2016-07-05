os.dataset.enumerations     <- fromJSON("../manifests/os.dataset.enumerations.json" )

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
  
}
#---------------------------------------------------------
save.collection <- function(Manifest, dataset, dataType,result, parent, 
							process,processName, outputDirectory){

  index <- get.new.collection.index(Manifest, dataset, dataType)
  result$id <- index
  outputFile <- paste(dataset, dataType, index, processName , sep="_")
  newCollection <- data.frame(id=index,
                              date = date,
                              directory= outputDirectory,
                              file=outputFile)
  newCollection$process <- process
  newCollection$parent <- parent
  Manifest <- add.new.collection(Manifest, dataset, dataType, newCollection)
  
  # Save Data Frame
  os.data.save(
    df = result,
    directory=outputDirectory,
    file= outputFile,
    format = "JSON") 
    
	return(Manifest)
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
subset.collections <- function(collections, process=NA){

    coll.subset <- collections
    if(!all(is.na(process))){
      matchColl <- apply(coll.subset,1, function(coll){ 
        all( apply(process,1, function(param){
          coll$process[param["field"]] == param["value"] }) )
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
  