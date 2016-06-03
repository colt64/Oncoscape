###
#
#       This Script Executes Basic Processing On TCGA Files
#       Specifically It Types, Uppercases and In Cases Enforces Enumeration Types
#       
###

# Library Imports ---------------------------------------------------------
library(jsonlite)

# Configuration -----------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

os.data.batch.inputFile    <- "os.tcga.ucsc.filename.manifest.json"

# IO Utility Functions :: [Batch, Load, Save]  -------------------------------------------------------

### Save Function Takes A matrix/data.frame + Base File Path (w/o extension) & Writes to Disk In Multiple (optionally specified) Formats
os.json.save <- function(df, directory, file){
  
  outFile = paste(directory, file, sep="")
  write(toJSON(df, pretty=TRUE), file=paste(outFile,".json", sep = "") )
  
  # Return DataFrame For Chaining
  return(df)
}

### Load Function Takes An Import File + Column List & Returns A DataFrame
os.data.load <- function(inputFile, variable, provenance){
  
  mtx<- read.delim(inputFile)
  
  noGeneSymbol <- which(is.na(mtx[,1]))
  if(length(noGeneSymbol) > 0)
	  mtx <- mtx[-noGeneSymbol,]
  rownames(mtx) <- mtx[,1]
  mtx <- as.matrix(mtx[,-1])
  
  category = ""; subcategory = "";
  entity.type = "gene"; feature.type = "patient"
  entity.count = nrow(mtx); feature.count = ncol(mtx);
  minVal = min(mtx); maxVal = max(mtx); 

  if(grepl("cn", variable)){ 
  	category = "copy number"; subcategory = "gistic scores"; 
  }else if(grepl("mut", variable)){ category = "mutations"; subcategory = "indicator" }
  
  return(list(df=mtx, metaData=list("file"=paste(variable, ".RData", sep=""), "variable"=variable, "class"="matrix", "category"=category,
         "subcategory"=subcategory, "entity.count"=entity.count, "feature.count"=feature.count,
         "entity.type"=entity.type, "feature.type" = feature.type, "minValue"=minVal, "maxValue"=maxVal, 
         "provenance"= provenance)))

 
}

### Batch Is Used To Process Multiple TCGA Files Defined 
os.data.batch <- function(inputFile){
        
    # Load Input File 
	#   - an array of objects each specifying an input file with metadata
    dataFiles <- fromJSON(inputFile)

                        
		# Loop for each disease type
		for (i in 1:length(dataFiles))
		{
			dataObj <- dataFiles[[i]]
		  currentDirectory <- dataFiles[[i]][["directory"]]
		  currentFile <- dataFiles[[i]][["input_file"]]
		  outputDirectory <- dataFiles[[i]][["output_dir"]]
  
		  cat(diseaseName, currentTable,"\n")
			inputFile <- paste(currentDirectory, currentFiles[[currentTable]], sep = "")
			outputFile <- paste(outputDirectory, currentTable, sep="")
			provenance <- currentFiles[[currentTable]]
			
			# Load Data Frame - map and filter by named columns
			result <- os.data.load( inputFile = inputFile, variable = currentTable, provenance = provenance )

			df <- result$df
			metaData <- list("directory"=currentDirectory, "file"=currentTable, "data"=result$metaData)
	
			# Save Data Frame
			os.data.save(
					df = df,
					variable = currentTable,
					directory = outputDirectory,
					file = currentTable,
					format = c("RData", "manifest"),
					metaData = metaData)
			
			# Remove Df From Memory
			rm(df)
        }
}

# Run Block  -------------------------------------------------------

os.data.batch(  inputFile = os.data.batch.inputFile)
