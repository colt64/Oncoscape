###
#
#       This Script Executes Basic Processing On TCGA molecular Files downloaded from UCSC
#       Specifically it takes 
#       
###

# Library Imports ---------------------------------------------------------
library(jsonlite)

# Configuration -----------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

os.data.batch.inputFile    <- "os.tcga.ucsc.filename.manifest.json"

# IO Utility Functions :: [Batch, Load, Save]  -------------------------------------------------------

### Save Function Takes A matrix/data.frame + Base File Path (w/o extension) & Writes to Disk in JSON
os.json.save <- function(df, directory, file){
  
  if(!dir.exists(directory))
    dir.create(file.path(directory), recursive=TRUE)
  
  if(!grepl("/$", directory)) directory <- paste(directory, "/", sep="")
  outFile = paste(directory, file, sep="")
  write(toJSON(df, pretty=TRUE), file=paste(outFile,".json", sep = "") )
  
  # Return DataFrame For Chaining
  return(df)
}

### Takes matrix, returns list of row names, col names, and data with NA labels removed from mtx
get.processed.mtx <- function(mtx, dimension){
  if("row" %in% dimension){
    noName_row <- which(is.na(mtx[1,]))
    if(length(noName_row) > 0)
      mtx <- mtx[-noName_row,]
    
    rownames <- mtx[-1,1]
    mtx <- as.matrix(mtx[,-1])
  } 
  if("col" %in% dimension){
    noName_col <- which(is.na(mtx[,1]))
    if(length(noName_col) > 0)
      mtx <- mtx[,-noName_col]
    
    colnames <- mtx[1,]
    mtx <- as.matrix(mtx[-1,])
  }
  
  return(list(rownames=rownames, colnames=colnames, data=mtx))
}


### Load Function Takes An Import File + Column List & Returns A DataFrame
os.data.load <- function(inputFile){
  
  mtx<- read.delim(inputFile, header=F)
  
  #orient mtx so row: gene, col: patient/sample
  if(all(grepl("^TCGA", mtx[-1,1]))) { mtx <- t(mtx)}
  col_type <- "patient"; row_type <- "gene"

  if(all(grepl("TCGA-\\w{2}-\\w{4}-\\w{2}", mtx[1,-1])))
    col_type <- "sample"

  mtx.Data<- get.processed.mtx(mtx, dimension= c("row", "col"))
 
  return(list(row_type=row_type, col_type=col_type, rows=mtx.Data$rownames, cols=mtx.Data$colnames, data=mtx.Data$data))
 
}

### Batch Is Used To Process Multiple TCGA Files Defined 
os.data.batch <- function(inputFile){
        
    # Load Input File 
	#   - an array of objects each specifying an input file with metadata
    dataFiles <- fromJSON(inputFile)

                        
		# Loop for each disease type
		for (i in 1:nrow(dataFiles))
		{
			dataObj <- dataFiles[i,]
		  currentDirectory <- dataObj$directory
		  currentFile <- dataObj$input_file
		  outputDirectory <- dataObj$output_dir
		  outputFile_tokens <- paste(dataObj$source, dataObj$disease, dataObj$molecular_type, 
		                             dataObj$process, dataObj$date, sep="_")
  
		  cat(dataObj$disease, dataObj$molecular_type,"\n")
			inputFile <- paste(currentDirectory, currentFile, sep = "")

			# Load Data Frame - map and filter by named columns
			result <- os.data.load( inputFile = inputFile )

			if(dataObj$molecular_type %in%  c("cnv", "mutation01"))
			    result$data <- apply(result$data, 2, as.integer)
			dataObj$directory <- NULL
			dataObj$row_type <- result$row_type; dataObj$col_type <- result$col_type;
			dataObj$rows <- list(result$rows); dataObj$cols <- list(result$cols);
			dataObj$data <- list(result$data)
#			fullData <- c(as.list(dataObj), result)
	
			# Save Data Frame
			os.json.save(
					df= dataObj,
					directory = outputDirectory,
					file = outputFile_tokens)

    }
}

# Run Block  -------------------------------------------------------

os.data.batch(  inputFile = os.data.batch.inputFile)
