###
#
#       This Script Executes Basic Processing On TCGA Files
#       Specifically It Types, Uppercases and In Cases Enforces Enumeration Types
#       
###

# Library Imports ---------------------------------------------------------
library(RUnit)
library(R.utils)
library(stringr)
library(plyr)
library(jsonlite)

# Configuration -----------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

os.tcga.batch.inputFile    <- "os.tcga.ucsc.filename.manifest.json"

# IO Utility Functions :: [Batch, Load, Save]  -------------------------------------------------------

### Save Function Takes A matrix/data.frame + Base File Path (w/o extension) & Writes to Disk In Multiple (optionally specified) Formats
os.data.save <- function(df, variable, directory, file, format = c("tsv", "csv", "RData", "manifest"), metaData=NA){
  
  outFile = paste(directory, file, sep="")
  
  # Write Tab Delimited
  if("tsv" %in% format)
    write.table(df, file=paste(outFile,".tsv", sep = ""), quote=F, sep="\t")
  
  # Write CSV Delimited
  if("csv" %in% format)
    write.csv(df, file=paste(outFile,".csv",sep = ""), quote = F)
  
  # Write RData File
  if("RData" %in% format)
    save(list=variable, file=paste(outFile,".RData", sep = "") )

  # Update Manifest File
  if("manifest" %in% format)
    write.table(metaData$data, file=paste(directory, "manifest.txt", sep = "", append=TRUE) )
  
  # Write JSON File
  if(!missing(metaData))
    write(toJSON(metaData, pretty=TRUE), file=paste(outFile,"_metadata.json", sep = "") )
  
  # Return DataFrame For Chaining
  return(df)
}

### Load Function Takes An Import File + Column List & Returns A DataFrame
os.data.load <- function(inputFile, variable, provenance){
  
  mtx <- read.delim(inputFile)
  category = ""; subcategory = "";
  entity.type = "gene"; feature.type = "patient"
  entity.count = 0; feature.count = 0;

  if(grep("cn", variable){ 
  	category = "copy number"; subcategory = "gistic scores";
  else if(grep("mut", variable){ category = "mutations"; subcategory = "indicator"
  
  return(df=mtx, metadata=list("variable"=variable, "class"="matrix", "category"=category,
         "subcategory"=subcategory, "entity.count"=nrow(mtx), "feature.count"=ncol(mtx),
         "entity.type"=entity.type, "feature.type" = feature.type, "minValue"=min(mtx), "maxValue"=max(mtx), 
         "provenance"= provenance))

 
}

### Batch Is Used To Process Multiple TCGA Files Defined 
os.data.batch <- function(inputFile, ...){
        
    # Load Input File 
    dataFiles <- fromJSON(inputFile)
                        
		# Loop for each disease type
		for (diseaseName in names(dataFiles))
		{
		  currentDirectory <- dataFiles[[diseaseName]][["directory"]]
		  currentFiles <- dataFiles[[diseaseName]][["files"]]
		  outputDirectory <- dataFiles[[diseaseName]][["output.directory"]]
		  
		    # Loop Column Wise: for each file type
				for (currentTable in names(currentFiles))
				{
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
}

# Run Block  -------------------------------------------------------
os.data.batch(
  inputFile = os.data.batch.inputFile,
  checkEnumerations = TRUE,
  checkClassType = "os.class.tcgaCharacter")
