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

#commands <- c("molecular", "clinical")
commands <- c("molecular")
args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

os.molecular.ucsc.batch   <- fromJSON("../manifests/os.ucsc.molecular.manifest.json")
os.clinical.tcga.batch    <- fromJSON("../manifests/os.tcga.clinical.manifest.json")

outputDir_molecular <- "../data/molecular/clean/"
outputDir_clinical  <- "../data/clinical/clean/"

os.tcga.field.enumerations  <- fromJSON(paste("../manifests","os.tcga.field.enumerations.json" , sep="/"))
os.tcga.column.enumerations <- fromJSON(paste("../manifests","os.tcga.column.enumerations.json", sep="/"))
os.dataset.enumerations     <- fromJSON(paste("../manifests","os.dataset.enumerations.json"    , sep="/"))

date <- Sys.Date()
process <- "import"

# Class Definitions :: Enumerations -------------------------------------------------------
os.enum.na <- c("", "NA", "[NOTAVAILABLE]","[UNKNOWN]","[NOT AVAILABLE]","[NOT EVALUATED]","UKNOWN","[DISCREPANCY]",
                "NOT LISTED IN MEDICAL RECORD","[NOT APPLICABLE]","[PENDING]","PENDING", "[NOT AVAILABLE]","[PENDING]",
                "[NOTAVAILABLE]","NOT SPECIFIED","[NOT AVAILABLE]|[NOT AVAILABLE]",
                "[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]",
                "[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILABLE]",
                "[NOT AVAILABLE]|[NOT AVAILABLE]|[NOT AVAILosABLE]",
                "[NOT AVAILABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT AVAILABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]|[NOT AVAILABLE]|[NOT APPLICABLE]|[NOT APPLICABLE]","N/A")
os.enum.logical.true  <- c("TRUE","YES","1","Y")
os.enum.logical.false <- c("FALSE","NO","0","N")
os.tcga.ignore.columns <- c("bcr_patient_uuid", 
                            "bcr_drug_uuid","bcr_drug_barcode",
                            "bcr_followup_uuid","bcr_followup_barcode",
                            "bcr_radiation_uuid","bcr_radiation_barcode", 
                            "bcr_omf_uuid", "bcr_omf_barcode",
                            "informed_consent_verified", "form_completion_date", 
                            "project_code", "patient_id")

# aggregate list of unmapped data & cde id mapping
unmapped.List <- list()
cde.df <- data.frame()

# Data Processing Functions :: [Map, Clean, Filter]  -------------------------------------------------------

Map( function(key, value, env=parent.frame()){
  setClass(key)
  setAs("character", key, function(from){ 
    # Convert To Upper + Set NAs  
    from<-toupper(from) 
    from.na<-which(from %in% os.enum.na)
    from[from.na]<-NA    
    
    from.clean <- rep(NA, length(from))
    
    # Return Enum or NA
    standardVals <- names(os.tcga.field.enumerations[[key]])
    for(fieldName in standardVals){
      values <-os.tcga.field.enumerations[[key]][[fieldName]]
      from.clean[ which(from %in% values)] <- paste(from.clean[which(from %in% values)], fieldName, sep=";")
    }
    from.clean <- gsub("^NA;", "", from.clean)
 #   from.clean[from.clean==""] <- NA
    
    if(all(unlist(sapply(from.clean, function(val){strsplit(val, ";")})) %in% c(standardVals, NA)))
      return(from.clean)
    
    # Kill If Not In Enum or Na
    stop(paste(key, " not set due to: ", paste(setdiff(from.clean,c(standardVals, NA)), collapse="..."), " not belonging to ", paste(standardVals, collapse=";")))
  })
}, names(os.tcga.field.enumerations), os.tcga.field.enumerations);

# Class Definitions :: TCGA [ID | DATE | CHAR | NUM | BOOL] -------------------------------------------------------

### TCGA ID
setClass("os.class.tcgaId")
setAs("character","os.class.tcgaId", function(from) {
  as.character(str_replace_all(from,"-","." )) 
})

### TCGA Date
setClass("os.class.tcgaDate");
setAs("character","os.class.tcgaDate", function(from){
  
  # Convert Input Character Vector To Uppercase
  from<-toupper(from) 
  
  # Validate Format + Convert Day-Month to 1-1
  if ((str_length(from)==4) && !is.na(as.integer(from) ) ){
    return(as.numeric(as.POSIXct(paste(from, "-1-1", sep=""), format="%Y-%m-%d")))
    #    return(format(as.Date(paste(from, "-1-1", sep=""), "%Y-%m-%d"), "%m/%d/%Y"))
  }
  
  # Return NA If Validation Fails
  return(NA)
})

### TCGA Character
setClass("os.class.tcgaCharacter");
setAs("character","os.class.tcgaCharacter", function(from){
  
  # Convert Input Character Vector To Uppercase
  from<-toupper(from) 
  
  # Get Indexes Of Fram Where Value Is In NA
  from.na<-which(from %in% os.enum.na)
  
  # Set From Indexes Values To NA
  from[from.na]<-NA 
  
  return(from)
})

### TCGA Numeric Radiation
setClass("os.class.tcgaNumeric.radiation");
setAs("character","os.class.tcgaNumeric.radiation", function(from){
  
  # Convert Input Character Vector To Uppercase
  from<-toupper(from) 
  
  # Get Indexes Of Fram Where Value Is In NA
  from.na<-which(from %in% os.enum.na)
  
  # Set From Indexes Values To NA
  from[from.na]<-NA 
  
  from<- gsub("MCI|MILLICURIES|-MILLICURIE|MCI (3730 MBQ)|MILLICURIES 131-IODINE", "", from)
  trim(from)
  
  from <- as.numeric(from)
  
  if(all(is.numeric(from))) return (from)
  
  # Kill If Not In Enum or Na
  stop(paste("os.class.tcgaNumeric.radiation not properly set: ", from[!is.numeric(from)], collapse=";"))
  
})


### TCGA Numeric
setClass("os.class.tcgaNumeric");
setAs("character","os.class.tcgaNumeric", function(from){
  
  # Convert Input Character Vector To Uppercase
  from<-toupper(from) 
  
  # Get Indexes Of Fram Where Value Is In NA
  from.na<-which(from %in% os.enum.na)
  
  # Set From Indexes Values To NA
  from[from.na]<-NA 
  
  from <- as.numeric(from)
  
  if(all(is.numeric(from))) return (from)
  
  # Kill If Not In Enum or Na
  stop(paste("os.class.tcgaNumeric not properly set: ", from[!is.numeric(from)], collapse=";"))
  
})

### TCGA Boolean
setClass("os.class.tcgaBoolean");
setAs("character","os.class.tcgaBoolean", function(from){
  
  from<-toupper(from) 
  
  from.na<-which(from %in% os.enum.na)
  from[from.na]<-NA  
  
  from.true <- which( from %in% os.enum.logical.true )
  from[from.true] <- "TRUE"
  
  from.false <- which(from %in% os.enum.logical.false )
  from[from.false] <- "FALSE"
  
  from <- as.logical(from)
  
  # Return Enum or NA        
  if( all(from %in% c( TRUE, FALSE, NA))) return( from )
  
  # Kill If Not In Enum or Na
  stop(paste("os.class.tcgaBoolean not properly set: ", setdiff(from,c( TRUE, FALSE, NA )), collapse=";"))
})

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


# IO Utility Functions :: [Batch, Load, Save]  -------------------------------------------------------

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

### Load Function Takes An Import File + Column List & Returns A DataFrame
os.data.load.molecular <- function(inputFile){
  
  mtx<- read.delim(inputFile, header=F)
  
  #orient mtx so row: gene, col: patient/sample
  if(all(grepl("^TCGA", mtx[-1,1]))) { mtx <- t(mtx)}
  colType <- "patient"; rowType <- "gene"

  if(all(grepl("TCGA-\\w{2}-\\w{4}-\\w{2}", mtx[1,-1])))
    colType <- "sample"

  mtx.Data<- get.processed.mtx(mtx, dimension= c("row", "col"))
 
  return(list(rowType=rowType, colType=colType, rows=mtx.Data$rownames, cols=mtx.Data$colnames, data=mtx.Data$data))
 
}

### Load Function Takes An Import File + Column List & Returns A DataFrame
os.data.load.clinical <- function(inputFile, checkEnumerations=FALSE, checkClassType = "character"){
  
  # Columns :: Create List From Url
  header <- readLines(inputFile, n=3)
  columns <- unlist(strsplit(header[1],'\t'));
  cde_ids <- unlist(strsplit(header[3],'\t'));
  cde_ids <- gsub("CDE_ID:", "", cde_ids)

  unMappedData <- list();
  tcga_columns <- columns
  
  if(grepl("clinical_patient_skcm.txt",inputFile)){
  	columns[match("submitted_tumor_site", columns)] = "skcm_tissue_site"
  	columns[match("submitted_tumor_site", columns)] = "skcm_tumor_type"
  }
  if(grepl("follow_up_v2.0_skcm.txt",inputFile)){
  	columns[match("new_tumor_event_type", columns)] = "skcm_tumor_event_type"
  }
  if(grepl("clinical_patient_thca.txt",inputFile)){
    columns[columns=="metastatic_dx_confirmed_by_other"] = "thca_metastatic_dx_confirmed_by_other"
  }
  if(grepl("clinical_patient_kirp.txt",inputFile)){
    columns[columns=="tumor_type"] = "disease_subtype"
  }
  
  # if checkEnumerations - all columns will be read in and assigned 'character' class by default
  # otherwise only classes with defined enumerations will be stored in the mapped table
  if(checkEnumerations) { column_type <- rep("character", length(columns))}
  else                  { column_type <- rep("NULL", length(columns)) }
  

  # assign class types for recognized columns
  #   for each enumerated class type, 
  #     rename matching column to mapped name and assign appropriate type
  os.tcga.classes <- names(os.tcga.column.enumerations)
  for(class.type in os.tcga.classes){
    for(colName in names(os.tcga.column.enumerations[[class.type]])){
      values <-os.tcga.column.enumerations[[class.type]][[colName]]
      matching.values <- which(columns %in% values)
      columns[matching.values ] <- colName
      column_type[ matching.values] <- class.type
    }
  }
  
  # Table :: Read Table From URL
  mappedTable<-read.delim(inputFile,
                          header = FALSE, 
                          skip = 3,
                          dec = ".", 
                          sep = "\t",
                          strip.white = TRUE,
                          check.names=FALSE,
                          numerals = "warn.loss",
                          col.names = columns,
                          colClasses = column_type
  );
  
  if(checkEnumerations) {
    
    # Grab columns matching class type and remove those within the ignore list
    headerWithData <- columns[column_type == checkClassType]
    ignoreCols <- which(headerWithData %in% os.tcga.ignore.columns)
    if(length(ignoreCols > 0))       headerWithData <- headerWithData[- ignoreCols ]
    if(length(headerWithData) == 0)  return(list(mapped=mappedTable, unmapped=unMappedData, "cde"=cbind(tcga_columns,columns,cde_ids, column_type)));
    
    # Discard columns where all values are NA
    DataIndicator <- sapply(headerWithData, function(colName){!all(toupper(mappedTable[,colName]) %in% os.enum.na)})
    headerWithData <- headerWithData[DataIndicator]
    if(length(headerWithData) == 0) return(list(mapped=mappedTable, unmapped=unMappedData, "cde"=cbind(tcga_columns,columns,cde_ids, column_type)));
    
    # Print list of unique values for each column
    unMappedData <- lapply(headerWithData, function(colName){ unique(toupper(mappedTable[,colName]))})
    names(unMappedData) <- headerWithData
    print("---Unused columns")
    print(unMappedData)

  }
  #return(list(mapped=mappedTable, unmapped = unMappedData))
  return(list("mapped"=mappedTable, "unmapped" = unMappedData, "cde"=cbind(tcga_columns,columns,cde_ids, column_type)))
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
    Manifest[Manifest$dataset==datasetName & Manifest$dataType ==dataTypeName,"collections"] <- list(rbind(dataObj$collections[[1]],collection))
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
mapProcess <- function(process){
  
	processFound <-	sapply(os.dataset.enumerations$dataType, function(typeMap){ process %in% unlist(typeMap) })
	numMatches <- length(which(processFound))
	if(numMatches==1)
		return (names(os.dataset.enumerations$dataType)[which(processFound)])

	stop(printf("mapProcess found %d matches for process %s", numMatches, process))
	return(NA)
}
#---------------------------------------------------------
### Batch Is Used To Process Multiple TCGA Files Defined 
os.data.batch <- function(manifest, outputDirectory, ...){
           
    # From Input File: dataframe of datasets, datatypes and list of collections
    datasets <- manifest
    
    Manifest <- data.frame()
      
		# Loop for each file to load
		for (i in 1:nrow(datasets))
		{
			sourceObj <- datasets[i,]
			stopifnot(all(c("dataset", "dataType", "collections") %in% names(sourceObj)))
			cat(sourceObj$dataset, sourceObj$dataType,"\n")

			dataset <- sourceObj$dataset
			collections <- sourceObj$collections[[1]]
			
			for(j in 1:nrow(collections)){

				dataObj <- collections[j,]
				stopifnot(all(c("id", "process", "directory", "file") %in% names(dataObj)))
				
				inputDirectory <- dataObj$directory
				if(!grepl("/$", inputDirectory)) inputDirectory <- paste(inputDirectory, "/", sep="")	
				inputFile <- paste(inputDirectory, dataObj$file, sep = "")

				dataType <- mapProcess(dataObj$process)
				resultObj <- data.frame(dataset = sourceObj$dataset, dataType = dataType)
				
				if(dataType %in%  c("cnv","mut01")){
					# Load Data Frame - map and filter by named columns
					result <- os.data.load.molecular( inputFile = inputFile)

					result$data <- apply(result$data, 2, as.integer)

					resultObj$rowType <- result$rowType; resultObj$colType <- result$colType;
					resultObj$rows <- list(result$rows); resultObj$cols <- list(result$cols);
					resultObj$data <- list(result$data)
				}
				if(dataType %in%  c("clinical")){
					# Load Data Frame - map and filter by named columns
					result <- os.data.load.clinical( inputFile = inputFile, ...)

					resultObj$data <- list(result$data)
				}

				index <- get.new.collection.index(Manifest, dataset, dataType)
				resultObj$id <- index
				outputFile <- paste(dataset, dataType, index, process , sep="_")
				parent <- list(c(sourceObj$dataset, sourceObj$dataType, dataObj$id))

	      newCollection <- data.frame(id=index,
	          								process=process,
	          								date = date,
	          								directory= outputDirectory,
	          								file=outputFile)
        newCollection$parent <- parent
				Manifest <- add.new.collection(Manifest, dataset, dataType, newCollection)
				
				# Save Data Frame
				os.data.save(
						df = resultObj,
						directory=outputDirectory,
						file= outputFile,
						format = "JSON") 

		   }  # collection
		}  # dataset
    return(Manifest)
}


# Run Block  -------------------------------------------------------

if("molecular" %in% commands){
	Manifest <- os.data.batch(  manifest = os.molecular.ucsc.batch, 
							outputDirectory = outputDir_molecular		 )
	os.data.save(
	  df = Manifest,
	  directory="../manifests",
	  file= "os.import.molecular.ucsc.manifest",
	  format = "JSON") 
}
if("clinical" %in% commands){
	Manifest <- os.data.batch(
	  manifest = os.clinical.tcga.batch,
	  outputDirectory = outputDir_clinical,
	  checkEnumerations = FALSE,
	  checkClassType = "os.class.tcgaCharacter")

	os.data.save(
	  df = Manifest,
	  directory="../manifests",
	  file= "os.import.clinical.tcga.manifest",
	  format = "JSON") 

}
