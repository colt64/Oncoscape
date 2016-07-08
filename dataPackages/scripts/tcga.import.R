###
#
#       This Script Executes Basic Processing On TCGA Files
#       Specifically It Types, Uppercases and In Cases Enforces Enumeration Types
#       
###


# Configuration -----------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

source("common.R")
source("os.tcga.mappings.R")

#commands <- c("molecular", "clinical", "categories")
commands <- c("molecular")
#commands <- c("clinical")

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

os.molecular.batch   <- fromJSON("../manifests/os.molecular.manifest.json")
os.clinical.tcga.batch    <- fromJSON("../manifests/os.tcga.clinical.manifest.json")

output.molecular.dir   <- "../data/molecular/clean/"
output.clinical.dir    <- "../data/clinical/clean/"
output.categories.dir  <- "../data/categories/"
output.manifest.dir    <- "../manifests/"

output.molecular.manifest <- "os.import.molecular.ucsc.manifest"
output.clinical.manifest <- "os.import.clinical.tcga.manifest"


date <- as.character(Sys.Date())
process <- "import"

# -------------------------------------------------------
# aggregate list of unmapped data & cde id mapping
unmapped.List <- list()
cde.df <- data.frame()

# Data Processing Functions :: [Map, Clean, Filter]  -------------------------------------------------------
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

### Load Function Takes An Import File + Column List & Returns A DataFrame
os.data.load.molecular <- function(inputFile){
  
  mtx <- matrix();
  
  if(grepl("\\.RData$",inputFile)){
    mtx <- get(load(inputFile))
    if(all(grepl("^TCGA", rownames(mtx)))) { mtx <- t(mtx)}
    colType <- "patient"; rowType <- "gene"

    if(all(grepl("TCGA-\\w{2}-\\w{4}-\\w{2}", colnames(mtx))))
      colType <- "sample"
    
    rownames <- rownames(mtx); colnames <- colnames(mtx)
    colnames <- gsub("\\.", "-", colnames); 
    dimnames(mtx) <- NULL
    mtx.Data<- list(rownames=rownames, colnames=colnames, data=mtx)
    
        
  } else{ 
    mtx<- read.delim(inputFile, header=F) 
    #orient mtx so row: gene, col: patient/sample
    if(all(grepl("^TCGA", mtx[-1,1]))) { mtx <- t(mtx)}
    colType <- "patient"; rowType <- "gene"

    if(all(grepl("TCGA-\\w{2}-\\w{4}-\\w{2}", mtx[1,-1])))
      colType <- "sample"
    
    mtx.Data<- get.processed.mtx(mtx, dimension= c("row", "col"))
    
  }
 
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
  
  if(grepl("../archive/clinical/nationwidechildrens.org_clinical_patient_skcm.txt",inputFile)){
  	columns[match("submitted_tumor_site", columns)] = "skcm_tissue_site"
  	columns[match("submitted_tumor_site", columns)] = "skcm_tumor_type"
  }
  if(grepl("../archive/clinical/nationwidechildrens.org_follow_up_v2.0_skcm.txt",inputFile)){
  	columns[match("new_tumor_event_type", columns)] = "skcm_tumor_event_type"
  }
  if(grepl("../archive/clinical/nationwidechildrens.org_clinical_patient_thca.txt",inputFile)){
    columns[columns=="metastatic_dx_confirmed_by_other"] = "thca_metastatic_dx_confirmed_by_other"
  }
  if(grepl("../archive/clinical/nationwidechildrens.org_clinical_patient_kirp.txt",inputFile)){
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
  return(list("mapped"=mappedTable, "unmapped" = unMappedData, "cde"=cbind(tcga_columns,columns,cde_ids, column_type)))
}
#---------------------------------------------------------
### Batch Is Used To Process Multiple TCGA Files Defined 
os.data.batch <- function(manifest, outputDirectory, ...){
           
    # From Input File: dataframe of datasets, datatypes and list of collections
    datasets <- manifest
    Manifest <- data.frame()
      
		# Loop for each file to load
		for (i in 1:nrow(datasets)){
			sourceObj <- datasets[i,]
			stopifnot(all(c("dataset", "dataType", "collections") %in% names(sourceObj)))
			cat(sourceObj$dataset, sourceObj$dataType,"\n")

			collections <- sourceObj$collections[[1]]
			
			for(j in 1:nrow(collections)){

				dataObj <- collections[j,]
				stopifnot(all(c("id", "process", "directory", "file") %in% names(dataObj)))
				
				inputDirectory <- dataObj$directory
				if(!grepl("/$", inputDirectory)) inputDirectory <- paste(inputDirectory, "/", sep="")	
				inputFile <- paste(inputDirectory, dataObj$file, sep = "")

				dataType <- mapProcess(dataObj$process)
				resultObj <- data.frame(dataset = sourceObj$dataset, dataType = dataType)
				
				if(dataType %in%  c("cnv","mut01", "mut", "rna", "protein", "methyl")){
					# Load Data Frame - map and filter by named columns
					result <- os.data.load.molecular( inputFile = inputFile)
        
					if(dataType != "mut"){
					  result$data <- apply(result$data, 2, as.numeric)
					}
					resultObj$rowType <- result$rowType; resultObj$colType <- result$colType;
					resultObj$rows <- list(result$rows); resultObj$cols <- list(result$cols);
					resultObj$data <- list(result$data)
				}
				if(dataType %in%  c("patient", "drug", "radiation", "otherMalignancy", "followUp", "newTumor")){
					# Load Data Frame - map and filter by named columns
					result <- os.data.load.clinical( inputFile = inputFile, ...)
					resultObj$data <- list(result$mapped)
					resultObj$cde <- list(result$cde)
					
					}

				parent <- list(c(sourceObj$dataset, sourceObj$dataType, dataObj$id))
				
				Manifest <- save.collection(Manifest=Manifest, dataset=sourceObj$dataset, dataType=dataType, result=resultObj,
				                parent=parent, process=process,processName=process, outputDirectory=outputDirectory)
				
		   }  # collection
		}  # dataset
    return(Manifest)
}

#----------------------------------------------------------------------------------------------------
get.category.data<- function(name, table, cat.col.name, color.col.name){
  
  catNames <- unique(table[,cat.col.name])
  categories.type.list <- lapply(catNames, function(cat.name){
    matches <- table[table[,cat.col.name]==cat.name,]
    data <- data.frame(	name=cat.name, 
                        color=unique(matches$color))
    data$values = list(rownames(matches))
    return(data)
  })
  return (categories.type.list)
}
#----------------------------------------------------------------------------------------------------
add.category.fromFile <- function(file, name, col.name, dataset, datatype){
  
  tbl <- get(load(file))
  categories.list <- get.category.data(name=name, table=tbl, cat.col.name=col.name, color.col.name="color")
  df <- data.frame(dataset=dataset, datatype=datatype, name=name)
  df$data=list(categories.list)
  return(df)
}

#----------------------------------------------------------------------------------------------------
os.save.categories <- function(output.dir, datasets = c("gbm")){
  
  color.categories <- list()
  datatype= "colorCategory"
  
  if("gbm" %in% datasets){  

  ## Patient Colors by Diagnosis, glioma8, tumorGrade, verhaak
  color.categories <- list(
    add.category.fromFile(file='../archive/categories/brain/tumorDiagnosis.RData', name="diagnosis", col.name="diagnosis", dataset="gbm", datatype=datatype) ,
    add.category.fromFile(file='../archive/categories/brain/ericsEightGliomaClusters.RData', name="glioma8", col.name="cluster", dataset="gbm", datatype=datatype) ,
    add.category.fromFile(file='../archive/categories/brain/metabolicExpressionStemness.RData', name="metabolicExpressionStemness", col.name="cluster", dataset="gbm", datatype=datatype) ,
    add.category.fromFile(file='../archive/categories/brain/tumorGrade.RData', name="tumorGrade", col.name="cluster", dataset="gbm", datatype=datatype) ,
    add.category.fromFile(file='../archive/categories/brain/verhaakGbmClustersAugmented.RData', name="verhaakPlus1", col.name="cluster", dataset="gbm", datatype=datatype) 
    )
  }
  if("brca" %in% datasets){
    categories.list <- fromJSON("../archive/categories/brca/colorCategories.json")
    color.categories <- c(color.categories, list(categories.list))
  }
  os.data.save(color.categories, output.dir, "os.categories.color.data", format="JSON")
  
}

# Run Block  -------------------------------------------------------
if("categories" %in% commands) 
  os.save.categories(output.dir=output.categories.dir, datasets=c("gbm", "brca"))

if("molecular" %in% commands){
	Manifest <- os.data.batch(  manifest = os.molecular.batch, 
							outputDirectory = output.molecular.dir		 )
		os.data.save(
	  df = Manifest,
	  directory=output.manifest.dir,
	  file= output.molecular.manifest,
	  format = "JSON") 
}
if("clinical" %in% commands){
	Manifest <- os.data.batch(
	  manifest = os.clinical.tcga.batch,
	  outputDirectory = output.clinical.dir,
	  checkEnumerations = FALSE,
	  checkClassType = "os.class.tcgaCharacter")

	os.data.save(
	  df = Manifest,
	  directory=output.manifest.dir,
	  file= output.clinical.manifest,
	  format = "JSON") 

}
