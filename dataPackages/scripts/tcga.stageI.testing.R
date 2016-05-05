# TCGA Pipeline Stage I Testing 
# Testing Aspects: tables existance & usage; 
#                  Field collection;
#                  Field Data Type consistency.
# Date: May 3rd, 2016
source("tcga.import.R")

library(RUnit)
# for a specified study
os.test.study = ""
os.test.stageI.dir = "../tcga.clean"
os.test.stageI.tables = dir(path=os.test.stageI.dir, pattern=paste(os.test.study,".*RData", sep=""))


inputFiles <- read.delim(inputFile, sep="\t", header=TRUE)
# check tables
os.test.stageI.tables.getnames <- function(inputFiles, rowIndex){
  os.test.study <- inputFiles[rowIndex, os.data.batch.inputFile.studyCol]
  os.test.stageI.dir = "../tcga.clean"
  os.test.stageI.tables = dir(path=os.test.stageI.dir, pattern=paste(os.test.study,".*RData", sep=""))
  os.test.stageI.tables.names <- gsub(paste(os.test.study,"_|.RData", sep=""), "", os.test.stageI.tables)
  return(os.test.stageI.tables.names)
}


os.test.stageI.tables.check <- function(inputFiles, rowIndex){
  target.files <- inputFiles[rowIndex,][which(!is.na(inputFiles[rowIndex,]))]
  target <- names(target.files)[-c(1,2)]
  os.test.study <- target.files[os.data.batch.inputFile.studyCol]
  os.test.stageI.tables.names <- os.test.stageI.tables.getnames(inputFiles, rowIndex)
  checkEquals(sort(target), sort(os.test.stageI.tables.names))
}

# check fields
  # Get used raw Col Names
  rawColNames = c()
  for(i in 1:length(os.tcga.column.enumerations)){
    rawColNames = c(rawColNames, names(os.tcga.column.enumerations[[i]]))
  }
  # Get converted Col Names
  #pattern = paste(paste(names(os.tcga.column.enumerations), collapse=".|"), ".", sep="")
  stageIcol = unlist(os.tcga.column.enumerations, recursive=FALSE)
  for(i in 1:length(os.tcga.column.enumerations)) {
    names(stageIcol) <- gsub(paste(names(os.tcga.column.enumerations)[i], ".", sep=""), "", names(stageIcol))
  }
  stageIColNames <- c(names(stageIcol), "bcr_patient_uuid","bcr_followup_barcode","bcr_followup_uuid","form_completion_date")

   

  # load .RData
  url = paste(target.files$directory, target.files[[target[i]]], sep="/") 
  for(i in 1:length(os.test.stageI.tables)){
    # processed col names
    load(paste(os.test.stageI.dir, os.test.stageI.tables[i], sep="/"))
    checkTrue(all(names(df) %in% stageIColNames))

    # raw col capture checking: there are no missing columns 
    rawTable <- read.delim(url)
    rawTable <- rawTable[-c(2,3),]
    length(which(names(rawTable) %in% rawColNames)) 
  }

  # 

os.test.stageI.fields.check <- function(table, rowIndex){
  
}

# check Data Type



# batch testing
os.test.stageI.batch <- function(inputFile){
  inputFiles <- read.delim(inputFile, sep="\t", header=TRUE)
  studies <- inputFiles[os.data.batch.inputFile.studyCol]
  for (rowIndex in 1:nrow(inputFiles)){  
    # checking table
    res <- try({
      os.test.stageI.tables.check(inputFiles,rowIndex)
      })
    # checking 
  }
}

# execution 
os.test.stageI.batch(inputFile = os.data.batch.inputFile)


