# TCGA Pipeline Stage I Testing 
# Testing Aspects: tables existance & usage; 
#                  Field collection;
#                  Field Data Type consistency.
# Date: May 3rd, 2016

library(RUnit)


source("tcga.import.R")

# for a specified study
os.test.study = "TCGAgbm"
os.test.stageI.dir = "../tcga.clean"
os.test.stageI.tables = dir(path=os.test.stageI.dir, pattern=paste(os.test.study,".*RData", sep=""))
#lapply(paste(os.test.stageI.dir, os.test.stageI.tables, sep="/"), load)
os.test.stageI.tables.names <- gsub("TCGAgbm_|.RData", "", os.test.stageI.tables)


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


# check Data Type



# batch testing
os.test.stageI.batch <- function(inputFile){
  inputFiles <- read.delim(inputFile, sep="\t", header=TRUE)
  studies <- inputFiles[os.data.batch.inputFile.studyCol]
  for (rowIndex in 1:nrow(inputFiles)){  
    res <- try({
      os.test.stageI.tables.check(inputFiles,rowIndex)
      })
  }
}

# execution 
os.test.stageI.batch(inputFile = os.data.batch.inputFile)


