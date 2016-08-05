###
#       
###

# Library Imports ---------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

source("common.R")

# Configuration -----------------------------------------------------------

date <- as.character(Sys.Date())
scaleFactor = 100000

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 0)
	commands <- args

#----------------------------------------------------------------------------------------------------
os.save.ptLayouts <- function(scaleFactor=100000){

	datatypeName= "cluster"
	mds_colls <- mongo.find.all(mongo, "oncoscape.manifest", query=list(dataType="mds"))
	
	for(collection in mds_colls){
	  scale <- collection$process[[1]]$scale
	  if(is.null(scale) || scale != scaleFactor) next;
	  data_coll <- mongo.find.one(mongo, paste("oncoscape", collection$collection, sep="."))
    if(length(data_coll)==0){
      print(paste("ERROR: collection not found - ", collection$collection, sep=""))
      next;
    }
	   mongo.insert(mongo, "oncoscape.render_patient", data_coll)
	}
}

##----------------------------
os.copy.chromosome.layout <- function(scaleFactor=100000){
  
  datatypeName= "cluster"
  collection <- mongo.find.all(mongo, "oncoscape.manifest", 
                              query=list(dataset="hg19", dataType="chromosome"))
  scaled <- which(sapply(collection, function(coll){ 
    if(coll$process == "length") return(FALSE)
    coll$process[[1]][["scale"]]==scaleFactor}))
  
  data_coll <- mongo.find.one(mongo, paste("oncoscape", collection[[scaled]]$collection, sep="."))
  mongo.insert(mongo, "oncoscape.render_chromosome", data_coll)

  collection <- mongo.find.all(mongo, "oncoscape.manifest", 
                               query=list(dataset="hg19", dataType="genesets", process=list(scale=scaleFactor)))[[1]]
  
  data_coll <- mongo.find.one(mongo, paste("oncoscape", collection$collection, sep="."))
  mongo.insert(mongo, "oncoscape.render_chromosome", data_coll)
  
}
#----------------------------------------------------------------------------------------------------
os.save.pca <- function(scaleFactor=NA){
  
  datatypeName= "cluster"
  pca_colls <- mongo.find.all(mongo, "oncoscape.manifest", 
                              query=list(dataType="pcaScores"))
  
  for(collection in pca_colls){
    scale <- collection$process[[1]]$scale
    if(is.na(scaleFactor)){
      if(is.null(scale)){
        data_coll <- mongo.find.one(mongo, paste("oncoscape", collection$collection, sep="."))
        mongo.insert(mongo, "oncoscape.render_pca", data_coll)
      }
      next;
    }
    if(is.null(scale) ||is.na(scale) || scale != scaleFactor) next;
    data_coll <- mongo.find.one(mongo, paste("oncoscape", collection$collection, sep="."))
    mongo.insert(mongo, "oncoscape.render_pca", data_coll)
  }
}

##----------------------------
#commands <- c("patient", "pca", "chromosome")
commands <- c("pca")

mongo <- connect.to.mongo()

if("patient" %in% commands) 
 os.save.ptLayouts()
if("chromosome" %in% commands) 
  os.copy.chromosome.layout()
if("pca" %in% commands) 
  os.save.pca()

close.mongo(mongo)