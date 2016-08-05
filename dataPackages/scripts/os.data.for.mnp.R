###
#       
###

# Library Imports ---------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

source("common.R")

# Configuration -----------------------------------------------------------

date <- as.character(Sys.Date())
<<<<<<< HEAD
scaleFactor = 100000
=======
scaleFactor = 10000
>>>>>>> HoBo

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
<<<<<<< HEAD
	  data_coll <- mongo.find.one(mongo, paste("oncoscape", collection$collection, sep="."))
=======
	  data_coll <- mongo.find.all(mongo, paste("oncoscape", collection$collection, sep="."))
>>>>>>> HoBo
    if(length(data_coll)==0){
      print(paste("ERROR: collection not found - ", collection$collection, sep=""))
      next;
    }
<<<<<<< HEAD
	   mongo.insert(mongo, "oncoscape.render_patient", data_coll)
	}

	datatypeName= "color"
#	mds_colls <- mongo.find.all(mongo, "oncoscape.manifest", query=list(dataType="color"))
	
#	for(collection in mds_colls){
#	  scale <- collection$process[[1]]$scale
#	  if(is.null(scale) || scale != scaleFactor) next;
#	  data_coll <- mongo.find.one(mongo, paste("oncoscape", collection$collection, sep="."))
#	  if(length(data_coll)==0){
#	    print(paste("ERROR: collection not found - ", collection$collection, sep=""))
#	    next;
#	  }
#	  mongo.insert(mongo, "oncoscape.render_patient", data_coll)
#	}
	
=======
	   mongo.insert(mongo, "oncoscape.render_patient", data_coll[[1]])
	}
	
	
	cat_colls <- mongo.find.all(mongo, "oncoscape.manifest", query=list(dataType="colorCategory"))
	
	for(collection in cat_colls){
	  data_coll <- mongo.find.one(mongo, paste("oncoscape", collection$collection, sep="."))
	  if(length(data_coll)==0){
	    print(paste("ERROR: collection not found - ", collection$collection, sep=""))
	    next;
	  }
	  mongo.insert(mongo, "oncoscape.render_patient", data_coll)
	}
>>>>>>> HoBo
	
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

<<<<<<< HEAD
  collection <- mongo.find.all(mongo, "oncoscape.manifest", 
                               query=list(dataset="hg19", dataType="genesets", process=list(scale=scaleFactor)))[[1]]
  
  data_coll <- mongo.find.one(mongo, paste("oncoscape", collection$collection, sep="."))
  mongo.insert(mongo, "oncoscape.render_chromosome", data_coll)
  
=======
  genesets <- mongo.find.all(mongo, "oncoscape.manifest", 
                               query=list(dataset="hg19", dataType="genesets", process=list(scale=scaleFactor)))[[1]]
  geneset_coll <- mongo.find.all(mongo, paste("oncoscape", genesets$collection, sep="."))
  
  for(collection in geneset_coll){
    mongo.insert(mongo, "oncoscape.render_chromosome", collection)
  }
>>>>>>> HoBo
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
<<<<<<< HEAD
commands <- c("pca")

mongo <- connect.to.mongo()

if("patient" %in% commands) 
 os.save.ptLayouts()
if("chromosome" %in% commands) 
  os.copy.chromosome.layout()
if("pca" %in% commands) 
  os.save.pca()
=======
commands <- c("patient", "chromosome")

mongo <- connect.to.mongo()

#if("patient" %in% commands) 
 os.save.ptLayouts(scaleFactor=100000)
#if("chromosome" %in% commands) 
#  os.copy.chromosome.layout(scaleFactor=100000)
#if("pca" %in% commands) 
#  os.save.pca()
>>>>>>> HoBo

close.mongo(mongo)