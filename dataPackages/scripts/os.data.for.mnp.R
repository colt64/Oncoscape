###
#       
###

# Library Imports ---------------------------------------------------------
library(jsonlite)

# Configuration -----------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

os.mds.manifest   <- fromJSON("../manifests/os.mds.2016-06-27.manifest.json")
os.edges.manifest   <- fromJSON("../manifests/os.edges.2016-06-27.manifest.json")

geneset_file <- "../data/molecular/hg19/hg19_genesets_1_hgnc.json"
genesets <- fromJSON(geneset_file)

output.dir <- "../data/tools/markerPatient/"

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
os.save.categories <- function(){

	dataset = "gbm"
	datatype= "colorCategory"
	MnP.categories <- list()
	# 
	## Patient Colors by Diagnosis
	load('~/Desktop/OncoGit/Oncoscape/r_modules/dataPackages/TCGAbrain/inst/extdata/tumorDiagnosis.RData')
	name="diagnosis"
	categories.list <- get.category.data(name=name, table=diagnosis, cat.col.name="diagnosis", color.col.name="color")
	MnP.categories[[1]] <- data.frame(dataset=dataset, datatype=datatype, name=name)
	MnP.categories[[1]]$data=list(categories.list)
	rm(categories.list)
	
	## Patient Colors by Glioma8
	load('~/Desktop/OncoGit/Oncoscape/r_modules/dataPackages/TCGAbrain/inst/extdata/ericsEightGliomaClusters.RData')
	name="glioma8"
	categories.list <- get.category.data(name=name, table=tbl.glioma8, cat.col.name="cluster", color.col.name="color")
	MnP.categories[[2]] <- data.frame(dataset=dataset, datatype=datatype, name=name)
	MnP.categories[[2]]$data=list(categories.list)
	rm(categories.list)
	
	## Patient Colors by Glioma8
	load('~/Desktop/OncoGit/Oncoscape/r_modules/dataPackages/TCGAbrain/inst/extdata/metabolicExpressionStemness.RData')
	name="metabolicExpressionStemness"
	categories.list <- get.category.data(name=name, table=tbl.expression, cat.col.name="cluster", color.col.name="color")
	MnP.categories[[3]] <- data.frame(dataset=dataset, datatype=datatype, name=name)
	MnP.categories[[3]]$data=list(categories.list)
  rm(categories.list)
	
	## Patient Colors by Glioma8
	load('~/Desktop/OncoGit/Oncoscape/r_modules/dataPackages/TCGAbrain/inst/extdata/tumorGrade.RData')
	name="tumorGrade"
	categories.list <- get.category.data(name=name, table=tbl.grade, cat.col.name="cluster", color.col.name="color")
	MnP.categories[[4]] <- data.frame(dataset=dataset, datatype=datatype, name=name)
	MnP.categories[[4]]$data=list(categories.list)
	rm(categories.list)
	
	## Patient Colors by Glioma8
	load('~/Desktop/OncoGit/Oncoscape/r_modules/dataPackages/TCGAbrain/inst/extdata/verhaakGbmClustersAugmented.RData')
	name="verhaakPlus1"
	categories.list <- get.category.data(name=name, table=tbl.glioma8, cat.col.name="cluster", color.col.name="color")
	MnP.categories[[5]] <- data.frame(dataset=dataset, datatype=datatype, name=name)
	MnP.categories[[5]]$data=list(categories.list)
	rm(categories.list)
	
	os.data.save(MnP.categories, output.dir, "os.MnP.categories.data", format="JSON")
	
}

#----------------------------------------------------------------------------------------------------
os.save.ptLayouts <- function(manifest){

	datasetNames = c("gbm", "brca")
	datatypeName= "ptLayout"
	MnP.ptLayouts <- list()
	# 
	## Patient Layout by MDS CNV/SNV OV
	name = "MDS-CNV;SNV-OV"
	
	for(i in 1:length(datasetNames)){
	  datasetName = datasetNames[i]
  	mdsCollections <- subset(manifest, dataset==datasetName & dataType=="mds")$collections[[1]]
	  ptLayoutObj <- subset(mdsCollections, 
  		process$calculation=="mds" & 
	  	process$geneset =="oncoVogel274" &
		  all(c("cnv", "mut01") %in% unlist(process$input)))
		
	  ptLayout.data <- fromJSON(paste(ptLayoutObj$directory, ptLayoutObj$file,".json", sep=""))
	  MnP.ptLayouts[[i]] <- data.frame(dataset=datasetName, datatype=datatypeName, name=name)
	  MnP.ptLayouts[[i]]$data=list(ptLayout.data)
  }
	os.data.save(MnP.ptLayouts, output.dir, "os.MnP.ptLayouts.data", format="JSON")
		
}

#----------------------------------------------------------------------------------------------------
os.save.network.edges <- function(manifest){
  
  datasetNames = c("brca")
  datatypeName= "networkEdges"
  MnP.edges <- list()
  # 
  ## Patient Layout by MDS CNV/SNV OV
  
  for(i in 1:length(datasetNames)){
    datasetName = datasetNames[i]
    edgeCollections <- subset(manifest, dataset==datasetName & dataType=="network")$collections[[1]]
    
    # get Edge sets for each Gene set
    for(j in 1:length(genesets)){			
      genesetName <- names(genesets)[j]
      edgeObj <- subset(edgeCollections, 
                            process$geneset ==genesetName &
                            all(c("cnv", "mut01") %in% unlist(process$edgeType)))
      edgeObj <- edgeObj[unlist(lapply(edgeObj$process$edgeType, function(eType) all(eType== c("cnv", "mut01")))),]
    
      edge.data <- list(
          edges = fromJSON(paste(edgeObj$directory, edgeObj$file[[1]][1],".json", sep="")),
          ptDegree = fromJSON(paste(edgeObj$directory, edgeObj$file[[1]][2],".json", sep="")),
          geneDegree = fromJSON(paste(edgeObj$directory, edgeObj$file[[1]][3],".json", sep=""))
      )          
      MnP.edges[[j]] <- data.frame(dataset=datasetName, datatype=datatypeName, name=genesetName)
      MnP.edges[[j]]$data=list(edge.data)
    }
  }
  os.data.save(MnP.edges, output.dir, "os.MnP.edges.data", format="JSON")
}

##----------------------------

# os.save.categories()
# os.save.ptLayouts(os.mds.manifest)
  os.save.network.edges(os.edges.manifest)
