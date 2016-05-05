library(NetworkMaker)
library(SttrDataPackage)
library(org.Hs.eg.db)

printf = function (...) print (noquote (sprintf (...)))
options(stringsAsFactors=FALSE)

#install.packages("/Volumes/homes/Lisa/oncoscape/OncoGit/Oncoscape/dataPackages/RCyjs", type="source", repos=NULL)
#browserFile = "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/RCyjs/scripts/rcyjs.html"
#http://oncoscape-static.s3-website-us-west-2.amazonaws.com/

#--------------------------------- make plot data -----------------------------#
diseaseAbbr <-c("BRCA", "LUNG", "LUAD","PRAD","LGG","GBM","LGG.GBM", "PAAD", "COADREAD")
diseaseDataP <- c("TCGAbrca", "TCGAlung","TCGAluad","TCGAprad","TCGAlgg","TCGAgbm","TCGAbrain", "TCGApaad", "TCGAcoadread")

oncoVogel274 <- get(load(paste("../extdata/oncoVogel274.RData", sep="/")))

#----------------------------------------------------------------------------------------------------
create.and.display <- function(includeUnpositionedSamples=TRUE, threshold = NA, regex= NA)
{
#   load(system.file(package="NetworkMaker", "extdata", "genesets.RData"))
#   goi <- sort(unique(c(genesets$tcga.GBM.classifiers, genesets$marker.genes.545)))
#   db <- org.Hs.eg.db
#   tbl <- select(db, columns=c("SYMBOL", "MAP", "CHRLOC"), keytype="SYMBOL", keys=goi)
#   goi <- unique(tbl[!is.na(tbl$MAP),"SYMBOL"]);
#   goi <- goi[-which(goi=="MAPT")]
#	goi <- getAlteredGeneNames(netMaker)


   gistic.scores <-c(-2,-1,0,1, 2)
  
   goi = oncoVogel274
   if(!is.na(regex))  samples <- NA
   else{              samples <- get.filtered.sampleIDs(regex) }
   
   calculateSampleSimilarityMatrix(netMaker, copyNumberValues=gistic.scores, genes = goi, samples=samples)
   #filename <- "MDS.SNV.CNV.tsv"
   #usePrecalculatedSampleSimilarityMatrix(netMaker, filename)

   g <- getSamplesGraph(netMaker, includeUnpositionedSamples)
   rcy <- RCyjs(portRange=6047:6100, quiet=TRUE, graph=g, hideEdges=TRUE)
   httpSetStyle(rcy, system.file(package="NetworkMaker", "extdata", "style.js"))
   tbl.pos <- getSimilarityScreenCoordinates(netMaker, xOrigin=0, yOrigin=0, xMax=6000, yMax=6000)
   setPosition(rcy, tbl.pos)    
   fit(rcy, 100)

   g.chrom <- getChromosomeGraph(netMaker, goi)
   httpAddGraph(rcy, g.chrom)
   httpSetStyle(rcy, system.file(package="NetworkMaker", "extdata", "style.js"))
   tbl.pos <- getChromosomeScreenCoordinates(netMaker, xOrigin=3400, yOrigin=0, yMax=3000, chromDelta=200)
   setPosition(rcy, tbl.pos)
   fit(rcy, 100)

   poi <- names(which(noa(g, "positioned")))
   g.mut <- getMutationGraph(netMaker, goi, poi)
   httpAddGraph(rcy, g.mut)
   hideAllEdges(rcy)

   g.cn <- getCopyNumberGraph(netMaker, goi, poi, gistic.scores)
   httpAddGraph(rcy, g.cn)
   hideAllEdges(rcy)

#   g.splice <- getSplicingGraph(netMaker, goi)
#   httpAddGraph(rcy, g.splice)
#   hideAllEdges(rcy)
   showEdges(rcy, "chromosome")
   fit(rcy)

   httpSetStyle(rcy, system.file(package="NetworkMaker", "extdata", "style.js")) 
   # temporary fix, accomodating orphan genes (not mapped to chromosomes):
   
   unpositioned.nodes <- names(which(!noa(g, "positioned")))
#   selectNodes(rcy, unpositioned.nodes)
   layoutSelectionInGrid(rcy, x=-2000, y=3300, w=1400, h=400)
   fit(rcy)

   return(list(rcy=rcy, g=g, g.mut=g.mut, g.cn=g.cn, unpositioned=unpositioned.nodes))

} # create.and.display
#----------------------------------------------------------------------------------------------------
saveGraph <- function(rcy)
{
   g.markers.json <- getJSON(rcy)   # about 1M
   filename <- "../../../extdata/markers.json.RData"
   printf("saving as %s, %d nodes, %d edges", filename, getNodeCount(rcy), getEdgeCount(rcy))
   save(g.markers.json, file=filename)

} # saveGraph
#----------------------------------------------------------------------------------------------------

dataPackage_dir = "../../../../../../dataPackages/Raw_UCSC";
markerFolder = "inst/import";


for(i in 1:length(diseaseAbbr)){
	diseaseName= diseaseAbbr[i]
	dataFolderName = diseaseDataP[i]
	
	print(diseaseName)

	setwd(paste(dataPackage_dir,dataFolderName, markerFolder, sep="/"))

	eval(parse(text=sprintf("library(%s)", dataFolderName)))
	eval(parse(text=sprintf("dz <- %s()" , dataFolderName)))

	netMaker <- NetworkMaker(dz)
	
	regex = ".01$"; threshold = NA;
	if(diseaseName == "AML"){       regex = ".03$|.09$" }
	else if(diseaseName == "LUAD"){ regex = "TCGA.(17)^.\d{4}.01$" }

	if(diseaseName == "BRCA" | diseaseName == "LGG.GBM")  threshold = -1e-04

	x <- create.and.display(includeUnpositionedSamples=FALSE, regex, threshold) 
	rcy <- x$rcy

	hideAllEdges(rcy)  # deselect any selections
	showEdges(rcy, "chromosome")
	saveGraph(rcy)

}
