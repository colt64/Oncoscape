#! /bin/bash

#### Data Pipeline
#  	Create Files
#		Clinical Files
#		Molecular Tables (UCSC, cBio)
# 		Genomic Coordinates
#		Genesets
#	Calculate Values
#		Hobo MDS
#		Edges
#		Genomic Coordinates Scaled
#	Import into MongoDB
####

#----------------------------------------------------

#-----  Create Files

writeTCGAfiles:
	Rscript tcga.import.R
		# - reads ucsc.filename.manifest
		# - writes mol collections metatable
		# - writes UCSC molecular data json
		# - TO DO: cbio

saveGenomicCoordinates:
	Rscript hg19.save.genomic.coordinates.R
		# genomic coordinates

	#### TO DO: grab more genesets

#----- Calculate Values

calculateMDS:
	Rscript networks.calculate.mds.edges.R
		# hobo MDS
		# edges

scaleNetworkPlots:
	Rscript < networks.scale.mds.genomic.coordinates.R
		# genomic coordinates scaled
		# MDS values scaled
		
#----- Import into MongoDB

setupMongo:
	npm install co-mongodb co underscore assert fs readline sync-request xmldom
	
syncMongo:
	sudo node sync.js
	