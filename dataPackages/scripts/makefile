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

import:
	Rscript tcga.import.R
		# - reads ucsc.filename.manifest
		# - writes mol collections metatable
		# - writes UCSC molecular data json
		# - TO DO: cbio

mapCoordinates:
	Rscript hg19.save.genomic.coordinates.R
		# genomic coordinates

	#### TO DO: grab more genesets

#----- Calculate Values

mds:
	Rscript networks.calculate.mds.edges.R mds cnvsnv
		# hobo MDS
edges:
	Rscript networks.calculate.mds.edges.R edges
		# edges

scale10k:
	Rscript networks.scale.mds.genomic.coordinates.R mdsScaled geneScaled chromosomeScaled
		# MDS values scaled
		# geneset coordinates scaled
		# genomic coordinates scaled
		
#----- Import into MongoDB

setupMongo:
	npm install co-mongodb co underscore assert fs readline sync-request xmldom
	
syncMongo:
	sudo node sync.js
	