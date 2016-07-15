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

default: setup run

data: mapCoordinates import mds edges scale tools

setup: mapCoordinates setupMongo

run: import mds edges scale tools syncMongo

#-----  Create Files

import:
	Rscript tcga.import.R categories molecular clinical
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
	Rscript networks.calculate.mds.edges.R mds
		# hobo MDS
edges:
	Rscript networks.calculate.mds.edges.R edges
		# edges

scale:
	Rscript networks.scale.mds.genomic.coordinates.R mdsScaled geneScaled chromosomeScaled
		# MDS values scaled
		# geneset coordinates scaled
		# genomic coordinates scaled

#----- Create Tool Specific Views

tools: MnP

MnP:
	Rscript os.data.for.mnp.R
		# hobo MDS

		
#----- Import into MongoDB

setupMongo:
	npm install co-mongodb co underscore assert fs readline sync-request xmldom
	
syncMongo:
	sudo node sync.js
	