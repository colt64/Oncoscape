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

<<<<<<< HEAD
default: setup run

setup: mapCoordinates

run: import cluster tools

data: mapCoordinates import mds edges scale tools


=======
>>>>>>> develop
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

<<<<<<< HEAD
cluster:
	Rscript networks.calculate.mds.edges.R
		# hobo MDS

#----- Create Tool Specific Views

tools: MnP

MnP:
	Rscript os.data.for.mnp.R
		# hobo MDS

=======
mds:
	Rscript networks.calculate.mds.edges.R cnvsnv
		# hobo MDS
edges:
	Rscript networks.calculate.mds.edges.R mds edges
		# edges

scale10k:
	Rscript networks.scale.mds.genomic.coordinates.R mdsScaled
		# MDS values scaled
	Rscript networks.scale.mds.genomic.coordinates.R geneScaled
		# genomic coordinates scaled
>>>>>>> develop
		
#----- Import into MongoDB

setupMongo:
	npm install co-mongodb co underscore assert fs readline sync-request xmldom
	
syncMongo:
	sudo node sync.js
	