####---------------------
##  Compare CNV/SNV Input


## Hamid
LGG.GBM.mut <- get(load("~/Desktop/OncoGit/Oncoscape/dataPackages/archive/UCSC/TCGA_BRAIN/validation/HoBo.LGG.GBM.mut.RData"))
# 39356 x 746
LGG.GBM.cnv <- get(load("~/Desktop/OncoGit/Oncoscape/dataPackages/archive/UCSC/TCGA_BRAIN/validation/HoBo.LGG.GBM.cnv.RData"))
# 24174 x 1034

## Oncoscape: join of lgg & gbm tables
brain.cnv <- get(load("~/Desktop/OncoGit/Oncoscape/dataPackages/archive/UCSC/TCGA_BRAIN/mtx.ucsc.cnv.RData"))
# 24776 x 1090
brain.mut <- get(load("~/Desktop/OncoGit/Oncoscape/dataPackages/archive/UCSC/TCGA_BRAIN/mtx.ucsc.mut01.RData"))
# 39867 x 818

> dim(brain.mut) - dim(LGG.GBM.mut)
[1] 511  72
> dim(brain.cnv) - dim(LGG.GBM.cnv)
[1] 602  56


colnames(brain.mut) <- gsub("\\.\\d{2}$","", colnames(brain.mut))
colnames(brain.cnv) <- gsub("\\.\\d{2}$","", colnames(brain.cnv))

setdiff(rownames(LGG.GBM.mut), rownames(brain.mut))

setdiff(colnames(LGG.GBM.mut), colnames(brain.mut))


rownames(LGG.GBM.cnv) <- gsub("\\|.+", "", toupper(rownames(LGG.GBM.cnv)))
rownames(brain.cnv) <- gsub("\\|.+", "", toupper(rownames(brain.cnv)))

setdiff(rownames(LGG.GBM.cnv), rownames(brain.cnv))



####---------------------
##  Compare MDS output

coll <- mongo.find.all(mongo, "oncoscape.brain_mds_ucsc-hobo_mds-allgenes-cnv-mut01")
mtx <- t(sapply(coll[[1]]$data, function(ptData) c(ptData$x, ptData$y)))

mtx <- mtx *-1

plot(mtx, pch=19)

load("/Volumes/homes/HollandLabShared/Hamid/HoBo/autoHoboPlotData/HoBoPlotData7Apr2015.RData")
h.plot <- plotData$joint.SNA.CNA
rownames(h.plot) <- gsub("\\.", "\\-", rownames(h.plot))


diff = mtx[rownames(h.plot),] - h.plot


coll <- mongo.find.all(mongo, "oncoscape.brain_mds_ucsc-hobo_mds-allgenes-cnv-mut01")
