library(RUnit)

load(file="../../../../TCGAgbm/inst/extdata/mtx.mut.RData")
tbl.TCGAgbm <- mtx.mut
checkEquals(dim(tbl.TCGAgbm), c(291, 6698))

load(file="../../../../TCGAlgg/inst/extdata/mtx.mut.RData")
tbl.TCGAlgg <- mtx.mut
checkEquals(dim(tbl.TCGAlgg), c(289, 6154))

genes.gbm <- colnames(tbl.TCGAgbm)
genes.lgg <- colnames(tbl.TCGAlgg)

Addgene <- setdiff(genes.gbm, genes.lgg)
if(length(Addgene)>0) tbl.TCGAlgg <- cbind(tbl.TCGAlgg, matrix(NA, nrow=nrow(tbl.TCGAlgg), ncol=length(Addgene), dimnames = list(rownames(tbl.TCGAlgg),Addgene)))
Addgene <- setdiff(genes.lgg, genes.gbm)
if(length(Addgene)>0) tbl.TCGAgbm <- cbind(tbl.TCGAgbm, matrix(NA, nrow=nrow(tbl.TCGAgbm), ncol=length(Addgene), dimnames = list(rownames(tbl.TCGAgbm),Addgene)))

HugoGenes <- unique(c(genes.gbm, genes.lgg))
checkEquals(length(HugoGenes), ncol(tbl.TCGAgbm))
checkEquals(length(HugoGenes), ncol(tbl.TCGAlgg))

tbl.TCGAlgg <- tbl.TCGAlgg[, colnames(tbl.TCGAgbm)]

mtx.mut <- rbind(tbl.TCGAgbm, tbl.TCGAlgg)

gene = "MUC16"
sample = "TCGA.DU.8168.01"
checkEquals(mtx.mut[sample, gene], "P9537T")
checkEquals(nrow(mtx.mut), nrow(tbl.TCGAgbm) + nrow(tbl.TCGAlgg) )

save(mtx.mut, file="../../extdata/mtx.mut.RData")




# indx <- match(rownames(tcgaLGG.cnv), rownames(tcgaGBM.cnv))
# iLGG <- which(!is.na(indx))
# iGBM <- indx[!is.na(indx)]
# LGG.GBM.cnv <- cbind(tcgaLGG.cnv[iLGG, ], tcgaGBM.cnv[iGBM, ])			# 24174 x 1034

# indx <- match(rownames(tcgaLGG.mut), rownames(tcgaGBM.mut))
# iLGG <- which(!is.na(indx))
# iGBM <- indx[!is.na(indx)]
# LGG.GBM.mut <- cbind(tcgaLGG.mut[iLGG, ], tcgaGBM.mut[iGBM, ])			# 39356 x 746

