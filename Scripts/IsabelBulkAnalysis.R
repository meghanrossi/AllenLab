
read10x <- function(dir) {
  names <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
  # for (i in 1:length(names)) {
  #   R.utils::gunzip(paste0(dir, "/", names[i], ".gz"))
  # }
  file.copy(paste0(dir, "/features.tsv"), paste0(dir, "/genes.tsv"))
  mat <- Seurat::Read10X(dir)
  file.remove(paste0(dir, "/genes.tsv"))
  for (i in 1:length(names)) {
    R.utils::gzip(paste0(dir, "/", names[i]))
  }
  mat
}
library(Seurat)
library(dplyr)
# setwd("/gpfs/analyses/maxsh/Allen/Isabel/190308_10x/Bulk")

###### LOADING DATA #########
# Tell the program where the data directory is
data<-read10x("Data/Source")
so <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 500, 
                           project = "Isabel10xBulk")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = so@data), value = TRUE)
percent.mito <- Matrix::colSums(so@raw.data[mito.genes, ])/Matrix::colSums(so@raw.data)
so <- AddMetaData(object = so, metadata = percent.mito, col.name = "percent.mito")

######## BASIC QC ###########
pdf(file = "QC.pdf",width = 12,height = 12,pointsize=10)
VlnPlot(object = so, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, point.size.use=FALSE)
dev.off()
pdf(file = "QC2.pdf",width = 12,height = 12,pointsize=10)
par(mfrow = c(1, 2))
GenePlot(object = so, gene1 = "nUMI", gene2 = "percent.mito",cex.use = 0.5)
GenePlot(object = so, gene1 = "nUMI", gene2 = "nGene",cex.us=0.5)
dev.off()

######## FILTERING ##########
so <- FilterCells(object = so, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(1000, -Inf), high.thresholds = c(4000, 0.08))
pdf(file = "QCafterFilter.pdf",width = 12,height = 12,pointsize=10)
VlnPlot(object = so, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,point.size.use = FALSE)
dev.off()

######## NORMALIZATION #####
so <- NormalizeData(object = so, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

######## VAR GENE IDENTIFICATION ####
pdf(file = "VarGenes.pdf",width = 12,height = 12,pointsize=10)
so<- FindVariableGenes(object = so, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.025, x.high.cutoff = 4, y.cutoff = 0.6)
dev.off()
######## REMOVING UNWANTED VARIATION ###
so <- ScaleData(object = so, vars.to.regress = c("nUMI", "percent.mito"))

######## PCA #######
so <- RunPCA(object = so, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
pdf(file = "PC1-4.pdf",width = 12,height = 12,pointsize=10)
VizPCA(object = so, pcs.use = 1:4)
dev.off()

pdf(file = "PCA.pdf",width = 12,height = 12,pointsize=10)
PCAPlot(object = so, dim.1 = 1, dim.2 = 2,cols.use=rgb(0, 0, 255, max = 255, alpha = 50))
dev.off()
pdf(file = "PCHeatmap.pdf",width = 12,height = 12,pointsize=10)
PCHeatmap(object = so, pc.use = 1:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
dev.off()
pdf(file = "PCElbow.pdf",width = 12,height = 12,pointsize=10)
PCElbowPlot(object = so)    
dev.off()

########## Clustering ###########
so <- FindClusters(object = so, reduction.type = "pca", dims.use = 1:15, 
                     resolution = 1.0, print.output = 0, save.SNN = TRUE)
so <- RunTSNE(object = so, dims.use = 1:15, do.fast = TRUE)
pdf(file = "Tsne13clusters.pdf",width = 12,height = 12,pointsize=10)
TSNEPlot(object = so)
dev.off()

########## Coloring by average marker expression #######

## lets define the markers
markers=list()
markers$microglia=c("Tmem119","C1qa","C1qc","C1qb","Csf1r","Laptm5","Cd68","Olfml3","Tyrobp","Fcgr3","Alox5ap","Cx3cr1","P2ry13","Hexb","Lst1","Ctss","Coro1a","Ly86","Rtn4rl1","Rgs10","Rnase4","Kcnk12","Sla","Egr2","Ptgs1","Trpv2","Arhgdib","Gpr34","P2ry12","Rgs14","Ptk2b","Rassf5","Spp1")
markers$neuron=c("Reln","Slc17a6","Sst","5330417C22Rik","Dlx1","Npy","Nos1","Islr2","Dlx1as","Pnoc","Grem2","Clstn2","Npas4","Grm2","Nell1","Tmem130","Glra2","Slc32a1","Tbr1","Trim66","Lhx6","Bcl11a","Dpysl5","L1cam","Dync1i1","Tac2","Gpr83","Stmn2","Syt1","Crhbp","Tubb3","Rian","Kcnc2","Vgf","Bcl11b","Syt4","Cdh8","Kctd8","Gpr12","Elavl2","Synpr","Celsr3","Lrrc3b","Crh","Cpne5","Neurod1","Crmp1","Myt1l","Ina","Rgs8","Cpne4","Calb2","Gap43","Cacng2","Kcnip4","Eomes","Scg2","Adcyap1","Ablim3","Gpr22","Gabrg2","Rims3","Car10","Gpr26","Atp1a3","Rem2","Slc6a7","Gad1","Usp29","Cyp4x1","Nrn1","Kcnq5","Cbln4","Kcns2","Efna3","AW551984","Gabra5","Baiap3","Cacna1b","Zcchc12","Zic1","Sv2b","Hspa12a","Rhov","Gria1","Kcna4","Cntnap4","Hcn1","Trhde","Cacna2d1","Doc2b","A830018L16Rik")
markers$opc=c("Lhfpl3","3110035E14Rik","Dcn","Lrrtm3","Sulf1","3632451O06Rik","Gria3","Prkg2","Tmem132d","C1ql3","Nxph1","Gpr17","Cntn4","Cobl","Hrasls","Slitrk3")
markers$oligo=c("Mog","Cldn11","Mobp","Mbp","Plp1","Mag","Mal","Ppp1r14a","Fa2h","Gsn","Gpr62","Tnni1","Kndc1","Ndrg1","Ugt8a","Tspan2","Adamts4","1700047M11Rik","Adssl1","Tmem117","Tmeff2","Sh3gl3","Plekhh1","Kcna1","Bcas1","Cpm")
markers$endo=c("Car4","Itm2a","Ly6a","Egfl7","Cd34","Kitl","Sgpp2","Cdc42ep3","Ebf1","Ankrd37","Arhgef15","Akap2","Myo1b","Mgp","Gcnt2","Tiam1","Plk2","Slc39a10","Smoc2","Abcg2","Grb10","Pde5a","Bok","Myl4","Rasgrp3","Vtn","Tmsb10","9430020K01Rik","Anxa3","Rassf3")
markers$astrocyte=c("Aqp4","Bmpr1b","Plcd4","Ppp1r3c","Slc30a10","Itih3","Slc14a1","Atp13a4","Grm3","Paqr6","Fgfr3","Slc4a4","Itga7","Aldh1l1","St6galnac5","Cbs","Phkg1","Slc6a11","Slc25a34","Sorcs2","Cyp4f15","Entpd2","Grhl1","Dio2","Sox9","Slc15a2","Slc7a2","AI464131","Slc7a10","Slc39a12","Mlc1","Ppp1r3g","Ephx2","Slc1a3","2900052N01Rik","Vcam1","Adhfe1","Rorb","Daam2","Rgs20","Abcd2","Lrig1","Slc7a11","Cldn10","Elovl2","Igsf1","Slc25a18","F3","Fabp7","Emp2","Pdk4","Acsbg1","Kcne1l","Igfbp2","Gpc4","Hes5","Slc1a2","Tlcd1","Dzip1l","Gfap","Tst","Rfx4","Grin2c","Slc13a3","Fjx1","Megf10","Slc27a1","Pla2g7","Mfge8","Pax6","Rdh5","Cyp4f14","Hsd11b1","Lcat","Gpam","Agt","Olfml1","Dbx2","Hapln1","Amot","Gjb6","Mtmr11","Ccnd2","Mapk4","Acsl6","Sfxn5","Ptch1","Htra1","Rasa2","Gabrg1","Id4","Scara3","Atp1b2","Acot11","Klf15","Tpbg","Ptprz1","Fgfr1","Vegfa","Fxyd1","Lonrf3","Nhsl1")
so=AddMetaData(so,Matrix::colMeans(so@scale.data[rownames(so@scale.data) %in% markers$microglia, ]),col.name = "Avg. Microglia")
so=AddMetaData(so,Matrix::colMeans(so@scale.data[rownames(so@scale.data) %in% markers$neuron, ]),col.name = "Avg. Neuron")
so=AddMetaData(so,Matrix::colMeans(so@scale.data[rownames(so@scale.data) %in% markers$opc, ]),col.name = "Avg. OPC")
so=AddMetaData(so,Matrix::colMeans(so@scale.data[rownames(so@scale.data) %in% markers$oligo, ]),col.name = "Avg. Oligo")
so=AddMetaData(so,Matrix::colMeans(so@scale.data[rownames(so@scale.data) %in% markers$endo, ]),col.name = "Avg. Endo")
so=AddMetaData(so,Matrix::colMeans(so@scale.data[rownames(so@scale.data) %in% markers$astrocyte, ]),col.name = "Avg. Astrocyte")
pdf(file = "TsneCellTypes.pdf",width = 12,height = 12,pointsize=10)
FeaturePlot(object=so,features.plot=c("Avg. Microglia","Avg. Neuron","Avg. OPC","Avg. Oligo","Avg. Endo","Avg. Astrocyte"),cols.use=c("grey","blue"),reduction.use = "tsne")
dev.off()

pdf(file = "TsneAstrocyteGenes.pdf",width = 12,height = 120,pointsize=10)
FeaturePlot(object=so,features.plot=markers$astrocyte,cols.use=c("grey","blue"),reduction.use = "tsne")
dev.off()

######## Saving the object #####
saveRDS(so, file = "Version1.rds")
######## FINDING MARKERS #######
so.markers <- FindAllMarkers(object = so, only.pos = TRUE, min.pct = 0.1)
top10 <- so.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf(file = "MarkerGenesInClusters.pdf",width = 20,height = 20,pointsize=10)
DoHeatmap(object = so, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()


