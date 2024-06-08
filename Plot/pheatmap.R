# install.packages("R.matlab")
# install.packages("pheatmap")

library(Seurat)
library(R.matlab)
library(pheatmap)
library(dplyr)

bench <- read.csv("./datasets/simulated_data1/other_res/ground_truth.csv")
rownames(bench) <- bench[,1]; bench[,1] <- NULL
bench <- bench/rowSums(bench)
bench <- bench %>% select(B.cells,EPCAM..cells.and.cholangiocytes,Hepatocytes,LSEC,NK.NKT.and.T.cells,Stellate.cells.and.myofibroblasts)
res <- bench

sp_data <- read.csv("./datasets/simulated_data1/X.csv")
rownames(sp_data) <- sp_data[,1]; sp_data[,1] <- NULL
celltype_data <- read.csv("./datasets/simulated_data1/Y.csv")
rownames(celltype_data) <- celltype_data[,1]; celltype_data[,1] <- NULL
res <- read.csv("./datasets/simulated_data1/F.csv") # celltype*spot
res <- as.data.frame(t(res))
res <- res/rowSums(res)
res[res<0.05] <- 0
res <- res/rowSums(res)
colnames(res) <- colnames(celltype_data); rownames(res) <- colnames(sp_data)

## cell2location
cell2location <- read.csv("./datasets/simulated_data1/other_res/cell2location_res.csv")
rownames(cell2location) <- cell2location[,1]; cell2location[,1] <- NULL
res <- cell2location
res <- res/rowSums(res)
new_col <- gsub("q05cell_abundance_w_sf_","",colnames(res))
colnames(res) <- new_col

# RCTD barcode XX.1 -> XX-1
RCTD <- read.csv("./datasets/simulated_data1/other_res/RCTD_res.csv")
rownames(RCTD) <- RCTD[,1]; RCTD[,1] <- NULL
res <- RCTD
new_row <- gsub("[.]","-",rownames(res))
rownames(res) <- new_row

# SPOTlight
SPOTlight <- read.csv("./datasets/simulated_data1/other_res/SPOTlight_res.csv")
rownames(SPOTlight) <- SPOTlight[,1]; SPOTlight[,1] <- NULL
res <- SPOTlight
res <- res/rowSums(res)
new_row <- gsub("[.]","-",rownames(res))
rownames(res) <- new_row

# Redeconve
Redeconve <- read.csv("./datasets/simulated_data1/other_res/Redeconve_res.csv")
rownames(Redeconve) <- Redeconve[,1]; Redeconve[,1] <- NULL
res <- Redeconve
res <- res/rowSums(res)
new_row <- gsub("[.]","-",rownames(res))
rownames(res) <- new_row

# CARD
CARD <- read.csv("./datasets/simulated_data1/other_res/CARD_res.csv")
rownames(CARD) <- CARD[,1]; CARD[,1] <- NULL
CARD <- CARD %>% select(B.cells,EPCAM..cells.and.cholangiocytes,Hepatocytes,LSEC,NK.NKT.and.T.cells,Stellate.cells.and.myofibroblasts)
res <- CARD
res <- res/rowSums(res)

# DestVI
DestVI <- read.csv("./datasets/simulated_data1/other_res/DestVI_res.csv")
rownames(DestVI) <- DestVI[,1]; DestVI[,1] <- NULL
res <- DestVI
res <- res/rowSums(res)

# SpatialDWLS
SpatialDWLS <- read.csv("./datasets/simulated_data1/other_res/SpatialDWLS_res.csv")
SpatialDWLS[,1] <- NULL;rownames(SpatialDWLS) <- SpatialDWLS[,1]; SpatialDWLS[,1] <- NULL
SpatialDWLS <- SpatialDWLS %>% select(B.cells,EPCAM..cells.and.cholangiocytes,Hepatocytes,LSEC,NK.NKT.and.T.cells,Stellate.cells.and.myofibroblasts)
res <- SpatialDWLS
res <- res/rowSums(res)

# DSTG
DSTG <- read.csv("./datasets/simulated_data1/other_res/DSTG_res.csv",header = FALSE)
DSTG_label <- read.csv("./datasets/simulated_data1/other_res/DSTG_label.csv")
DSTG_label[,1] <- NULL;colnames(DSTG) <- colnames(DSTG_label)
DSTG <- DSTG %>% select(B.cells,EPCAM..cells.and.cholangiocytes,Hepatocytes,LSEC,NK.NKT.and.T.cells,Stellate.cells.and.myofibroblasts)
res <- DSTG
res <- res/rowSums(res)


## select spot:963,964,965,988,989,990
show_res = as.data.frame(res)
select_show_res <- show_res[c(963,964,965,988,989,990),]
rownames(select_show_res) <- c('spot_963','spot_964','spot_965','spot_988','spot_989','spot_990')
colnames(select_show_res) <- c("celltype1","celltype2","celltype3","celltype4","celltype5","celltype6")


ph <- pheatmap(mat=select_show_res ,scale="none",cluster_cols = FALSE,cluster_rows=FALSE,
               display_numbers = TRUE,show_rownames=TRUE,angle_col = 45,show_colnames = TRUE,
             cellheight =20,cellwidth = 30,color = colorRampPalette(colors = c("white","red"))(100),
             border_color = "black",fontsize_number=10,fontsize = 12)