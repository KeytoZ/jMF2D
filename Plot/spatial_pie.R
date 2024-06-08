library(Seurat) 
library(ggplot2)
library(dplyr)
library(tibble)
library(R.matlab)

set.seed(321)

sp_data <- read.csv("./datasets/breast_cancer/X.csv")
rownames(sp_data) <- sp_data[,1]; sp_data[,1] <- NULL
celltype_data <- read.csv("./datasets/breast_cancer/Y.csv")
rownames(celltype_data) <- celltype_data[,1]; celltype_data[,1] <- NUL
res <- read.csv("./datasets/breast_cancer/F.csv") # celltype*spot
res <- as.data.frame(t(res))
res <- res/rowSums(res)
colnames(res) <- colnames(celltype_data); rownames(res) <- colnames(sp_data)

# spatial location
location <- read.csv("./datasets/breast_cancer/spatial/tissue_positions_list.csv")
rownames(location) <- location[,1]; location[,1] <- NULL
location <- location[rownames(location) %in% rownames(res), ]
location <- location[,c(4,5)]

# image
image_meta_data <- Read10X_Image("./datasets/breast_cancer/spatial/")
image_row_col <- image_meta_data@coordinates[c("imagerow","imagecol")]
image_row_col <- dplyr::mutate(image_row_col,imagerow_scaled = imagerow * image_meta_data@scale.factors$lowres, imagecol_scaled = imagecol * image_meta_data@scale.factors$lowres)

res['barcodes'] <- rownames(res)
image_row_col['barcodes'] <- rownames(image_row_col)
spatial_coord <- dplyr::inner_join(res,image_row_col,by="barcodes")
rownames(spatial_coord) <- spatial_coord[,"barcodes"]; spatial_coord[,"barcodes"] <- NULL

img <- png::readPNG("./datasets/breast_cancer/spatial/tissue_lowres_image.png")
img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))


# plot
defined_cols = c("#466791","#60bf37","#5a51dc","#d49f36","#db37aa","#6b2940","#417d61", "#da8a6d","#a79cd4"
                 ,"#507f2d","#db37aa","#84b67a","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb",
                 "#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367",
                 "#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
temp_spatial_coord = spatial_coord
show_celltypes = colnames(celltype_data)

scatterpie_pie <- suppressMessages(ggplot2::ggplot() + ggplot2::annotation_custom(grob = img_grob,
                  xmin = 0   , xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
                  scatterpie::geom_scatterpie(data = temp_spatial_coord, ggplot2::aes(x = imagecol_scaled
                 , y = imagerow_scaled), cols =show_celltypes, color = NA, alpha = 1, pie_scale = 0.5)
                 + ggplot2::scale_y_reverse() + ggplot2::ylim(nrow(img), 0) + ggplot2::xlim(0, ncol(img)) 
                 + cowplot::theme_half_open(11, rel_small = 1) + ggplot2::theme_void() 
                 + ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") 
                 + ggplot2::theme(legend.key.size = unit(15,'pt'),legend.title = element_text(size = 10),legend.text = element_text(size = 10))
                 + scale_fill_manual(values = defined_cols[1:length(show_celltypes)]))
print(scatterpie_pie)