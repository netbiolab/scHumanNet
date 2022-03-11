#Perform DEG analysis on celltype and on celltype vs condition
library(Seurat)
library(future)
library(magrittr)
library(dplyr)
library(ggpubr)

#load data and convert to sce#############################################################################
setwd('/home3/junhacha/HumanNetv3/Velemeshev2019/merged_celltype/')
counts <- Seurat::Read10X('../10X_out/')
meta <- read.table('../meta.txt', header = T, sep = '\t')
############################################################################################################

#check if barcodes match
rownames(meta) <- meta$cell
meta$cell <- NULL
identical(colnames(counts), rownames(meta))

#merge annotated celltypes to larger granularity
#neu_mat, NRGN neurons are seperated and will be excluded because it is either similar to Excitatoy neurons based on UMAP analysis and is thus considered ambiguous
meta$celltypes_merged <- ifelse(meta$cluster %in% c('AST-FB','AST-PP'), 'Astrocyte',
                                ifelse(meta$cluster %in% c('IN-PV', 'IN-SST','IN-SV2C', 'IN-VIP'), 'Inhibitory',
                                       ifelse(meta$cluster %in% c('L2/3', 'L4', 'L5/6','L5/6-CC'), 'Excitatory',
                                              ifelse(meta$cluster %in% c('Neu-mat','Neu-NRGN-I', 'Neu-NRGN-II'), 'Others', 
                                                     as.character(meta$cluster)))))





meta$celltype_condition <- paste(meta$diagnosis, meta$celltypes_merged, sep = '_')

#make seurat
seurat <- CreateSeuratObject(counts=counts, meta.data = meta)


#normalize using standard method
#==================================================
v.features <- as.numeric(3000)
#===================================================

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize") 
seurat <- ScaleData(seurat, features = rownames(seurat)) 
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = v.features)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))


dims <- as.numeric(20) #use jackstraw for guidance
#======================================================================================
seurat <- RunUMAP(seurat, dims = 1:dims)
p1 <- DimPlot(seurat, reduction = 'umap', group.by = 'celltypes_merged', label = T) #by celltype
p2 <- DimPlot(seurat, reduction = 'umap', group.by = 'cluster', label = T) #by condition

#Neu_mat seems to be smeared across different excitatory neurons..this usuall means they have weired celltype markers in DEG..lets check it out
#all cells are contorled for 5 percent of mito and ribo genes


library(patchwork)
p2 + p1
saveRDS(seurat, '/home3/junhacha/HumanNetv3/Velemeshev2019/seurat_all.rds')

#check expression MECP2, for each condition Violin plot
seurat@active.ident <- as.factor(seurat$celltypes_merged)
VlnPlot(seurat, features='MECP2',split.by = 'diagnosis', pt.size = 0)

#draw violin plot # do this in 
pdf('/home3/junhacha/HumanNetv3/Graphs/SplitViolinn_Brain_CACNA1A.pdf', height=4, width=5)
VlnPlot(seurat, features = 'CACNA1A', split.by ='diagnosis', split.plot = T,pt.size = 0) +
  scale_y_continuous(expand = c(0, 0))
dev.off()


GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ],
                     newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1,
                                                                      'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data,
                                                           draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data),
                                                          c("x", "y")),
                         drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...),
                                      quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin",
                       GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}
give.n <- function(x){return(y = -2.6, label = length(x))}

seurat$CACNA1A <- as.numeric(seurat@assays$RNA['CACNA1A',])
data <- seurat@meta.data
vlnplot <- ggplot(data,aes(celltypes_merged, CACNA1A, fill = diagnosis)) +
  geom_split_violin(color="black", trim=FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.25, notch=F, notchwidth = .4, outlier.shape = NA, coef=0)+
  # Add labels:
  
  labs(x = "merged_celltypes", y = "Normalized expression",
       fill = "HPV") +
  ggtitle("CACNA1 expression") + 
  
  # Define variable colours and theme:
  
  #scale_fill_manual(values = c('#1EB742','#c13744')) +
  #scale_color_manual(values = c('#1EB742','#c13744')) +
  theme_light() +
  RotatedAxis() +
  #stats
  stat_summary(fun.data = give.n, geom = "text") + 
  stat_compare_means(method = "wilcox.test", paired = FALSE) + 
  stat_summary(fun.data = give.n, geom = "text")


vlnplot






#perform DEG for each celltype
#set mullticore
plan()
plan('multiprocess', workers = 16)



markers <- FindAllMarkers(seurat, only.pos = T)
markers <- markers[markers$p_val_adj < 0.05,]

#since humanNetv3 consists of only conding genes, we filter out DEGs to coding genes to make fair comparision
hnv3 <- read.table('~/HumanNetv3/HumanNetv3_networks/HumanNetv3-XC_symbol_linkOnly.tsv', sep = '\t')
coverage.genes <- union(as.character(hnv3[,1]), as.character(hnv3[,2]))

write.table(markers, file='/home3/junhacha/HumanNetv3/Velemeshev2019/DEG/markers_MergedCelltype_logcounts.tsv', 
            sep = '\t' ,quote = F, col.names = T, row.names = T)

markers.f <- markers[markers$gene %in% coverage.genes, ]

write.table(markers.f, file='/home3/junhacha/HumanNetv3/Velemeshev2019/DEG/markers_MergedCelltype_logcounts_CodingGenes.tsv', 
            sep = '\t' ,quote = F, col.names = T, row.names = T)

#check markers for maturing neuron
seurat@active.ident <- as.factor(seurat$cluster)
markers.mat <- FindMarkers(object = seurat, ident.1 = 'Neu-mat', verbose = T)

saveRDS(markers.mat, '/home3/junhacha/HumanNetv3/Velemeshev2019/DEG/markers_celltype_MatNeuron.rds')


#perform DEG Individually for each celltype between condition
setwd('/home3/junhacha/HumanNetv3/Velemeshev2019/DEG/')
seurat@active.ident <- as.factor(seurat$celltypes_merged)
for (cluster in levels(seurat@active.ident)){
  data <- seurat[,seurat@active.ident == cluster]
  markers <- FindMarkers(data, group.by = 'diagnosis', ident.1 = 'ASD', ident.2 = 'Control', test.use = 'wilcox') 
  markers <- markers[markers$p_val_adj < 0.05,]
  markers <- markers[order(markers$avg_logFC, decreasing = T),]
  write.table(markers, quote=F, row.names = T, col.names = T, file = paste0('./markers_ASDvsCTRL_', cluster), sep = '\t')
}

#filter each conditional markers for coding genes
setwd('/home3/junhacha/HumanNetv3/Velemeshev2019/DEG/')
data1 <- read.table('markers_ASDvsCTRL_Astrocyte', header=T, sep='\t')
data1 <- data1[data1$gene %in% coverage.genes ]
write.table(data1, quote=F, row.names = T, col.names = T, file = 'markers_ASDvsCTRL_Astrocyte_filtered', sep = '\t')

data2 <- read.table('markers_ASDvsCTRL_Endothelial', header=T, sep='\t')
data2 <- data2[data2$gene %in% coverage.genes ]
write.table(data2, quote=F, row.names = T, col.names = T, file = 'markers_ASDvsCTRL_Endothelial_filtered', sep = '\t')

data3 <- read.table('markers_ASDvsCTRL_Excitatory', header=T, sep='\t')
data3 <- data1[data3$gene %in% coverage.genes ]
write.table(data3, quote=F, row.names = T, col.names = T, file = 'markers_ASDvsCTRL_Excitatory_filtered', sep = '\t')

data4 <- read.table('markers_ASDvsCTRL_Inhibitory', header=T, sep='\t')
data4 <- data4[data4$gene %in% coverage.genes ]
write.table(data4, quote=F, row.names = T, col.names = T, file = 'markers_ASDvsCTRL_Inhibitory_filtered', sep = '\t')

data5 <- read.table('markers_ASDvsCTRL_Maturing', header=T, sep='\t')
data5 <- data1[data5$gene %in% coverage.genes ]
write.table(data5, quote=F, row.names = T, col.names = T, file = 'markers_ASDvsCTRL_Maturing_filtered', sep = '\t')

data6 <- read.table('markers_ASDvsCTRL_Microglia', header=T, sep='\t')
data6 <- data1[data6$gene %in% coverage.genes ]
write.table(data6, quote=F, row.names = T, col.names = T, file = 'markers_ASDvsCTRL_Microglia_filtered', sep = '\t')

data7 <- read.table('markers_ASDvsCTRL_Oligodendrocytes', header=T, sep='\t')
data7 <- data7[data7$gene %in% coverage.genes ]
write.table(data7, quote=F, row.names = T, col.names = T, file = 'markers_ASDvsCTRL_Oligodendrocytes_filtered', sep = '\t')

data8 <- read.table('markers_ASDvsCTRL_OPC', header=T, sep='\t')
data8 <- data8[data8$gene %in% coverage.genes ]
write.table(data8, quote=F, row.names = T, col.names = T, file = 'markers_ASDvsCTRL_OPC_filtered', sep = '\t')



#Venn diagram, get overlap of celltype markers NET vs DEG 
markers.all <- read.table('coding_gene_filtered_DEGs/markers_MergedCelltype_logcounts_CodingGenes.tsv', header=T, sep = '\t')
marker.df.list <- split.data.frame(markers.all, markers.all$cluster)

#get merged Celltype network hubs
lls.rank.df <- read.table('/home3/junhacha/HumanNetv3/Velemeshev2019/merged_celltype/Brain_LLS_strength_rank.tsv', header = T, sep = '\t')

#make list of hubs for each cellype in descending strength
strength.list <- list()
for (i in 1:length(colnames(lls.rank.df))){
  strength.list[[i]] <- rownames(lls.rank.df)[order(lls.rank.df[,i], decreasing = T)]
}

names(strength.list) <- colnames(lls.rank.df)
lapply(strength.list, head)


#no Others Neurons have DEGs so lets remove them from the list
strength.list[['Others']] <- NULL


#for each number of signifcant DEGs get the same number of Network Hub and draw venndiagram
library(ggvenn)
library(gridExtra)

gene.list <- list()
plot.list <- list()
for (i in 1:length(strength.list)){
  deg.genes <- as.character(marker.df.list[[i]]$gene)
  net.genes <- as.character(strength.list[[i]][1:length(deg.genes)])
  gene.list[[i]] <- list(LLS_top=net.genes, DEG_pos=deg.genes)
  names(gene.list)[[i]] <- names(marker.df.list)[[i]]
  
  p <- ggvenn(
    gene.list[[i]], columns = c("LLS_top", "DEG_pos"), #the names have to match the names of list 
    stroke_size = 0.5) +
    ggtitle(names(gene.list)[[i]])
  
  plot.list[[i]] <- p

}

ggsave(
  filename = '/home3/junhacha/HumanNetv3/Velemeshev2019/VennDiagram_celltype_CodingGeneFiltered.pdf', 
  plot = marrangeGrob(plot.list, nrow=3, ncol=3), 
  width = 12, height = 12
)
















