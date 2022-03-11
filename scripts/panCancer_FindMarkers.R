# this code reads top toplogically important Genes and and find the percentile rank of these genes in DEGs
library(Seurat)
library(future)


#save each data as seperate Seurat Object (version3)
count <- Read10X("~/public_Data/lambrechts_TME/10x_counts/OvC_counts/")
meta <- read.table('~/public_Data/lambrechts_TME/OvC_metadata.csv', sep = ',', header = T)
rownames(meta) <- meta$Cell
all.equal(colnames(count), rownames(meta))
meta$Cell <- NULL

seurat <- CreateSeuratObject(counts = count, meta.data = meta, project = 'OvC', min.cells = 0, min.features = 0)
saveRDS(seurat, '~/HumanNetv3/data/seurat_OvC.rds')


#perform DEG on log normalized counts (performs analysis on data slot)
plan()
plan('multiprocess', workers = 8)
options(future.globals.maxSize= 10485760000)

#for BC
###########################################################
seurat <- readRDS('~/HumanNetv3/data/seurat_BC.rds')
###########################################################

#log normalize and scale
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize") 
seurat <- ScaleData(seurat, features = rownames(seurat))

seurat@active.ident <- seurat$CellType
table(seurat$CellType)


markers <- FindAllMarkers(seurat)
markers <- markers[markers$p_val_adj < 0.05,]

write.table(markers, file='/home3/junhacha/HumanNetv3/data/DEGlist_BC_byCelltype_all_logcounts.tsv', 
            sep = '\t' ,quote = F, col.names = T, row.names = T)


#for LC
###########################################################
seurat <- readRDS('~/HumanNetv3/data/seurat_LC.rds')
###########################################################

#log normalize and scale
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize") 
seurat <- ScaleData(seurat, features = rownames(seurat))

seurat@active.ident <- seurat$CellType
table(seurat$CellType)

markers <- FindAllMarkers(seurat)
markers <- markers[markers$p_val_adj < 0.05,]

write.table(markers, file='/home3/junhacha/HumanNetv3/data/DEGlist_LC_byCelltype_all_logcounts.tsv', 
            sep = '\t' ,quote = F, col.names = T, row.names = T)


#for CRC
###########################################################
seurat <- readRDS('~/HumanNetv3/data/seurat_CRC.rds')
###########################################################

#log normalize and scale
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize") 
seurat <- ScaleData(seurat, features = rownames(seurat))

seurat@active.ident <- seurat$CellType
table(seurat$CellType)

markers <- FindAllMarkers(seurat)
markers <- markers[markers$p_val_adj < 0.05,]

write.table(markers, file='/home3/junhacha/HumanNetv3/data/DEGlist_CRC_byCelltype_all_logcounts.tsv', 
            sep = '\t' ,quote = F, col.names = T, row.names = T)




#for OvC
###########################################################
seurat <- readRDS('~/HumanNetv3/data/seurat_OvC.rds')
###########################################################

#log normalize and scale
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize") 
seurat <- ScaleData(seurat, features = rownames(seurat))

seurat@active.ident <- seurat$CellType
table(seurat$CellType)

markers <- FindAllMarkers(seurat)
markers <- markers[markers$p_val_adj < 0.05,]

write.table(markers, file='/home3/junhacha/HumanNetv3/data/DEGlist_OvC_byCelltype_all_logcounts.tsv', 
            sep = '\t' ,quote = F, col.names = T, row.names = T)



#draw grouped histogram of comparable perc rank values of top 100? genes
#read lls.rank data
setwd('~/HumanNetv3/table/')
bc.rank <- read.table('./BC_LLS_strength_rank.tsv', sep = '\t', header = T, stringsAsFactors = F)
lc.rank <- read.table('./LC_LLS_strength_rank.tsv', sep = '\t', header = T, stringsAsFactors = F)
ovc.rank <- read.table('./OvC_LLS_strength_rank.tsv', sep = '\t', header = T, stringsAsFactors = F)
crc.rank <- read.table('./CRC_LLS_strength_rank.tsv', sep = '\t', header = T, stringsAsFactors = F)


#bc
####################################################
colnames(bc.rank)
celltype <- colnames(bc.rank)[2:length(colnames(bc.rank))]
###################################################

#calculate perc.rank of LLS and DEG for each geneset 
test <- bc.rank[,c("rowname","T_cell")]
test <- test[order(test$T_cell, decreasing = T),]
head(test)

genes.to.test <- test$rowname[1:1000]
test1<- test[1:1000, ]



bc.deg <- read.table('~/HumanNetv3/data/DEGlist_BC_byCelltype_all_logcounts.tsv', sep = '\t', header = T)
deg.tcell <- bc.deg[bc.deg$cluster == 'T_cell',]

#sort by pvalue or logfold
deg.tcell$per.rank.deg.logfc <- dplyr::percent_rank(deg.tcell$avg_logFC)
deg.tcell$per.rank.deg.pval <- 1 - dplyr::percent_rank(deg.tcell$p_val_adj) #reverse percent rank

deg.tcell.pos <- deg.tcell[deg.tcell$avg_logFC > 0, ]

#draw venndigaram as example
#see overlap
library(VennDiagram)
legends = c("LLS_top1000", "DEG_pos")
###############################################
inter <- intersect(genes.to.test,as.character(deg.tcell.pos$gene))


p1 <-draw.pairwise.venn(
  area1 = length(genes.to.test),
  area2 = length(as.character(deg.tcell.pos$gene)),
  cross.area = length(inter),
  
  #colors
  fill = c('#E5B3BB', '#7B92AA'),
  #lty = "blank",
  col = 'black',
  
  #legends
  category = legends,
  cat.pos = c(0,0),
  ind=FALSE
)
grid.draw(p1)
grid.newpage()




test2 <- deg.tcell.pos[deg.tcell.pos$gene %in% genes.to.test, c("per.rank.deg.logfc","gene")]
colnames(test1) <- c('gene', 'per.rank.lls')

tmp <- dplyr::left_join(test1, test2, by='gene') %>% mutate_all(tidyr::replace_na, 0)


# library
library(ggplot2)
library(dplyr)
library(hrbrthemes)

# Build dataset with different distributions
data <- reshape2::melt(tmp)


# Represent it
data %>%
  ggplot( aes(x=value, fill=variable)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 0.1) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="") +
  ggtitle("Percentile Rank of top 1000 Network genes in DEG analysis")


head(tmp)

inter.data <- tmp[tmp$gene %in% inter,]

library(ggpubr)
ggpaired(inter.data, cond1 = 'per.rank.lls', cond2 = 'per.rank.deg.logfc', 
         fill='condition', palette = 'jco', line.color = 'gray')

genes.diff <- inter.data[abs(inter.data$per.rank.lls - inter.data$per.rank.deg) > 0.2 ,]


#make column based on absolute value change
#get interdata from melted df
inter.data.2 <- data[data$gene %in% inter.data$gene, ]

inter.data.2 %>%
  tidyr::spread(variable, value) %>%
  dplyr::mutate(is_diff = abs(per.rank.lls - per.rank.deg.logfc) > 0.2) %>%
  tidyr::gather("variable", "value", 2:3) %>%
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable), alpha = 0.5, col = "grey") +
  geom_point() +
  geom_line(aes(group = gene, col = is_diff)) +
  scale_colour_manual(values = c("gray", "red")) +
  theme(panel.background = element_blank())






