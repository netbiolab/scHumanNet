#conda activate HNv3
library(ACTIONet)
library(SCINET)
library(igraph)
library(SingleCellExperiment)
library(Seurat)


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



#convert to sce
data <- SingleCellExperiment(assays = list(logcounts = counts), colData = meta)

#if in seurat object
#data <- SingleCellExperiment(assays = list(logcounts = seurat@assays$RNA@counts), colData = seurat@meta.data)

#step 1 ACTIONet framework, annotate cells
#reduce data
ace = reduce.ace(data)

#input of run.Actionet is the reduced sce class
#ace = run.ACTIONet(data, thread_no=24) 


#mark cell annotation in the sce dataset
#*********************************************************
ace[['Labels']] <- meta$celltypes_merged
#*********************************************************


#check that the given celltype lanbels distinctively divides cells in the reduced dimensions
#plot.ACTIONet(ace, ace$Labels, trans_attr = ace$node_centrality)


#get celltype specificity
ace = compute.cluster.feature.specificity(ace, ace$Labels, "celltype_specificity_scores")

#run SCINET
Celltype.specific.networks = run.SCINET.clusters(ace, specificity.slot.name = "celltype_specificity_scores_feature_specificity")

#save the list of celltype specific HumaNetv3, each element in igraph object
saveRDS(Celltype.specific.networks, './celltypeNet_list.rds')


#load entire HNv3
data('HNv3_XC_LLS')

#add LLS weight and sort each genepair alphabetically element is stored as edgelist
sorted.net.list <- lapply(Celltype.specific.networks, function(cell.net){
  int.net <- intersection(graph.hn3, cell.net, keep.all.vertices=F)
  names(edge_attr(int.net)) <- c('LLS', 'scinet_weight')
  tmp <- t(apply(get.data.frame(int.net)[,1:2],1,sort))
  x <- cbind(tmp, get.data.frame(int.net)[,3:4])
})


#save net list, where each element is a sorted, LLS tagged edgelist for across network comparison
saveRDS(sorted.net.list, './sorted_el_list.rds')

#check if they have been made without error
lapply(sorted.net.list, head)

#This only works on R3 envrionment...dont ask why...who knows
#come back to Rstudio
#setwd
#********************************************************************
net.brain <- readRDS('sorted_el_list.rds')
#********************************************************************

#when calculating percentile ranks, Ribo genes MT genes were excluded from the list
strength.list <- list()
cell.list <- list()
for (cell.net in (1:length(net.brain))){
  LLS.net <- net.brain[[cell.net]][,c('1','2','LLS')]
  colnames(LLS.net) <- c('gene1', 'gene2', 'weight')
  net.graph <- graph_from_data_frame(LLS.net, directed=FALSE) 
  genes.strength <- strength(net.graph)[order(strength(net.graph), decreasing = T)]
  
  #remove ribo.genes and mito genes
  genes.all <- V(net.graph)$name
  ribo.genes <-  genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
  mito.genes <- genes.all[grep('^MRPS[0-9]*|^MRPL[0-9]*', genes.all)]
  genes.strength.ribo.r <- genes.strength[!(names(genes.strength) %in% c(ribo.genes, mito.genes))]
  
  #to normalize for network size, we use percentile rank
  final.rank <- dplyr::percent_rank(genes.strength.ribo.r)
  cell.list[[cell.net]] <- final.rank
}

strength.list <- cell.list
names(strength.list) <- names(net.brain)

#make each element of list to dataframe
#draw heat map from strength.list
#combine each named vector to a 0 filled dataframe with unique gene name rows
#Somehow...this works in R3 version...not HNv3 R4
library(purrr)
library(dplyr)
brain <- strength.list
brain <- lapply(brain, data.frame)
brain <- brain %>%
  imap(~setNames(.x, .y)) %>%
  map(tibble::rownames_to_column) %>%
  reduce(full_join, by = "rowname") %>%
  mutate_all(tidyr::replace_na, 0)


rownames(brain) <- brain$rowname

#add rest of the genes with 0 values
coverage.genes <- readRDS('/home3/junhacha/HumanNetv3/Benchmark/HNv3_coverage_genes.rds')
add.genes <- coverage.genes[!(coverage.genes %in% brain$rowname)]

add.df <- data.frame(matrix(0, nrow = length(add.genes), ncol = ncol(brain)))
add.df[,1] <- add.genes
colnames(add.df) <- colnames(brain)
rank.df.final <- rbind(brain, add.df)

#write brain centrality to file
write.table(rank.df.final, quote = F, col.names = T, row.names = F, sep = '\t',file = './Brain_LLS_strength_rank.tsv')


rownames(rank.df.final) <- rank.df.final$rowname
rank.df.final$rowname <- NULL

#get top 50 high rank genes for each brain celltype
top.genes <- list()
for (cell in colnames(rank.df.final)){
  perc.rank <- unlist(rank.df.final[cell])
  names(perc.rank) <- rownames(rank.df.final)
  perc.rank.top <- head(names(perc.rank[order(perc.rank, decreasing = T)]),50)
  top.genes[[cell]] <- perc.rank.top 
}

test <- data.frame(matrix(ncol = 1, nrow = 50))
for (genes in top.genes){
  test <- cbind(test,genes)
}
test[,1] <- NULL

#write top 50 genes for each merged celltype
colnames(test) <- names(rank.df.final)
write.table(test, 'top50_genes_Brain_celltype_MTexcluded.tsv', sep='\t',quote=F, row.names=F, col.names=T)


#draw heatmap for merged celltype

#get top 15 high rank genes of brain cells and plot in heatmap
top.genes <- vector()
for (cell in colnames(rank.df.final)){
  print(cell)
  perc.rank <- unlist(rank.df.final[cell])
  names(perc.rank) <- rownames(rank.df.final)
  perc.rank.top <- head(names(perc.rank[order(perc.rank, decreasing = T)]),15)
  top.genes <- union(top.genes, perc.rank.top) 
}



#subset data if needed (exxlude Others)
celltype <- c('Endothelial','Microglia','Astrocyte','Oligodendrocytes','OPC','Inhibitory', 'Excitatory')
data <- brain[rownames(brain) %in% top.genes, colnames(brain) %in% celltype]   

#remove rows with all zero values
data <- data[rowSums(data) !=0, ]

#which heatmap looks pretty?
pheatmap::pheatmap(data, scale = 'row', cluster_cols = T)

#manual ordering
###############################################################################
data <- data[,c('Endothelial','Oligodendrocytes','Astrocyte','Microglia','Inhibitory', 'OPC', 'Excitatory')] #bc top 15


pheatmap::pheatmap(data, scale = 'row', cluster_cols = F)

########################################################################
#save plot
pdf('~/HumanNetv3/Velemeshev2019/merged_celltype/Heatmap_top15_strength_genes.pdf', 7,14)
pheatmap::pheatmap(data, scale = 'row', cluster_cols = F)
dev.off()
dev.off()
#######################################################################



#for merged_celltype
colnames(test) <- colnames(brain)
write.table(test, 'top50_genes_Brain_celltype_MTexcluded.tsv', quote = F, col.names = T, row.names = F, sep = '\t')




#for each celltype out of top 50 genes find which gene has the biggest difference in strength
rank.list <- list()
for (celltype in names(table(meta$celltypes_merged))){
  print(celltype)
  #subset brain percentile rank
  df <- brain[,c(paste0('Control_', celltype), paste0('ASD_', celltype))]
  #its percentile rank...ASD-CTL lets higher value have higer centrality in ASD
  df$diff.rank <- df[,2] - df[,1]
  
  #make a list for each celltype get genes that have the most rank differential
  rank.list[[celltype]] <- df$diff.rank
  names(rank.list[[celltype]]) <- rownames(df)
  rank.list[[celltype]] <- rank.list[[celltype]][order(abs(rank.list[[celltype]]), decreasing = T)]
}

#check that all element of list contain same genes(all genes)

df.final <- data.frame(matrix(nrow=length(rank.list[[1]]), ncol=2*length(rank.list)))

for (i in seq(1,8)){
  t <- 2*i - 1
  df.final[,t] <- names(rank.list[[i]])
  df.final[,t+1] <- rank.list[[i]]
}

colnames(df.final) <- c('Astrocyte','Ast_ASD-Ctrl' ,'Endothelial','Endo_ASD-Ctrl','Excitatory','Exc_ASD-Ctrl', 'Inhibitory','Inh_ASD-Ctrl',
                        'Microglia','Mic_ASD-Ctrl','Oligo','Olg_ASD-Ctrl','OPC','OPC_ASD-Ctrl','Others','Others_ASD-Ctrl')

write.table(df.final, '/home3/junhacha/HumanNetv3/Velemeshev2019/condition_merged_celltype/Diff_rank_genes.tsv', sep='\t', quote=F, row.names = F, col.names = T)



#get only the Sfari genes DF
sfari.genes <- read.table('/home3/junhacha/HumanNetv3/Velemeshev2019/SFARI-Gene_genes_09-02-2021release_01-10-2022export.csv', header = T, sep = ',')
sfari.genes <- as.character(sfari.genes$gene.symbol)


#filter rank.list
rank.list.sfari <- lapply(rank.list, function(x){x[names(x) %in% sfari.genes]})

#check that all element of list contain same genes(all genes)
#out of 1023 download genes 936 were in the dataset
#which genes are not in the gene?
lapply(rank.list.sfari, length)
sfari.genes[!(sfari.genes %in% rownames(brain))]
df.final.sfari <- data.frame(matrix(nrow=length(rank.list.sfari[[1]]), ncol=2*length(rank.list.sfari)))
for (i in seq(1,8)){
  t <- 2*i - 1
  df.final.sfari[,t] <- names(rank.list.sfari[[i]])
  df.final.sfari[,t+1] <- rank.list.sfari[[i]]
}

colnames(df.final.sfari) <- colnames(df.final)

write.table(df.final.sfari, '/home3/junhacha/HumanNetv3/Velemeshev2019/condition_merged_celltype/Diff_rank_SFARIgenes.tsv', sep='\t', quote=F, row.names = F, col.names = T)




#test significance of difference....
test <- read.table('/home3/junhacha/HumanNetv3/Velemeshev2019/condition_merged_celltype/Diff_rank_genes.tsv', sep = '\t', header = T)

#first take the top 2.5 percent from each pos:high in ASD, neg: high in Ctrl

asd.highg.genes <- test$Inhibitory[which(test$Inh_ASD.Ctrl[test$Inh_ASD.Ctrl > quantile(test$Inh_ASD.Ctrl, 0.975)])]
hist(test$Inh_ASD.Ctrl, col = 'grey')





