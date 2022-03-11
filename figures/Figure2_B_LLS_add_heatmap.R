library(igraph)
library(SCINET)
library(ACTIONet)

#load entire HNv3
data('HNv3_XC_LLS')

#add LLS to all networks
#################################
cancertype <- 'BC'
#################################

setwd(paste0('~/HumanNetv3/XC_nobatch/',cancertype,'/'))
Celltype.specific.networks <- readRDS(paste0('./',cancertype,'_celltypeNet_list.rds'))


#add LLS weight and sort each genepair alphabetically element is stored as edgelist
sorted.net.list <- SortAddLLS(Celltype.specific.networks, reference.network = graph.hn3)

saveRDS(sorted.net.list, paste0('./',cancertype,'_celltypeNet_list_LLS.rds'))



net.bc <- readRDS('~/HumanNetv3/XC_nobatch/BC/BC_celltypeNet_list_LLS.rds')
net.lc <- readRDS('~/HumanNetv3/XC_nobatch/LC/LC_celltypeNet_list_LLS.rds')
net.crc <- readRDS('~/HumanNetv3/XC_nobatch/CRC/CRC_celltypeNet_list_LLS.rds')
net.ovc <- readRDS('~/HumanNetv3/XC_nobatch/OvC/OvC_celltypeNet_list_LLS.rds')

net.list.cancer <- list(net.bc, net.lc, net.crc, net.ovc)
strength.list <- list()
cell.list <- list()
for (net in seq_along(net.list.cancer)){
  for (cell.net in (1:length(net.list.cancer[[net]]))){
    LLS.net <- net.list.cancer[[net]][[cell.net]][,c('1','2','LLS')]
    colnames(LLS.net) <- c('gene1', 'gene2', 'weight')
    net.graph <- graph_from_data_frame(LLS.net, directed=FALSE) 
    genes.strength <- strength(net.graph)[order(strength(net.graph), decreasing = T)]
    
    #remove ribo.genes
    genes.all <- V(net.graph)$name
    ribo.genes <-  genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
    genes.strength.ribo.r <- genes.strength[!(names(genes.strength) %in% ribo.genes)]
    
    #to normalize for network size, we use percentile rank
    final.rank <- dplyr::percent_rank(genes.strength.ribo.r)
    cell.list[[cell.net]] <- final.rank
  }
  names(cell.list) <- names(net.list.cancer[[net]])
  strength.list[[net]] <- cell.list
}

names(strength.list) <- c('BC','LC','CRC','OvC')
lapply(strength.list, names)

#remove NAs from each cancer type
strength.list.final = list()
for (cancer in names(strength.list)){
  cancer.strength <- strength.list[[cancer]]
  na.slot <- is.na(names(cancer.strength))
  strength.list.final[[cancer]] <- cancer.strength[!na.slot]
}

lapply(strength.list.final, names)

#make each element of list to dataframe
#draw heat map from strength.list.final
#combine each named vector to a 0 filled dataframe with unique gene name rows
library(purrr)
library(dplyr)
bc <- strength.list.final$BC
bc <- lapply(bc, data.frame)
bc <- bc %>%
  imap(~setNames(.x, .y)) %>%
  map(tibble::rownames_to_column) %>%
  reduce(full_join, by = "rowname") %>%
  mutate_all(tidyr::replace_na, 0)

crc <- strength.list.final$CRC
crc <- lapply(crc, data.frame)
crc <- crc %>%
  imap(~setNames(.x, .y)) %>%
  map(tibble::rownames_to_column) %>%
  reduce(full_join, by = "rowname") %>%
  mutate_all(tidyr::replace_na, 0)

lc <- strength.list.final$LC
lc <- lapply(lc, data.frame)
lc <- lc %>%
  imap(~setNames(.x, .y)) %>%
  map(tibble::rownames_to_column) %>%
  reduce(full_join, by = "rowname") %>%
  mutate_all(tidyr::replace_na, 0)

ovc <- strength.list.final$OvC
ovc <- lapply(ovc, data.frame)
ovc <- ovc %>%
  imap(~setNames(.x, .y)) %>%
  map(tibble::rownames_to_column) %>%
  reduce(full_join, by = "rowname") %>%
  mutate_all(tidyr::replace_na, 0)


##################################
cancer.set <- as.data.frame(ovc)
##################################

celltype <- c('T_cell', "B_cell", "Myeloid", "Fibroblast", "EC", "Cancer")
rownames(cancer.set) <- cancer.set$rowname
cancer.set$rowname <- NULL

#get top 10 high rank genes of immune cells
top.genes <- vector()
for (cell in c('T_cell', "B_cell", "Myeloid", "Fibroblast", "EC")){
  print(cell)
  perc.rank <- unlist(cancer.set[cell])
  names(perc.rank) <- rownames(cancer.set)
  perc.rank.top <- head(names(perc.rank[order(perc.rank, decreasing = T)]),15)
  top.genes <- union(top.genes, perc.rank.top) 
}

#subset data
data <- cancer.set[rownames(cancer.set) %in% top.genes ,colnames(cancer.set) %in% celltype]   
                       
#which heatmap looks pretty?
pheatmap::pheatmap(data, scale = 'row', cluster_cols = F)

#manual ordering
###############################################################################
data <- data[,c('EC','Fibroblast','T_cell','B_cell','Myeloid','Cancer')] #bc 10
data <- data[,c('Fibroblast','EC','T_cell','B_cell','Myeloid','Cancer')] #lc 15
data <- data[,c('B_cell','Myeloid','T_cell','Fibroblast','EC','Cancer')] #crc 15
data <- data[,c('Fibroblast','EC','T_cell','B_cell','Myeloid','Cancer')] #ovc 15

data <- data[,c('Fibroblast','EC','T_cell','B_cell','Myeloid','Cancer')] #bc 15


pheatmap::pheatmap(data, scale = 'row', cluster_cols = F)

########################################################################

pdf('~/HumanNetv3/Graphs/OvC_heatmap_LLSstrength_15.pdf', 7,14)
pheatmap::pheatmap(data, scale = 'row', cluster_cols = F)
dev.off()
dev.off()
#######################################################################


#write percentile rank of each cancer
setwd('~/HumanNetv3/table/')
write.table(bc, quote = F, col.names = T, row.names = F, sep = '\t',file = './BC_LLS_strength_rank.tsv')
write.table(lc, quote = F, col.names = T, row.names = F, sep = '\t',file = './LC_LLS_strength_rank.tsv')
write.table(ovc, quote = F, col.names = T, row.names = F, sep = '\t',file = './OvC_LLS_strength_rank.tsv')
write.table(crc, quote = F, col.names = T, row.names = F, sep = '\t',file = './CRC_LLS_strength_rank.tsv')












