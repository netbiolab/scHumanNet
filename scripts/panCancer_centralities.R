library(igraph)

net.bc <- readRDS('~/HumanNetv3/XC_nobatch/BC/BC_celltypeNet_list_LLS.rds')
net.lc <- readRDS('~/HumanNetv3/XC_nobatch/LC/LC_celltypeNet_list_LLS.rds')
net.crc <- readRDS('~/HumanNetv3/XC_nobatch/CRC/CRC_celltypeNet_list_LLS.rds')
net.ovc <- readRDS('~/HumanNetv3/XC_nobatch/OvC/OvC_celltypeNet_list_LLS.rds')


get_centrality <- function(method = method){
  
  net.list.cancer <- list(BC=net.bc, LC=net.lc, CRC=net.crc, OvC=net.ovc)
  cent.list <- list()
  cell.list <- list()
  for (net in seq_along(net.list.cancer)){
    for (cell.net in (1:length(net.list.cancer[[net]]))){
      LLS.net <- net.list.cancer[[net]][[cell.net]][,c('1','2','LLS')]
      colnames(LLS.net) <- c('gene1', 'gene2', 'weight')
      net.graph <- graph_from_data_frame(LLS.net, directed=FALSE) 
      
      if (method == 'degree'){
      genes.cent <- strength(net.graph)[order(strength(net.graph), decreasing = T)]
      }
      else if (method == 'betweenness'){
        genes.cent <- betweenness(net.graph)[order(betweenness(net.graph), decreasing = T)]
      }
      else if (method == 'closeness'){
        genes.cent <- closeness(net.graph)[order(closeness(net.graph), decreasing = T)]
      }
      else if (method == 'eigenvector'){
        eigen.cent <- eigen_centrality(net.graph, weights = E(net.graph)$weight)$vector
        genes.cent <- eigen.cent[order(eigen.cent, decreasing = T)] # change this measure
      }
      else{
        print('wrong method input')
      }
      
      #remove ribo.genes
      genes.all <- V(net.graph)$name
      ribo.genes <-  genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
      genes.cent.ribo.r <- genes.cent[!(names(genes.cent) %in% ribo.genes)]
      
      #to normalize for network size, we use percentile rank
      final.rank <- dplyr::percent_rank(genes.cent.ribo.r)
      cell.list[[cell.net]] <- final.rank
    }
    names(cell.list) <- names(net.list.cancer[[net]])
    cent.list[[net]] <- cell.list
  }
  names(cent.list) <- names(net.list.cancer)
  
  #remove NAs from each cancer type
  cent.list.final = list()
  for (cancer in names(cent.list)){
    cancer.cent <- cent.list[[cancer]]
    na.slot <- is.na(names(cancer.cent))
    cent.list.final[[cancer]] <- cancer.cent[!na.slot]
  }
  
  return(cent.list.final)
}
  

strength.list <- get_centrality(method = 'degree')
betweenness.list <- get_centrality(method = 'betweenness')
closeness.list <- get_centrality(method = 'closeness')
eigenvector.list <- get_centrality(method = 'eigenvector')


#make each into a table and save
saveRDS(strength.list, '~/HumanNetv3/data/strength_list_pancacer.rds')
saveRDS(betweenness.list, '~/HumanNetv3/data/betweenness_list_pancacer.rds')
saveRDS(closeness.list, '~/HumanNetv3/data/closeness_list_pancacer.rds')
saveRDS(eigenvector.list, '~/HumanNetv3/data/eigenvector_list_pancacer.rds')


#make each element of list to dataframe
#draw heat map from strength.list.final
#combine each named vector to a 0 filled dataframe with unique gene name rows
library(purrr)
library(dplyr)
coverage.genes <- readRDS('/home3/junhacha/HumanNetv3/Benchmark/HNv3_coverage_genes.rds')

combine.percrank <- function(perc.rank.list = perc.rank.list){
  df <- lapply(perc.rank.list, data.frame)
  df <- df %>%
    imap(~setNames(.x, .y)) %>%
    map(tibble::rownames_to_column) %>%
    reduce(full_join, by = "rowname") %>%
    mutate_all(tidyr::replace_na, 0)
  
  add.genes <- coverage.genes[!(coverage.genes %in% df$rowname)]
  
  add.df <- data.frame(matrix(0, nrow = length(add.genes), ncol = ncol(df)))
  add.df[,1] <- add.genes
  colnames(add.df) <- colnames(df)
  rank.df.final <- rbind(df, add.df)
  
  return(rank.df.final)
}


makedf <- function(cent.list.final = cent.list.final){
  
  bc <- combine.percrank(cent.list.final$BC)
  lc <- combine.percrank(cent.list.final$LC)
  crc <- combine.percrank(cent.list.final$CRC)
  ovc <- combine.percrank(cent.list.final$OvC)
  
  df.list <- list(BC=as.data.frame(bc), LC=as.data.frame(lc), CRC=as.data.frame(crc), OvC=as.data.frame(ovc))
  
  return(df.list)
}


strength.df <- makedf(cent.list.final = strength.list)
betweenness.df <- makedf(cent.list.final = betweenness.list)
closeness.df <- makedf(cent.list.final = closeness.list)
eigenvector.df <- makedf(cent.list.final = eigenvector.list)


#write each centrality to file
library(openxlsx)

centrality.df.list <- list( 
           'BC_strength' = centrality.df.list$strength$BC,
           'BC_betweenness' = centrality.df.list$betweenness$BC,
           'BC_closeness' = centrality.df.list$closeness$BC,
           'BC_eigenvector' = centrality.df.list$eigenvector$BC,
           
           'LC_strength' = centrality.df.list$strength$LC,
           'LC_betweenness' = centrality.df.list$betweenness$LC,
           'LC_closeness' = centrality.df.list$closeness$LC,
           'LC_eigenvector' = centrality.df.list$eigenvector$LC,
           
           'CRC_strength' = centrality.df.list$strength$CRC,
           'CRC_betweenness' = centrality.df.list$betweenness$CRC,
           'CRC_closeness' = centrality.df.list$closeness$CRC,
           'CRC_eigenvector' = centrality.df.list$eigenvector$CRC,
           
           'OvC_strength' = centrality.df.list$strength$OvC,
           'OvC_betweenness' = centrality.df.list$betweenness$OvC,
           'OvC_closeness' = centrality.df.list$closeness$OvC,
           'OvC_eigenvector' = centrality.df.list$eigenvector$OvC
)

write.xlsx(centrality.df.list, file = '/home3/junhacha/HumanNetv3/data/All_centrality_PanCancer.xlsx')


#celltype <- c('T_cell', "B_cell", "Myeloid", "Fibroblast", "EC", "Cancer")
#rownames(cancer.set) <- cancer.set$rowname
#cancer.set$rowname <- NULL
#
##get top 10 high rank genes of immune cells
#top.genes <- vector()
#for (cell in c('T_cell', "B_cell", "Myeloid", "Fibroblast", "EC")){
#  print(cell)
#  perc.rank <- unlist(cancer.set[cell])
#  names(perc.rank) <- rownames(cancer.set)
#  perc.rank.top <- head(names(perc.rank[order(perc.rank, decreasing = T)]),15)
#  top.genes <- union(top.genes, perc.rank.top) 
#}
#
##subset data
#data <- cancer.set[rownames(cancer.set) %in% top.genes ,colnames(cancer.set) %in% celltype]   
#
##which heatmap looks pretty?
#pheatmap::pheatmap(data, scale = 'row', cluster_cols = F)
#
##manual ordering
################################################################################
#data <- data[,c('EC','Fibroblast','T_cell','B_cell','Myeloid','Cancer')] #bc 10
#data <- data[,c('Fibroblast','EC','T_cell','B_cell','Myeloid','Cancer')] #lc 15
#data <- data[,c('B_cell','Myeloid','T_cell','Fibroblast','EC','Cancer')] #crc 15
#data <- data[,c('Fibroblast','EC','T_cell','B_cell','Myeloid','Cancer')] #ovc 15
#
#data <- data[,c('Fibroblast','EC','T_cell','B_cell','Myeloid','Cancer')] #bc 15
#
#
#pheatmap::pheatmap(data, scale = 'row', cluster_cols = F)
#
#########################################################################
#
#pdf('~/HumanNetv3/Graphs/OvC_heatmap_LLSstrength_15.pdf', 7,14)
#pheatmap::pheatmap(data, scale = 'row', cluster_cols = F)
#dev.off()
#dev.off()
########################################################################
#

#write percentile rank of each cancer
setwd('~/HumanNetv3/table/')
write.table(bc, quote = F, col.names = T, row.names = F, sep = '\t',file = './BC_LLS_strength_rank.tsv')
write.table(lc, quote = F, col.names = T, row.names = F, sep = '\t',file = './LC_LLS_strength_rank.tsv')
write.table(ovc, quote = F, col.names = T, row.names = F, sep = '\t',file = './OvC_LLS_strength_rank.tsv')
write.table(crc, quote = F, col.names = T, row.names = F, sep = '\t',file = './CRC_LLS_strength_rank.tsv')

