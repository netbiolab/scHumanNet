#this code reads GRNboost2 PIDC and SAVER output and compare them with HumanNet BC celltype network
library(igraph)
library(dplyr)
library(RColorBrewer)
library(UpSetR)
library(enrichR)
library(viridis)
library(grid)
library(ggplot2)

#Tcell
t.grnboost2 <- read.table('/home3/junhacha/HumanNetv3/Benchmark/GRNboost2/TcellNet_f_GRNboost2.tsv', header = T, sep = '\t')
t.saver <- read.table('/home3/junhacha/HumanNetv3/Benchmark/SAVER/Tcell_SV_LLS_cut.tsv', header = F, sep = '\t')
t.hnv3 <- read.table('/home3/junhacha/HumanNetv3/XC_nobatch/BC/net_T_cell_sort.tsv', header = F, sep = '\t')
t.metacell <- read.table('/home3/junhacha/HumanNetv3/Benchmark/MetaCell/Tcell_f/BC_metacell_Tcell_PCCnet_possorted_LLScut.tsv', header=F, sep = '\t')
t.rawPCC <- read.table('/home3/junhacha/HumanNetv3/Benchmark/rawPCC/Tcell_rawPCCnet0.8_possorted', header = F, sep = '\t')
t.bigscale <- read.table('/home3/junhacha/HumanNetv3/Benchmark/bigSCale2/T/Tcell_BCnet_BS_PCCnet_possorted_LLScut.tsv', sep = ' ')

#Bcell
b.grnboost2 <- read.table('/home3/junhacha/HumanNetv3/Benchmark/GRNboost2/BcellNet_f_GRNboost2.tsv', header = T, sep = '\t')
b.saver <- read.table('/home3/junhacha/HumanNetv3/Benchmark/SAVER/Bcell_SV_LLS_cut.tsv', header = F, sep = '\t')
b.hnv3 <- read.table('/home3/junhacha/HumanNetv3/XC_nobatch/BC/net_B_cell_sort.tsv', header = F, sep = '\t')
b.metacell <- read.table('/home3/junhacha/HumanNetv3/Benchmark/MetaCell/Bcell_f/BC_metacell_Bcell_PCCnet_possorted_LLScut.tsv', header=F, sep = ' ')
b.rawPCC <- read.table('/home3/junhacha/HumanNetv3/Benchmark/rawPCC/Bcell_rawPCCnet0.8_possorted', header = F, sep = '\t')
b.bigscale <- read.table('/home3/junhacha/HumanNetv3/Benchmark/bigSCale2/B/Bcell_BCnet_BS_PCCnet_possorted_LLScut.tsv', sep = ' ')


#Myeloid
mye.grnboost2 <- read.table('/home3/junhacha/HumanNetv3/Benchmark/GRNboost2/MyeNet_f_GRNboost2.tsv', header = T, sep = '\t')
mye.saver <- read.table('/home3/junhacha/HumanNetv3/Benchmark/SAVER/Mye_SV_LLS_cut.tsv', header = F, sep = '\t')
mye.hnv3 <- read.table('/home3/junhacha/HumanNetv3/XC_nobatch/BC/net_Myeloid_sort.tsv', header = F, sep = '\t')
mye.metacell <- read.table('/home3/junhacha/HumanNetv3/Benchmark/MetaCell/Mye_f/BC_metacell_Mye_PCCnet_possorted_LLScut.tsv', header=F, sep = ' ')
mye.rawPCC <- read.table('/home3/junhacha/HumanNetv3/Benchmark/rawPCC/Myeloid_rawPCCnet0.8_possorted', header = F, sep = '\t')
mye.bigscale <- read.table('/home3/junhacha/HumanNetv3/Benchmark/bigSCale2/Mye/Mye_BCnet_BS_PCCnet_possorted_LLScut.tsv', sep = ' ')


#CAF
caf.grnboost2 <- read.table('/home3/junhacha/HumanNetv3/Benchmark/GRNboost2/CAFNet_f_GRNboost2.tsv', header = T, sep = '\t')
caf.saver <- read.table('/home3/junhacha/HumanNetv3/Benchmark/SAVER/CAF_SV_LLS_cut.tsv', header = F, sep = '\t')
caf.hnv3 <- read.table('/home3/junhacha/HumanNetv3/XC_nobatch/BC/net_Fibroblast_sort.tsv', header = F, sep = '\t')
caf.metacell <- read.table('/home3/junhacha/HumanNetv3/Benchmark/MetaCell/CAF_f/BC_metacell_CAF_PCCnet_possorted_LLScut.tsv', header=F, sep = ' ')
caf.rawPCC <- read.table('/home3/junhacha/HumanNetv3/Benchmark/rawPCC/CAF_rawPCCnet0.8_possorted', header = F, sep = '\t')
caf.bigscale <- read.table('/home3/junhacha/HumanNetv3/Benchmark/bigSCale2/CAF/CAF_BCnet_BS_PCCnet_possorted_LLScut.tsv', sep = ' ')


#EC
ec.grnboost2 <- read.table('/home3/junhacha/HumanNetv3/Benchmark/GRNboost2/ECNet_f_GRNboost2.tsv', header = T, sep = '\t')
ec.saver <- read.table('/home3/junhacha/HumanNetv3/Benchmark/SAVER/EC_SV_LLS_cut.tsv', header = F, sep = '\t')
ec.hnv3 <- read.table('/home3/junhacha/HumanNetv3/XC_nobatch/BC/net_EC_sort.tsv', header = F, sep = '\t')
ec.metacell <- read.table('/home3/junhacha/HumanNetv3/Benchmark/MetaCell/EC_f/BC_metacell_EC_PCCnet_possorted_LLScut.tsv', header=F, sep = ' ')
ec.rawPCC <- read.table('/home3/junhacha/HumanNetv3/Benchmark/rawPCC/EC_rawPCCnet0.8_possorted', header = F, sep = '\t')
ec.bigscale <- read.table('/home3/junhacha/HumanNetv3/Benchmark/bigSCale2/EC/EC_BCnet_BS_PCCnet_possorted_LLScut.tsv', sep = ' ')



#cut GRNBoost2 network to top 0.1 percent..we can't consider all genes and all TFs, (rawPCC was cut to 0.8), Metacell and SAVER was cut by LLS
t.grnboost2 <- t.grnboost2[t.grnboost2$importance > quantile(t.grnboost2$importance, 0.999),]
b.grnboost2 <- b.grnboost2[b.grnboost2$importance > quantile(b.grnboost2$importance, 0.999),]
mye.grnboost2 <- mye.grnboost2[mye.grnboost2$importance > quantile(mye.grnboost2$importance, 0.999),]
ec.grnboost2 <- ec.grnboost2[ec.grnboost2$importance > quantile(ec.grnboost2$importance, 0.999),]
caf.grnboost2 <- caf.grnboost2[caf.grnboost2$importance > quantile(caf.grnboost2$importance, 0.999),]




#make into igrpah and calculate degree and rank top 100
net.list.t <- list(hnv3 = t.hnv3, bigSCale2 = t.bigscale, saver = t.saver, grnboost2 = t.grnboost2, metacell = t.metacell, rawPCC = t.rawPCC)
net.list.b <- list(hnv3 = b.hnv3, bigSCale2 = b.bigscale, saver = b.saver, grnboost2 = b.grnboost2, metacell = b.metacell, rawPCC = b.rawPCC)
net.list.mye <- list(hnv3 = mye.hnv3, bigSCale2 = mye.bigscale, saver = mye.saver, grnboost2 = mye.grnboost2, metacell = mye.metacell, rawPCC = mye.rawPCC)
net.list.caf <- list(hnv3 = caf.hnv3, bigSCale2 = caf.bigscale, saver = caf.saver, grnboost2 = caf.grnboost2, metacell = caf.metacell, rawPCC = caf.rawPCC)
net.list.ec <- list(hnv3 = ec.hnv3, bigSCale2 = ec.bigscale, saver = ec.saver, grnboost2 = ec.grnboost2, metacell = ec.metacell, rawPCC = ec.rawPCC)

net.list.all <- list(Tcell=net.list.t,Bcell=net.list.b, Mye=net.list.mye, CAF=net.list.caf, EC=net.list.ec)

#get top 100 percentile rank genes for each celltype each method
top.hub <-lapply(net.list.all, function(net.list) {
  lapply(net.list, function(net){
  net.g <- graph_from_data_frame(net, directed = F)
  net.cent <- strength(net.g, weights = E(net.g)$weights)
  net.cent.rank <- dplyr::percent_rank(net.cent)
  net.cent.rank.ordered <- net.cent.rank[order(net.cent.rank, decreasing = T)]
  gene.names <- names(net.cent.rank.ordered)
  return(gene.names[1:100])
})})


##draw venndiagram of top 100
#library(Vennerable)
#Vdeg.t <- Venn(top.hub[[1]])
#Vdeg.b <- Venn(top.hub[[2]])
#Vdeg.mye <- Venn(top.hub[[3]])
#Vdeg.caf <- Venn(top.hub[[4]])
#Vdeg.ec <- Venn(top.hub[[5]])
#
#p1 <- plot(Vdeg.t, doWeights = TRUE)
#p2 <- plot(Vdeg.b, doWeights = TRUE)
#p3 <- plot(Vdeg.mye, doWeights = TRUE)
#p4 <- plot(Vdeg.caf, doWeights = TRUE)
#p5 <- plot(Vdeg.ec, doWeights = TRUE)
#
#pdf('/home3/junhacha/HumanNetv3/Benchmark/method_overlap_top100.pdf',4,4)
#plot(Vdeg.t, doWeights = TRUE)
#plot(Vdeg.b, doWeights = TRUE)
#plot(Vdeg.mye, doWeights = TRUE)
#plot(Vdeg.caf, doWeights = TRUE)
#plot(Vdeg.ec, doWeights = TRUE)
#dev.off()
#
#write gene.list
aaa <- as.data.frame(sapply(top.hub[[1]], unlist))

write.table(aaa, '/home3/junhacha/HumanNetv3/Benchmark/Tcell_methods.tsv', col.names = T, row.names = F, quote=F, sep='\t')



#getupset plot for top 100 genes in terms of degree

#set color pallete for 6 method
n <- 6
col <- brewer.pal(n, 'Accent')

plot.list.top100 <- list()
for (i in 1:length(top.hub)){
  celltype <- names(top.hub)[[i]]
  
  top100.hnv3 <- top.hub[[celltype]][['hnv3']]
  top100.saver <-  top.hub[[celltype]][['saver']]
  top100.grnboost2 <- top.hub[[celltype]][['grnboost2']]
  top100.metacell <-  top.hub[[celltype]][['metacell']]
  top100.rawPCC <- top.hub[[celltype]][['rawPCC']]
  top100.bigscale2 <- top.hub[[celltype]][['bigSCale2']]
  
  #make method list for each celltype
  net.array <- list(top100.hnv3, top100.bigscale2, top100.saver, top100.grnboost2, top100.metacell, top100.rawPCC)
  
  #see overlap
  set1 <- top100.hnv3
  set2 <- top100.bigscale2
  set3 <- top100.saver
  set4 <- top100.metacell
  set5 <- top100.grnboost2
  set6 <- top100.rawPCC
  
  set.list.v <-list('HNv3'=set1,'bigSCale2' = set2,'SAVER'=set3,'METACELL'=set4,'GRNBOOST2'=set5, 'rawPCC'=set6)
  
  m.node <- set.list.v %>% lapply(table) %>% lapply(as.list) %>% 
    lapply(data.frame) %>% bind_rows()
  
  rownames(m.node) <- c('HNv3', 'bigSCale2', 'SAVER', 'METACELL', 'GRNBOOST2', 'rawPCC')
  m.node[is.na(m.node)]<-0
  
  node.df <- as.data.frame(t(m.node))
  
  
  #draw upsetR
  p <- upset(node.df, sets = c('HNv3', 'bigSCale2', 'SAVER', 'METACELL', 'GRNBOOST2', 'rawPCC'), order.by = 'freq' , 
             point.size = 2, 
             sets.bar.color = col ,
             set_size.show = F, 
             text.scale = 2,
             show.numbers = "yes")
  
  plot.list.top100[[i]] <- p
}
saveRDS(plot.list.top100, '/home3/junhacha/HumanNetv3/Benchmark/node_upsetR_list_top100.rds')

#arrange upset plots for Node, T,B,Mye,CAF,EC...arrange doesn't work...T-T lets do it in illustrator
pdf('/home3/junhacha/HumanNetv3/Benchmark/upsetR_node_overlap_5methods_top100.pdf',6,6)
for (i in 1:length(plot.list.top100)){print(plot.list.top100[[i]])}
dev.off()




graph.list <-lapply(net.list.all, function(net.list) {
  lapply(net.list, function(net){
    net.g <- graph_from_data_frame(net, directed = F)
    names(net.g) <- names(net.list.all)
    return(net.g)
  })})


#get upsetR plot for all nodes and all edges
plot.list <- list()
for (i in 1:length(graph.list)){
  celltype <- names(graph.list)[[i]]
  hnv3 <- graph.list[[celltype]][['hnv3']]
  saver <-  graph.list[[celltype]][['saver']]
  grnboost2 <- graph.list[[celltype]][['grnboost2']]
  metacell <-  graph.list[[celltype]][['metacell']]
  rawPCC <- graph.list[[celltype]][['rawPCC']]
  bigscale2 <- graph.list[[celltype]][['bigSCale2']]
  
  #make method list for each celltype
  net.array <- list(hnv3, bigscale2, saver, metacell, grnboost2, rawPCC)
  
  #nodes overlap
  nodes.hnv3 <- V(hnv3)$name
  nodes.saver <- V(saver)$name
  nodes.grnboost2 <- V(grnboost2)$name
  nodes.metacell <- V(metacell)$name
  nodes.rawPCC <- V(rawPCC)$name
  nodes.bigscale2 <- V(bigscale2)$name
  
  #see overlap
  set1 <- nodes.hnv3
  set2 <- nodes.saver
  set3 <- nodes.grnboost2
  set4 <- nodes.metacell
  set5 <- nodes.rawPCC
  set6 <- nodes.bigscale2
  
  set.list.v <-list('HNv3'=set1,'bigSCale2' = set6,'SAVER'=set2,'METACELL'=set4,'GRNBOOST2'=set3, 'rawPCC'=set5)
  
  m.node <- set.list.v %>% lapply(table) %>% lapply(as.list) %>% 
    lapply(data.frame) %>% bind_rows()
  
  rownames(m.node) <- c('HNv3', 'bigSCale2', 'SAVER', 'METACELL', 'GRNBOOST2', 'rawPCC')
  m.node[is.na(m.node)]<-0
  
  node.df <- as.data.frame(t(m.node))
  
  
  #draw upsetR
  p <- upset(node.df, sets = c('HNv3', 'bigSCale2', 'SAVER', 'METACELL', 'GRNBOOST2', 'rawPCC'), order.by = 'freq' , 
             point.size = 2, 
             sets.bar.color = col ,
             set_size.show = F, 
             text.scale = 2,
             show.numbers = F)

  plot.list[[i]] <- p
}
saveRDS(plot.list, '/home3/junhacha/HumanNetv3/Benchmark/node_upsetR_list.rds')


#arrange upset plots for Node, T,B,Mye,CAF,EC
pdf('/home3/junhacha/HumanNetv3/Benchmark/node_overlap_5methods.pdf',6,6)
for (i in 1:length(plot.list)){print(plot.list[[i]])}
dev.off()





#get upset plot for edges
plot.list.edge <- list()
for (i in 1:length(graph.list)){
  celltype <- names(graph.list)[[i]]
  print(celltype)
  
  hnv3 <- graph.list[[celltype]][['hnv3']]
  saver <-  graph.list[[celltype]][['saver']]
  grnboost2 <- graph.list[[celltype]][['grnboost2']]
  metacell <-  graph.list[[celltype]][['metacell']]
  rawPCC <- graph.list[[celltype]][['rawPCC']]
  bigscale2 <- graph.list[[celltype]][['bigSCale2']]
  
  #make method list for each celltype
  net.array <- list(hnv3, bigscale2, saver, metacell, grnboost2, rawPCC)
  
  #sort edges
  net.sorted <- lapply(net.array, function(net){ 
    method.net <- get.data.frame(net)
    net.sorted <- as.data.frame(t(apply(method.net,1,sort)))
    net.final <- net.sorted[,c(2,3,1)]
    return(net.final)
  })
  names(net.sorted) <- c('HNv3','bigSCale2', 'SAVER','METACELL', 'GRNBOOST2','rawPCC')
  #lapply(net.sorted, head)
  lapply(net.sorted, dim)
  
  #get sorted edges
  links.hnv3 <- paste(net.sorted[["HNv3"]][,1], net.sorted[["HNv3"]][,2], sep = '-')
  links.bigscale2 <- paste(net.sorted[['bigSCale2']][,1], net.sorted[['bigSCale2']][,2], sep = '-')
  links.saver <- paste(net.sorted[['SAVER']][,1], net.sorted[['SAVER']][,2], sep = '-')
  links.grnboost2 <- paste(net.sorted[['GRNBOOST2']][,1], net.sorted[['GRNBOOST2']][,2], sep = '-')
  links.metacell <- paste(net.sorted[['METACELL']][,1], net.sorted[['METACELL']][,2], sep = '-')
  links.rawPCC <- paste(net.sorted[['rawPCC']][,1], net.sorted[['rawPCC']][,2], sep = '-')

  #see overlap
  set1 <- links.hnv3
  set2 <- links.saver
  set3 <- links.grnboost2
  set4 <- links.metacell
  set5 <- links.rawPCC
  set6 <- links.bigscale2
  
  set.list.l <-list('HNv3'=set1,'bigSCale2' = set6,'SAVER'=set2,'METACELL'=set4,'GRNBOOST2'=set3, 'rawPCC'=set5)
  
  m.edge <- set.list.l %>% lapply(table) %>% lapply(as.list) %>% 
    lapply(data.frame) %>% bind_rows()
  
  rownames(m.edge) <- c('HNv3', 'bigSCale2', 'SAVER', 'METACELL', 'GRNBOOST2', 'rawPCC')
  m.edge[is.na(m.edge)]<-0
  
  edge.df <- as.data.frame(t(m.edge))

  
  
  #draw upsetR
  p <- upset(edge.df, sets = c('HNv3', 'bigSCale2', 'SAVER', 'METACELL', 'GRNBOOST2', 'rawPCC'), order.by = 'freq' , point.size = 2, 
             sets.bar.color = col ,
             set_size.show = F, 
             text.scale = 2,
             show.numbers = F)
  
  plot.list.edge[[i]] <- p
}
saveRDS(plot.list.edge, '/home3/junhacha/HumanNetv3/Benchmark/edge_upsetR_list.rds')


#arrange upset plots for Node, T,B,Mye,CAF,EC...arrange doesn't work...T-T lets do it in illustrator
pdf('/home3/junhacha/HumanNetv3/Benchmark/edge_overlap_5methods.pdf',6,6)
for (i in 1:length(plot.list.edge)){print(plot.list.edge[[i]])}
dev.off()




#get enrichment overlap plot for azimuth for all celltypes top 100 nodes
top100.tcells <- as.data.frame(sapply(top.hub[[1]], unlist))
top100.bcells <- as.data.frame(sapply(top.hub[[2]], unlist))
top100.mye <- as.data.frame(sapply(top.hub[[3]], unlist))
top100.caf <- as.data.frame(sapply(top.hub[[4]], unlist))
top100.ec <- as.data.frame(sapply(top.hub[[5]], unlist))

#make into list of list
top100.df <- list(Tcells=top100.tcells, Bcells=top100.bcells, Myeloids=top100.mye, CAFs=top100.caf, ECs=top100.ec)


#get celltype for Tcell Bcells Myeloids
gsa.plot.list <- list()
for (i in 1:3){
  celltype = names(top100.df)[[i]]
  print(celltype)
  
  #GSA for each celltype each method
  dbs <- listEnrichrDbs()
  dbs <- "Azimuth_Cell_Types_2021"
  
  df.method <- top100.df[[celltype]]
  
  gsa.celltype <- list()
  for (k in 1:ncol(df.method)){
    method <- colnames(df.method)[[k]]
    genes <- as.character(df.method[,method])
    enrichr <- enrichr(genes, dbs)
  
    data.gsa <- enrichr$Azimuth_Cell_Types_2021
  
    aaa <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',1))
    bbb <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',2))
    
    #add column
    data.gsa$overlap_num <- aaa / bbb
    data.gsa$log10_qvalue <- -log10(data.gsa$Adjusted.P.value) 
    
    data.gsa.f <- data.gsa[order(data.gsa$log10_qvalue, decreasing = T),]
    
    #some have less than 15...
    if(nrow(data.gsa.f) < 15){
      data.gsa.f <- data.gsa.f
    }
    else{
      data.gsa.f <- data.gsa.f[1:15,]
    }
    
    p <- ggplot(data.gsa.f, aes(x = reorder(Term,log10_qvalue), log10_qvalue, y = log10_qvalue,
                                  fill = overlap_num)) + 
      geom_bar(stat = 'identity', width = 0.9, position = position_dodge(width = 0.1)) + 
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red", size=1) +
      coord_flip() + 
      theme(text = element_text(size = 20)) +
      #scale_fill_gradientn(colours = viridis(20), limits=c(0, 1)) +
      scale_fill_gradientn(colors = alpha(c(viridis(20)), alpha = .5), limits=c(0,1)) +
      labs(fill="Overlap\nPercentage") +
      theme_classic() +
      labs(
        title = paste("Top 100 genes in Azimuth DB:", celltype, method),
        y = '-log10(qvalue)',
        x = "Azimuth Cell Types 2021"
      )
    
    gsa.celltype[[k]] <- p
  }
  gsa.plot.list[[i]] <- gsa.celltype
}


#save and arrange plot
library(gridExtra)

#main figure is T B Mye for Hnv3
main.plot.list <- list(gsa.plot.list[[1]][[1]],gsa.plot.list[[2]][[1]],gsa.plot.list[[3]][[1]])

outputfolder <- '/home3/junhacha/HumanNetv3/Benchmark/'
ggsave(file = paste0(outputfolder,'HNV3_TBMye_Azimuth2021.pdf'), arrangeGrob(grobs = main.plot.list, ncol =3), width=16, height=4)  ## save plot


ggsave(file = paste0(outputfolder,'Tcells_6method_Azimuth2021.pdf'), arrangeGrob(grobs = gsa.plot.list[[1]][c(1:6)], ncol =3), width=24, height=10)  ## save plot
ggsave(file = paste0(outputfolder,'Bcells_6method_Azimuth2021.pdf'), arrangeGrob(grobs = gsa.plot.list[[2]][c(1:6)], ncol =3), width=24, height=10)  ## save plot
ggsave(file = paste0(outputfolder,'Myeloids_6method_Azimuth2021.pdf'), arrangeGrob(grobs = gsa.plot.list[[3]][c(1:6)], ncol =3), width=24, height=10)  ## save plot



##bar plot of housekeeping genes in the top 100
#hk.genes <- read.table('~/HumanNetv3/Benchmark/Human_House_Keeping_gene_list_1076.txt', header=F)
#hk.genes <- as.character(hk.genes[,1])
#
#
#all.nodes.5000 <-lapply(net.list.all, function(net.list) {
#  lapply(net.list, function(net){
#    net.g <- graph_from_data_frame(net, directed = F)
#    net.cent <- strength(net.g, weights = E(net.g)$weights)
#    net.cent.rank <- dplyr::percent_rank(net.cent)
#    net.cent.rank.ordered <- net.cent.rank[order(net.cent.rank, decreasing = T)]
#    gene.names <- names(net.cent.rank.ordered)
#    return(gene.names[1:5000])
#  })})
#
#
#lapply(all.nodes.5000, function(list){
#  lapply(list, function(genes){
#    table(genes %in% hk.genes)
#  })
#})
#
#
#hk.genes[grep('^RPS[0-9]|^RPL[0-9]', hk.genes)]
#
##what about ribosomal genes
#lapply(top.hub[[1]], function(genes){
#  ribo.genes <- genes[grep('^RPS[0-9]|^RPL[0-9]', genes)]
#  return(length(ribo.genes))
#})
#
#bad.genes <- var.features[grep('^RPS[0-9]|^RPL[0-9]|^TRAV[0-9]|^TRBV[0-9]|^IG[^(?!F)]', var.features)]
#
#
#