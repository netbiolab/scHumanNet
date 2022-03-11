#this code sees if upregulated DEGs and increased centrality can be comparable and draws Venndiagram for each

#there are very small DEGs...lets just see if these gene are actually in the non zero pos or neg diff ranked genes of networks
df.rank.cond <- read.table('/home3/junhacha/HumanNetv3/Velemeshev2019/condition_merged_celltype/Diff_rank_genes.tsv', header = T, sep = '\t')


#get differential hub genes for each celltype
strength.list.cond <- list()
head(df.rank.cond)
for (i in seq(1,16,2)){print(i)
  genes <- as.character(df.rank.cond[,i])
  diff.rank <- as.numeric(df.rank.cond[,i+1])
  
  strength.list.cond[[(i+1)/2]] <- diff.rank
  names(strength.list.cond[[(i+1)/2]]) <- genes
  names(strength.list.cond)[[(i+1)/2]] <- colnames(df.rank.cond)[i]
}
lapply(strength.list.cond, head)



#get logfold change value of ASD vs CTRL for normalized single cell and plot for each celltype
seurat <- readRDS('/home3/junhacha/HumanNetv3/Velemeshev2019/seurat_all.rds')
seurat@active.ident <- as.factor(seurat$celltypes_merged)


plot.list <- list()
diffK.logFC.list <- list()
for (i in seq(1,8,1)){
  celltype <- levels(seurat@active.ident)[[i]]
  if (celltype == 'Oligodendrocytes'){cell.name <- 'Oligo'}
  else {cell.name <- celltype}
  seurat.celltype <- seurat[,seurat@active.ident == celltype] 
  data <- FoldChange(seurat.celltype, ident.1 = colnames(seurat.celltype[,seurat.celltype$diagnosis == 'ASD']), 
                     ident.2 = colnames(seurat.celltype[,seurat.celltype$diagnosis == 'Control']), base=2)
  
  data.f <- data[rownames(data) %in% df.rank.cond[,cell.name],] 
  df.rank.cond.celltype <- df.rank.cond[, c(which(colnames(df.rank.cond) == cell.name),which(colnames(df.rank.cond) == cell.name) +1)]
  rownames(df.rank.cond.celltype) <- df.rank.cond.celltype[,1]
  data.f$diffK <- df.rank.cond.celltype[rownames(data.f),2]
  
  diffK.logFC.list[[i]] <- data.f
  names(diffK.logFC.list)[[i]] <- celltype
  
  p <- ggplot(data.f, aes(x=diffK, y=avg_log2FC)) + 
    geom_point() + 
    ggtitle(celltype) +
    theme(text=element_text(size=20))
    
  plot.list[[i]] <- p
  
}

#draw plot
library(gridExtra)
ggsave(
  filename = '/home3/junhacha/HumanNetv3/Velemeshev2019/DiffK_LogFC.pdf', 
  plot = marrangeGrob(plot.list, nrow=3, ncol=3), 
  width = 12, height = 12
)

saveRDS(diffK.logFC.list, '/home3/junhacha/HumanNetv3/Velemeshev2019/celltype_diffK_logFC.rds')




#divide strength.list.cond to pos and neg
#we probabaly want a threshold for each celltype differential hubs...and it probably wont be fair to put this in comparision with DEGs
#lets say that we want at least 30 percentile rank differences to say they are different
strength.list.cond.pos <- lapply(strength.list.cond, function(x){x[x>0]})
strength.list.cond.neg <- lapply(strength.list.cond, function(x){x[x<0]})

lapply(strength.list.cond.pos, length)
lapply(strength.list.cond.neg, length)


#For DEGs we need to impose a stricter pvalue because many of the differntial nodes are 0..!
setwd('/home3/junhacha/HumanNetv3/Velemeshev2019/DEG')
data1 <- read.table('markers_ASDvsCTRL_Astrocyte', header = T, sep = '\t')
data2 <- read.table('markers_ASDvsCTRL_Endothelial', header = T, sep = '\t')
data3 <- read.table('markers_ASDvsCTRL_Excitatory', header = T, sep = '\t')
data4 <- read.table('markers_ASDvsCTRL_Inhibitory', header = T, sep = '\t')
data5 <- read.table('markers_ASDvsCTRL_Maturing', header = T, sep = '\t')
data6 <- read.table('markers_ASDvsCTRL_Microglia', header = T, sep = '\t')
data7 <- read.table('markers_ASDvsCTRL_Oligodendrocytes', header = T, sep = '\t')
data8 <- read.table('markers_ASDvsCTRL_OPC', header = T, sep = '\t')

marker.list.cond <- list(data1,data2,data3,data4,data5,data6,data7,data8)
names(marker.list.cond) <- names(strength.list.cond) #i ordered it as the same as NET hub list

#divde by logfold change and get only logfold data
marker.list.cond.up <- lapply(marker.list.cond, function(x){
  up <- x[x$avg_logFC > 0,]
  up.ordered <- up[order(up$avg_logFC, decreasing = T),]
  up.ordered <- up.ordered$avg_logFC
  names(up.ordered) <- rownames(up)
  return(up.ordered)
})

marker.list.cond.down <- lapply(marker.list.cond, function(x){
  down <- x[x$avg_logFC < 0,]
  down.ordered <- down[order(down$avg_logFC, decreasing = F),]
  down.ordered <- down.ordered$avg_logFC
  names(down.ordered) <- rownames(down)
  return(down.ordered)
})

lapply(marker.list.cond.up, length)
lapply(marker.list.cond.down, length)


#lets see how many of DEG genes are in either pos or neg
for (i in 1:length(marker.list.cond.up)){
  deg.genes <- names(marker.list.cond.up[[i]])
  net.genes.pos <- names(strength.list.cond.pos[[i]])
  net.genes.neg <- names(strength.list.cond.neg[[i]])
  print(paste(names(marker.list.cond.up)[[i]], 'upDEGs in increased Hub:'))
  print(table(deg.genes %in% net.genes.pos))
  print(paste(names(marker.list.cond.up)[[i]], 'upDEGs in decreased Hub:'))
  print(table(deg.genes %in% net.genes.neg))
}


#it seems that up regulated degs are usually in increased hubs..but not all degs
for (i in 1:length(marker.list.cond.down)){
  deg.genes <- names(marker.list.cond.down[[i]])
  net.genes.pos <- names(strength.list.cond.pos[[i]])
  net.genes.neg <- names(strength.list.cond.neg[[i]])
  print(paste(names(marker.list.cond.up)[[i]], 'downDEGs in decreased Hub:'))
  print(table(deg.genes %in% net.genes.neg))
  print(paste(names(marker.list.cond.up)[[i]], 'downDEGs in increased Hub:'))
  print(table(deg.genes %in% net.genes.pos))
}



#perform fisher test to see upregulated DEGs and increased centraligy genes are associated
#make a contingency up  table
all.genes <- as.character(df.rank.cond$Astrocyote) #all comlumn contains every coding genes
for (i in 1:length(marker.list.cond.down)){
  deg.genes.up <- names(marker.list.cond.up[[i]])
  deg.genes.up.NOT <- all.genes[!(all.genes %in% deg.genes.up)]
  net.genes.pos <-  names(strength.list.cond.pos[[i]])
  net.genes.pos.NOT <- all.genes[!(all.genes %in% net.genes.pos)]
  
  conting <- matrix(c(length(intersect(deg.genes.up, net.genes.pos)),
                      length(intersect(deg.genes.down, net.genes.pos)),
                      length(intersect(deg.genes.up, net.genes.neg)),
                      length(intersect(deg.genes.down, net.genes.neg))),
                    nrow=2, dimnames = list(c("posNET", "!posNET"),
                                            c("upDEG", "!upDEG")))
  
  celltype <- names(marker.list.cond.up)[[i]]                  
  pvalue <- fisher.test(conting, alternative = "two.sided")$p.value
  print(paste(celltype, 'Fisher:', pvalue))
  print(conting)
}

#there seems to be no real associations..and this is expected
#because DEG and NET centrality offer different information

#how many DEGs are there in whole again?
lapply(marker.list.cond.down, length)
lapply(marker.list.cond.up, length)


#should we remove maturing...it seems sketchy??
#strength.list.cond[['Maturing']] <- NULL
#marker.list.cond[['Maturing']] <- NULL


#since there is no realy association with upDEG and pos NET we will not divide them to draw venn diagram


#check the target genes for several genes

#read scNET for each condition
#do target genes contatin SFARI genes?
sfari.genes <- read.table('/home3/junhacha/HumanNetv3/Velemeshev2019/SFARI-Gene_genes_09-02-2021release_01-10-2022export.csv', header = T, sep = ',')
sfari.genes <- as.character(sfari.genes$gene.symbol)


setwd('/home3/junhacha/HumanNetv3/Velemeshev2019/condition_merged_celltype/')
net.list <- readRDS('sorted_el_list.rds')

net.inh.asd <- net.list[['ASD_Inhibitory']]
net.inh.ctrl <- net.list[['Control_Inhibitory']]

net.inh.asd$scinet_weight <- NULL
net.inh.ctrl$scinet_weight <- NULL


#select hub gene
gene <- 'GRIN2B'
gene <- 'CACNA1A'
gene <- 'MECP2'


#filter to network that has the gene
net.sub.asd <- net.inh.asd[net.inh.asd[,1] %in% gene | net.inh.asd[,2] %in% gene,]
net.sub.ctrl <- net.inh.ctrl[net.inh.ctrl[,1] %in% gene | net.inh.ctrl[,2] %in% gene,]

genes.asd <- union(net.sub.asd[,1], net.sub.asd[,2])
asd.target <- genes.asd[!(genes.asd %in% gene)]

genes.ctrl <- union(net.sub.ctrl[,1], net.sub.ctrl[,2])
ctrl.target <- genes.ctrl[!(genes.ctrl %in% gene)]

#write network genes for functional enrichment gene
write.table(genes.ctrl, quote=F,row.names = F, col.names = F, file = '/home3/junhacha/HumanNetv3/Velemeshev2019/CACMA1A_ctrlGenes.tsv')
write.table(genes.asd, quote=F,row.names = F, col.names = F, file = '/home3/junhacha/HumanNetv3/Velemeshev2019/GRIN2B_ASDGenes.tsv')
write.table(genes.asd, quote=F,row.names = F, col.names = F, file = '/home3/junhacha/HumanNetv3/Velemeshev2019/MECP2_ASDGenes.tsv')



#draw functional enrichment plot for ctrl genes of CACNA1A in ctrl inhibitory network
library(enrichR)
GSAplot <- function(genes, database, title, top.term){
  dbs <- listEnrichrDbs()
  dbs <- database
  enrichr <- enrichr(genes, dbs)
  data.gsa <- enrichr[[dbs]]
  
  aaa <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',1))
  bbb <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',2))
  
  #add column
  data.gsa$overlap_num <- aaa / bbb
  data.gsa$log10_qvalue <- -log10(data.gsa$Adjusted.P.value) 
  
  data.gsa.f <- data.gsa[order(data.gsa$log10_qvalue, decreasing = T),]
  
  p <- ggplot(data.gsa.f[1:top.term,], aes(x = reorder(Term,log10_qvalue),log10_qvalue,  y = log10_qvalue,
                                           fill = overlap_num)) + 
    geom_bar(stat = 'identity', width = 0.9, position = position_dodge(width = 0.1)) + 
    geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'longdash') +
    coord_flip() + 
    theme(text = element_text(size = 14)) +
    scale_fill_gradient(name = 'overlap percentage',low = '#a1c4fd', high = '#ff9a9e') +
    theme_classic() +
    labs(
      title = title,
      y = '-log10(qvalue)',
      x = database
    )
  
  return(p)
}

#draw GSA of top n genes in rank df
pdf('/home3/junhacha/HumanNetv3/Graphs/GOBP_Brain_Inhib_CACNA1A_Ctrl.pdf',height=8,width=12)
GSAplot(genes.ctrl,'GO_Biological_Process_2021', 'CACNA1A Sub-network in Inhibitory neuron',20)
dev.off()

pdf('/home3/junhacha/HumanNetv3/Graphs/GOMF_Brain_Inhib_CACNA1A_Ctrl.pdf',height=8,width=12)
GSAplot(genes.ctrl,"GO_Molecular_Function_2021", 'CACNA1A Subnetwork in Inhibitory neuron',20)
dev.off()

GSAplot(top.genes,'KEGG_2021_Human', paste('GGI97 neigbor centrality top',length(top.genes),  'gene KEGG'),20)





#draw the two network
library(igraph)

#label Sfari Genes
ctrl.genes <- unique(c(net.sub.ctrl[,1],net.sub.ctrl[,2]))
asd.genes <- unique(c(net.sub.asd[,1],net.sub.asd[,2]))



#get subnet that contatins gene for Ctrl and ASD
inhib.net.sub.ctrl <- graph.data.frame(net.sub.ctrl, directed=F)
inhib.net.sub.asd <- graph.data.frame(net.sub.asd, directed = F)


#make SFARI genes label 
all.genes.asd <- data.frame(genes = asd.genes)
all.genes.asd$sfari.label <- ifelse(all.genes.asd$genes %in% sfari.genes, 1, 0)

all.genes.ctrl <- data.frame(genes = ctrl.genes)
all.genes.ctrl$sfari.label <- ifelse(all.genes.ctrl$genes %in% sfari.genes, 1, 0)


inhib.net.sub.ctrl <- set_vertex_attr(inhib.net.sub.ctrl, "SFARI", index = V(inhib.net.sub.ctrl), as.character(all.genes.ctrl$sfari.label))
inhib.net.sub.asd <- set_vertex_attr(inhib.net.sub.asd, "SFARI", index = V(inhib.net.sub.asd), as.character(all.genes.asd$sfari.label))



#color
V(inhib.net.sub.ctrl)$color = rgb(0.8,0.4,0.3,0.8)
V(inhib.net.sub.ctrl)$color = ifelse(V(inhib.net.sub.ctrl)$SFARI == '1', "red", V(inhib.net.sub.ctrl)$color)

V(inhib.net.sub.asd)$color = rgb(0.8,0.4,0.3,0.8)
V(inhib.net.sub.asd)$color = ifelse(V(inhib.net.sub.asd)$SFARI == '1', "red", V(inhib.net.sub.asd)$color)

#make GGI97 genes bigger in size
V(inhib.net.sub.ctrl)$size <- V(inhib.net.sub.ctrl)$SFARI
V(inhib.net.sub.asd)$size <- V(inhib.net.sub.asd)$SFARI

normalize_01 <- function(x) (x - min(x)) / (max(x) - min(x)) + 1.2

V(inhib.net.sub.ctrl)$size <- normalize_01(as.numeric(V(inhib.net.sub.ctrl)$SFARI)) * 5
V(inhib.net.sub.asd)$size <- normalize_01(as.numeric(V(inhib.net.sub.asd)$SFARI)) * 5



#visualize ctrl & ASD side by side
par(mfrow=c(1,2))
#visualize subnet
plot.igraph(inhib.net.sub.ctrl, asp = 0, main = paste("Inh Neu Subnet with", gene,'in Ctrl'),
                  ## colors =======================================
                  vertex.frame.color = "white",             ## border color
                  ## shapes =======================================
                  vertex.shape = "circle",                  ## none, circle, square, csquare, 
                  ## vrectangle, pie, raster, sphere
                  ## rectangle, crectangle
                  ## sizes =======================================
                  vertex.size = 10,                         ## size, default = 15
                  vertex.size2 = NA,                      ## second dimension size (for parallelograms)
                  ## color =======================================
                  vertex.label.color = "black",
                  ## font family =======================================
                  vertex.label.family = "Times",
                  ## font face =======================================
                  vertex.label.font = 2, ## 1 = plain, 2 = bold, 3 = italic, 4 = bold/italic
                  vertex.label.cex = normalize_01(as.numeric(V(inhib.net.sub.ctrl)$SFARI)) * 0.5
                  
)

#visualize subnet
plot.igraph(inhib.net.sub.asd, asp = 0, main = paste("Inh Neu Subnet with", gene, 'in ASD'),
                  ## colors =======================================     
                  vertex.frame.color = "white",             ## border color
                  ## shapes =======================================
                  vertex.shape = "circle",                  ## none, circle, square, csquare, 
                  ## vrectangle, pie, raster, sphere
                  ## rectangle, crectangle
                  ## sizes =======================================
                  vertex.size = 10,                         ## size, default = 15
                  vertex.size2 = NA,                      ## second dimension size (for parallelograms)
                  ## color =======================================
                  vertex.label.color = "black",
                  ## font family =======================================
                  vertex.label.family = "Times",
                  ## font face =======================================
                  vertex.label.font = 2, ## 1 = plain, 2 = bold, 3 = italic, 4 = bold/italic
                  vertex.label.cex = normalize_01(as.numeric(V(inhib.net.sub.asd)$SFARI)) * 0.5
                  
)













