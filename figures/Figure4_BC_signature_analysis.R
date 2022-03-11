#read signature data
sig.genes <- read.table('~/HumanNetv3/BC_signaturegenes_33union.csv', header = T, sep = ',')
sig.genes$X <- NULL
rownames(sig.genes) <- sig.genes$gene
sig.genes$gene <- NULL

bc.sig.list <- list()
for (i in seq_along(colnames(sig.genes))){
  sig.name <- colnames(sig.genes)[i]
  gene.list <- sig.genes[,colnames(sig.genes) == sig.name, drop=F]
  gene.list <- rownames(gene.list)[which(gene.list== 'yes')]
  bc.sig.list[[i]] <- gene.list
  names(bc.sig.list)[[i]] <- sig.name
}


#get number of genes for each signatures
lengths(bc.sig.list)

#load breast cancer networks
tnet <- read.table('~/HumanNetv3/table/BC_T_cell_net_LLS.tsv', sep = '\t', header = T)
bnet <- read.table('~/HumanNetv3/table/BC_B_cell_net_LLS.tsv', sep = '\t', header = T)
myenet <- read.table('~/HumanNetv3/table/BC_Myeloid_net_LLS.tsv', sep = '\t', header = T)
cafnet <- read.table('~/HumanNetv3/table/BC_Fibroblast_net_LLS.tsv', sep = '\t', header = T)
endonet <- read.table('~/HumanNetv3/table/BC_EC_net_LLS.tsv', sep = '\t', header = T)
cancernet <- read.table('~/HumanNetv3/table/BC_Cancer_net_LLS.tsv', sep = '\t', header = T)



#for each genesets, get connectivity for all bc networks
net.list <- list(tnet, bnet, myenet, cafnet, endonet, cancernet)
names(net.list) <- c('Tnet', 'Bnet', 'Myenet', 'CAFnet', 'ECnet', 'Cancernet')

#make gene node list for each network
node.list <- lapply(net.list, function(net){
  genes.all <- unique(c(as.character(net[,1]),as.character(net[,2])))
})
names(node.list) <- names(net.list)



sig.list <- list()
for (i in seq_along(bc.sig.list)){
  sig.genes <- bc.sig.list[[i]]
  connectivity.sig <- lapply(net.list, function(net){
    nrow(net[(net[,1] %in% sig.genes & net[,2] %in% sig.genes), ])
  })
  names(connectivity.sig) <- names(net.list)
  connectivity.sig <- dplyr::bind_rows(connectivity.sig)
  sig.list[[i]] <- connectivity.sig
}

connectivity.33sig <- as.data.frame(dplyr::bind_rows(sig.list))
rownames(connectivity.33sig) <- names(bc.sig.list)

df <- t(connectivity.33sig)




# Gathering data
#rearragne datraframe
data <- reshape::melt(df)
colnames(data) <- c('scNET', 'signature_name', 'connectivity')

#add gene sig length
data$signature_gene_num <- lengths(bc.sig.list)[as.character(data$signature_name)]


#add gene length detected in each celltype net
detected.sig <- list()
for (i in seq_along(table(data$signature_name))){
  signature <- bc.sig.list[[i]]
  sig.f <- lapply(node.list, function(node){length(signature[signature %in% node])})
  names(sig.f) <- names(node.list)
  detected.sig[[i]] <- sig.f
}
names(detected.sig) <- names(table(data$signature_name))
detected.sig <- unlist(detected.sig)
data$detected.sig.num <- detected.sig

#normalize connectivity count values by number of signature genes
data$connectivity.n <- data$connectivity / data$signature_gene_num
data$connectivity.n.detected <- data$connectivity / data$detected.sig.num

#reorder scNET
data$scNET <- factor(data$scNET, levels = c("ECnet", "CAFnet", 'Cancernet','Bnet','Myenet','Tnet'))

#remove signatures with all colSum 0, 6 signatures removed
low.signames <- colnames(as.data.frame(df[,colSums(df) == 0]))
data.f <- data[data$signature_name != low.signames,]


#draw stacked barplot
library(ggplot2)
p1 <- ggplot(data.f, aes(x=signature_name, y=connectivity, fill=scNET))+
  geom_bar(position = 'stack', stat = 'identity') +
  theme_classic()+
  theme(
    panel.grid=element_blank(),
    legend.text=element_text(size=10),
    text = element_text(size=12),
    legend.title = element_blank(),
    axis.title.x = element_blank()
  )+  
  ylab("Absolute # of within group connectivity")+
  xlab("Breast Cancer prognostic signatures") +
  rotate_x_text(45) +
  scale_y_continuous(expand = c(0, 0))


p2 <- ggplot(data.f, aes(x=signature_name, y=connectivity.n, fill=scNET))+
  geom_bar(position = 'stack', stat = 'identity') +
  theme_classic()+
  theme(
    panel.grid=element_blank(),
    legend.text=element_text(size=10),
    text = element_text(size=12),
    legend.title = element_blank(),
    axis.title.x = element_blank()
  )+  
  ylab("Normalized # of within group connectivity")+
  xlab("Breast Cancer prognostic signatures") +
  rotate_x_text(45) +
  scale_y_continuous(expand = c(0, 0))

#remove additional low signature with NAN values (0 / 0)
test <- data.f[!is.nan(data.f$connectivity.n.detected),]
data.f.2 <- data.f[!(data.f$signature_name %in% low.signames),]

p3 <- ggplot(data.f.2, aes(x=signature_name, y=connectivity.n.detected, fill=scNET))+
  geom_bar(position = 'stack', stat = 'identity') +
  theme_classic()+
  theme(
    panel.grid=element_blank(),
    legend.text=element_text(size=10),
    text = element_text(size=12),
    legend.title = element_blank(),
    axis.title.x = element_blank()
  )+  
  ylab("Normalized # of within group connectivity (detected sig)")+
  xlab("Breast Cancer prognostic signatures") +
  rotate_x_text(45) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

#sig gene number and humannet connectivity correlation
humannet <- read.table('~/HumanNetv3/HumanNetv3_networks/HumanNetv3-XC_symbol_LLS.tsv', sep = '\t')

hnv3.conn.list <- list()
for (i in seq_along(bc.sig.list)){
  sig.genes <- bc.sig.list[[i]]
  hnv3.connectivity.sig <- nrow(humannet[(humannet[,1] %in% sig.genes & humannet[,2] %in% sig.genes), ])
  hnv3.conn.list[[i]] <- hnv3.connectivity.sig
}


hnv3.genes <- unique(c(as.character(humannet[,1]), as.character(humannet[,2])))
hnv3.detected.sig <- list()

sig.detect <- lapply(bc.sig.list, function(sig){
  length(sig[sig %in% hnv3.genes])
})
names(sig.detect) <- names(bc.sig.list)
sig.detect <- unlist(sig.detect)


hnv3.connectivity.33sig <- as.data.frame(unlist(hnv3.conn.list))
rownames(hnv3.connectivity.33sig) <- names(bc.sig.list)
hnv3.connectivity.33sig$siggene_num <- lengths(bc.sig.list)
hnv3.connectivity.33sig$siggene_num_detected <- sig.detect
hnv3.connectivity.33sig$name <- rownames(hnv3.connectivity.33sig)
colnames(hnv3.connectivity.33sig)[1] <- 'connectivity'

#draw correlation with labels
library(ggpubr)
library(ggrepel)
options(ggrepel.max.overlaps = 10)

p4 <- ggscatter(hnv3.connectivity.33sig, x = "siggene_num", y = "connectivity", 
                label = "name", repel = TRUE,
                add = 'reg.line', conf.int = T,
                add.params = list(color='blue', fill = 'lightgray')) +
  stat_cor(method = "pearson", label.x = 200, label.y = 1000)

p5 <- ggscatter(hnv3.connectivity.33sig, x = "siggene_num_detected", y = "connectivity", 
                label = "name", repel = TRUE,
                add = 'reg.line', conf.int = T,
                add.params = list(color='blue', fill = 'lightgray')) +
  stat_cor(method = "pearson", label.x = 200, label.y = 1000)



#this is used
pdf('~/HumanNetv3/Graphs/supple3_humannet_connectivity_detected gene number norm.pdf',14,7)
p1 + p5
dev.off()

#not used
pdf('~/HumanNetv3/Graphs/supple3_humannet_connectivity_absolute gene number norm.pdf',14,7)
p2 + p4
dev.off()

#cacluate random probability of BC Tcell network in python


#functional enrichment of Tnet 1 neightbor genes with enrichR
ggi97 <- bc.sig.list[['GGI97']]

ggi97.tnet <- tnet[(tnet[,1] %in% ggi97 | tnet[,2] %in% ggi97),]
ggi97.tnet.genes <- unique(c(as.character(ggi97.tnet[,1]), as.character(ggi97.tnet[,2])))

#write.table(ggi97.tnet.genes,'~/HumanNetv3/table/BC_T_GGI97_neighboring_genes.tsv', quote = F, sep = '\t', row.names = F, col.names = F)

#find LLS strength percentile rank of these genes within Tnet
library(igraph)
tnet$scinet_weight <- NULL
tnet.igraph <- graph_from_data_frame(tnet, directed = F)

tnet.cent <- strength(tnet.igraph, weights = E(tnet.igraph)$LLS)
tnet.cent.rank <- dplyr::percent_rank(tnet.cent)

ggi97.tnet.rank <- tnet.cent.rank[ggi97.tnet.genes] #draw plot of 1order neighbors too
ggi97.tnet.rank <- ggi97.tnet.rank[order(ggi97.tnet.rank, decreasing = T)]


#percentile rank of GGI97 and it's 1 degree neighbor in BC_Tnet
rank.df <- as.data.frame(ggi97.tnet.rank)
rank.df$names <- rownames(rank.df)
rank.df$group <- rep('GGI97 neighboring genes in BC Tnet', nrow(rank.df))

rank.df$genes.to.label <- rank.df$names
rank.df[16:nrow(rank.df),]$genes.to.label <- NA

options(ggrepel.max.overlaps = Inf)
ggplot(rank.df, aes(x = group, y = ggi97.tnet.rank, fill='grey')) +
  geom_boxplot(data = rank.df, alpha = .3) +
  ggrepel::geom_text_repel(aes(label = genes.to.label), color = "black", size = 4, segment.color = "grey") +
  geom_point(size=0) +
  guides(color = "none", fill = "none") +
  theme_bw() +
  labs(
    title = "GGI97 neigbor 427 genes in BC Tnet",
    y = "Percentile Rank",
    x = ''
  )

#save hub genes for survival analysis
saveRDS(rank.df, '~/HumanNetv3/data/GGI97_neighbor_hubgenes_DF.rds')


#returns plot 
#plot1 <- GSAplot(top30.genes,'GO_Biological_Process_2021')
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
  scale_fill_gradient(name = 'overlap percentage',low = '#a1c4fd', high = '#ff9a9e') +
  theme_classic() +
  labs(
    title = title,
    y = '-log10(qvalue)',
    x = "GOBP Term"
  )

return(p)
}

#draw GSA of top n genes in rank df
n <- 30
top.genes <- rownames(head(rank.df, n))
GSAplot(top.genes,'GO_Biological_Process_2021', paste('GGI97 neigbor centrality top',length(top.genes),  'gene GOBP'),20)
GSAplot(top.genes,'KEGG_2021_Human', paste('GGI97 neigbor centrality top',length(top.genes),  'gene KEGG'),20)







#write node attribute table_tnet
tnet.genes <- data.frame(nodes = unique(c(as.character(tnet[,1]),as.character(tnet[,2]))))
cancer.genes <- data.frame(nodes = unique(c(as.character(cancernet[,1]),as.character(cancernet[,2]))))

tnet.genes$GGI97 <- ifelse(tnet.genes$nodes %in% ggi97, 1, 0)
cancer.genes$GGI97 <- ifelse(cancer.genes$nodes %in% ggi97, 1, 0)



write.table(tnet.genes, '~/HumanNetv3/table/HNv3_BC_T_nodetable.tsv', quote = F, col.names = T, row.names = F, sep = '\t')




#draw histogram of distribution and pvalue
tnet.ggi97 <- tnet[(tnet[,1] %in% ggi97 & tnet[,2] %in% ggi97),]


random.24seed <- read.table('~/HumanNetv3/connectivity_analysis/seed_24.cnt', header = F, sep = '\t')

connectivity <- nrow(tnet.ggi97)
ggplot(data=random.24seed, aes(x=V1, y=V2)) +
  scale_x_continuous(trans='log10') +
  geom_bar(stat="identity", fill="steelblue", width = 1, color='black')+
  theme_minimal() +
  ggtitle('Connectivity of 24 random genes in BC Tnet') +
  ylab('Occurence') + xlab('Number of links') +
  geom_vline(aes(xintercept=connectivity), colour="red", linetype="dashed")





#with EnrichR draw GOBP functional enrichment of 427 genes
library(enrichR)
library(dplyr)
library(ggpubr)



p.a <- GSAplot(ggi97.tnet.genes,'KEGG_2021_Human', '429 GGI97 neiboring genes in BC Tnet KEGG', 10)
p.b <- GSAplot(ggi97,'KEGG_2021_Human', 'GGI97 detected 75 genes KEGG Ontology',10)

pdf('~/HumanNetv3/Graphs/KEGG ontology detected GI97genes and neighbors.pdf', 16,8)
p.a + p.b
dev.off()

dbs <- listEnrichrDbs()
dbs <- 'GO_Biological_Process_2021'
enrichr <- enrichr(ggi97.tnet.genes, dbs)

data.gsa <- enrichr$GO_Biological_Process_2021

aaa <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',1))
bbb <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',2))

#add column
data.gsa$overlap_num <- aaa / bbb
data.gsa$log10_qvalue <- -log10(data.gsa$Adjusted.P.value) 

data.gsa.f <- data.gsa[order(data.gsa$log10_qvalue, decreasing = T),]

ggplot(data.gsa.f[1:20,], aes(x = reorder(Term,log10_qvalue),log10_qvalue,  y = log10_qvalue,
                    fill = overlap_num)) + 
  geom_bar(stat = 'identity', width = 0.9, position = position_dodge(width = 0.1)) + 
  coord_flip() + 
  scale_fill_gradient(name = 'overlap percentage',low = '#a1c4fd', high = '#ff9a9e') +
  theme_classic() +
  labs(
    title = "GGI97 neigbor 427 genes in GOBP",
    y = '-log10(qvalue)',
    x = "GOBP Term"
  )



#GSA for 22 genes
dbs <- listEnrichrDbs()
dbs <- 'GO_Biological_Process_2021'

genes <- ggi97[ggi97 %in% ggi97.tnet.genes]
enrichr <- enrichr(genes, dbs)

data.gsa <- enrichr$GO_Biological_Process_2021

aaa <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',1))
bbb <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',2))

#add column
data.gsa$overlap_num <- aaa / bbb
data.gsa$log10_qvalue <- -log10(data.gsa$Adjusted.P.value) 

data.gsa.f <- data.gsa[order(data.gsa$log10_qvalue, decreasing = T),]

ggplot(data.gsa.f[1:20,], aes(x = reorder(Term,log10_qvalue),log10_qvalue,  y = log10_qvalue,
                              fill = overlap_num)) + 
  geom_bar(stat = 'identity', width = 0.9, position = position_dodge(width = 0.1)) + 
  coord_flip() + 
  scale_fill_gradient(name = 'overlap percentage',low = '#a1c4fd', high = '#ff9a9e') +
  theme_classic() +
  labs(
    title = "GGI97 24 genes in GOBP",
    y = '-log10(qvalue)',
    x = "GOBP Term"
  )




