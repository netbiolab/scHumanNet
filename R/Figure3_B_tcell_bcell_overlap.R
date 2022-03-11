#get B anc T cell specific gene set from GO term

#readdata
data.tcell <- read.table('~/HumanNetv3/gene2go_tcell.tsv', header = F, sep = '\t')
data.bcell <- read.table('~/HumanNetv3/gene2go_bcell.tsv', header = F, sep = '\t')

#filter data by human first
data.tcell <- data.tcell[data.tcell$V1 == '9606', ]
data.bcell <- data.bcell[data.bcell$V1 == '9606', ]


table(data.tcell$V4)
table(data.bcell$V4)

data.tcell <- data.tcell[data.tcell$V4 == 'TAS' | 
                         data.tcell$V4 == 'IDA' |
                         data.tcell$V4 == 'IMP' |
                         data.tcell$V4 == 'IGI',]

data.bcell <- data.bcell[data.bcell$V4 == 'TAS' | 
                         data.bcell$V4 == 'IDA' |
                         data.bcell$V4 == 'IMP' |
                         data.bcell$V4 == 'IGI',]


tcell.genes <- unique(data.tcell$V2)
bcell.genes <- unique(data.bcell$V2)

length(tcell.genes)
length(bcell.genes)

#map gene ID to symbol
ncbi.ref <- read.table('~/public_Data/id_symbol_map', quote = '', sep = '\t')  #map from ncbi

ncbi.ref.mapped.t <- ncbi.ref[ncbi.ref$V1 %in% tcell.genes,]
tcell.genes.sym <- as.character(ncbi.ref.mapped.t[,2])
ncbi.ref.mapped.b <- ncbi.ref[ncbi.ref$V1 %in% bcell.genes,]
bcell.genes.sym <- as.character(ncbi.ref.mapped.b[,2])

#every gene is mapped
length(tcell.genes)
length(tcell.genes.sym)
length(bcell.genes)
length(bcell.genes.sym)


#save genelist
saveRDS(tcell.genes.sym, "~/HumanNetv3/data/tcell_GOgeneset.rds")
saveRDS(bcell.genes.sym, "~/HumanNetv3/data/bcell_GOgeneset.rds")



#get T/Bcell LLS ranked nodes and DEG
hnv3 <- read.table('~/HumanNetv3/HumanNetv3_networks/HumanNetv3-XC_symbol_linkOnly.tsv', sep = '\t')

setwd('~/HumanNetv3/table/')
bc.rank <- read.table('./BC_LLS_strength_rank.tsv', sep = '\t', header = T, stringsAsFactors = F)
lc.rank <- read.table('./LC_LLS_strength_rank.tsv', sep = '\t', header = T, stringsAsFactors = F)
ovc.rank <- read.table('./OvC_LLS_strength_rank.tsv', sep = '\t', header = T, stringsAsFactors = F)
crc.rank <- read.table('./CRC_LLS_strength_rank.tsv', sep = '\t', header = T, stringsAsFactors = F)

bc.deg <- read.table('~/HumanNetv3/data/DEGlist_BC_byCelltype_all_logcounts.tsv', sep = '\t', header = T)
crc.deg <- read.table('~/HumanNetv3/data/DEGlist_CRC_byCelltype_all_logcounts.tsv', sep = '\t', header = T)
lc.deg <- read.table('~/HumanNetv3/data/DEGlist_LC_byCelltype_all_logcounts.tsv', sep = '\t', header = T)
ovc.deg <- read.table('~/HumanNetv3/data/DEGlist_OvC_byCelltype_all_logcounts.tsv', sep = '\t', header = T)
#*********************************************************************************

###############################################################
celltype <- "B_cell" # Tcell / B_cell
rank.celltype <- lc.rank # bc/lc/ovc/crc
deg.celltype <- lc.deg #bc/lc/ovc/crc
filename <- 'LC_Bcell_CodingGeneFiltered'
###############################################################

#retrieve DEG gene list
deg.celltype <- deg.celltype[deg.celltype$cluster == celltype,] #filter by celltype
#sort by pvalue or logfold
deg.celltype <- deg.celltype[deg.celltype$avg_logFC > 0, ] #take only positive DEGs

#filter deg by coverage genes
coverage.genes <- union(as.character(hnv3[,1]), as.character(hnv3[,2]))
deg.celltype <- deg.celltype[deg.celltype$gene %in% coverage.genes, ]

deg.celltype$per.rank.deg.logfc <- dplyr::percent_rank(deg.celltype$avg_logFC) #rank the DEG genes by lof fold (qvalues were used to filter)
#deg.tcell$per.rank.deg.pval <- 1 - dplyr::percent_rank(deg.tcell$p_val_adj) #reverse percent rank



#retreive LLS genes by same amount of DEG
top.n <- nrow(deg.celltype)
rank.celltype <- rank.celltype[,c("rowname",celltype)]
rank.celltype <- rank.celltype[order(rank.celltype[,celltype], decreasing = T),]
top.genes <- rank.celltype$rowname[1:top.n]




#draw venndigaram
#see overlap
#library(VennDiagram)
#devtools::install_github("yanlinlin82/ggvenn")
set1 <- top.genes
set2 <- as.character(deg.celltype$gene)
inter <- intersect(set1,set2)

library(ggvenn)
setlist <- list(LLS_top = set1, DEG_pos=set2)
p1 <- ggvenn(setlist, text_size = 2, stroke_size = 0.4)


#legends = c(paste0("LLS_top",top.n), "DEG_pos")
################################################
#set1 <- top.genes
#set2 <- as.character(deg.celltype$gene)
#inter <- intersect(set1,set2)
#
#
#p1 <-draw.pairwise.venn(
#  area1 = length(set1),
#  area2 = length(set2),
#  cross.area = length(inter),
#  
#  #colors
#  fill = c('#E5B3BB', '#7B92AA'),
#  #lty = "blank",
#  col = 'black',
#  
#  #legends
#  category = legends,
#  cat.pos = c(0,0),
#  
#  ind=FALSE
#)
#grid.newpage()
#grid.draw(p1)
#grid.newpage()


#library(Vennerable)
#test <- list(set1,set2)
#names(test) <- c('LLS_top', 'DEG_pos')
#Vdeg <- Venn(test)
#plot(Vdeg, doWeights = TRUE)



#get geneA geneB geneC
geneA <- setdiff(set1, inter)
geneB <- inter
geneC <- setdiff(set2, inter)

#draw bar graph
library(ggplot2)

if (celltype == "T_cell"){
  print('tcell_calculated!')
  a <- length(geneA[geneA %in% tcell.genes.sym]) / length(geneA) * 100
  b <- length(geneB[geneB %in% tcell.genes.sym]) / length(geneB) * 100
  c <- length(geneC[geneC %in% tcell.genes.sym]) / length(geneC) * 100
}
if (celltype =="B_cell"){
  print('bcell_calculated!')
  a <- length(geneA[geneA %in% bcell.genes.sym]) / length(geneA) * 100
  b <- length(geneB[geneB %in% bcell.genes.sym]) / length(geneB) * 100
  c <- length(geneC[geneC %in% bcell.genes.sym]) / length(geneC) * 100
}



df <- data.frame(Gene_set = c('setA', 'setB','setC'),
                 Overlap_percentage = c(a, b, c))

p2 <- ggplot(data=df, aes(x=Gene_set, y=Overlap_percentage)) +
  geom_bar(stat="identity", width=0.5) + theme_minimal()

library(patchwork)
library(ggplotify)
pdf(paste0('~/HumanNetv3/Graphs/GO_overlap_',filename,'.pdf'), 3,6)
p1 / p2
dev.off()







