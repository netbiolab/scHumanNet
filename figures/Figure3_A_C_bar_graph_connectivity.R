#get bar graph of within connectivity of ICm geneset for each celltype network
library(ggpubr)


#read ICI gene list
ici.l <- readRDS('~/public_Data/ICI/immune_checkpoints_ligand.rds')
ici.r <- readRDS('~/public_Data/ICI/immune_checkpoints_receptor.rds')
ici.all <- unique(c(ici.r, ici.l)) #40


#read celltype specific networks for each cancer 
setwd('~/HumanNetv3/XC_nobatch/')
net.list.crc <- readRDS('./CRC/CRC_celltypeNet_list_LLS.rds')
net.list.bc <- readRDS('./BC/BC_celltypeNet_list_LLS.rds')
net.list.ovc <- readRDS('./OvC/OvC_celltypeNet_list_LLS.rds')
net.list.lc <- readRDS('./LC/LC_celltypeNet_list_LLS.rds')

#



#get connectivity of Tcell Net for each cancer type
ici.within.1 <- unlist(lapply(net.list.crc, function(net) nrow(net[(net[,1] %in% ici.all & net[,2] %in% ici.all), ])))
ici.within.2 <- unlist(lapply(net.list.bc, function(net) nrow(net[(net[,1] %in% ici.all & net[,2] %in% ici.all), ])))
ici.within.3 <- unlist(lapply(net.list.lc, function(net) nrow(net[(net[,1] %in% ici.all & net[,2] %in% ici.all), ])))
ici.within.4 <- unlist(lapply(net.list.ovc, function(net) nrow(net[(net[,1] %in% ici.all & net[,2] %in% ici.all), ])))

#filter to 5 common celltype
ici.within.1 <- ici.within.1[c('Cancer','Fibroblast','EC','Myeloid','B_cell','T_cell')]
ici.within.2 <- ici.within.2[c('Cancer','Fibroblast','EC','Myeloid','B_cell','T_cell')]
ici.within.3 <- ici.within.3[c('Cancer','Fibroblast','EC','Myeloid','B_cell','T_cell')]
ici.within.4 <- ici.within.4[c('Cancer','Fibroblast','EC','Myeloid','B_cell','T_cell')]

data1 <- data.frame(ici.within.1)
data1$cell <- rownames(data1)
p1 <- ggbarplot(data1, x = 'cell', y = 'ici.within.1', ylab ='Connectivity',xlab = '', title = 'Icm connectivity in CRC') + 
  rotate_x_text(-45) + 
  scale_y_continuous(expand = c(0, 0)) 

data2 <- data.frame(ici.within.2)
data2$cell <- rownames(data2)
p2 <- ggbarplot(data2, x = 'cell', y = 'ici.within.2', ylab ='Connectivity', xlab = '',title = 'Icm connectivity in BC') + 
  rotate_x_text(-45) +
  scale_y_continuous(expand = c(0, 0))

data3 <- data.frame(ici.within.3)
data3$cell <- rownames(data3)
p3 <- ggbarplot(data3, x = 'cell', y = 'ici.within.3', ylab ='Connectivity', xlab = '',title = 'Icm connectivity in LC') + 
  rotate_x_text(-45) +
  scale_y_continuous(expand = c(0, 0))

data4 <- data.frame(ici.within.4)
data4$cell <- rownames(data4)
p4 <- ggbarplot(data4, x = 'cell', y = 'ici.within.4', ylab ='Connectivity', xlab = '',title = 'Icm connectivity in OvC') + 
  rotate_x_text(-45) +
  scale_y_continuous(expand = c(0, 0))

library(patchwork)
pdf('~/HumanNetv3/Graphs/barplot_ICm_connectivity_filtered_celltypes.pdf', 6,6)
(p1+p2) / (p3 + p4)
dev.off()


#compare it with DEGs
setwd('~/HumanNetv3/DEGlists/')

## read data using loop
files.crc <- list.files(path = '~/HumanNetv3/DEGlists/', pattern="^CRC")
deg.crc <- NULL
for (f in files.crc) {
  dat.crc <- read.table(f, header=T, sep="\t", na.strings="", colClasses="character")
  deg.crc <- rbind(deg.crc, dat.crc)
}

files.bc <- list.files(path = '~/HumanNetv3/DEGlists/', pattern="^BC")
deg.bc <- NULL
for (f in files.bc) {
  dat.bc <- read.table(f, header=T, sep="\t", na.strings="", colClasses="character")
  deg.bc <- rbind(deg.bc, dat.bc)
}

files.lc <- list.files(path = '~/HumanNetv3/DEGlists/', pattern="^LC")
deg.lc <- NULL
for (f in files.lc) {
  dat.lc <- read.table(f, header=T, sep="\t", na.strings="", colClasses="character")
  deg.lc <- rbind(deg.lc, dat.lc)
}

files.ovc <- list.files(path = '~/HumanNetv3/DEGlists/', pattern="^OvC")
deg.ovc <- NULL
for (f in files.ovc) {
  dat.ovc <- read.table(f, header=T, sep="\t", na.strings="", colClasses="character")
  deg.ovc <- rbind(deg.ovc, dat.ovc)
}

#make cluster column as levels
deg.bc[, 'cluster'] <- as.factor(deg.bc[, 'cluster'])
deg.lc[, 'cluster'] <- as.factor(deg.lc[, 'cluster'])
deg.ovc[, 'cluster'] <- as.factor(deg.ovc[, 'cluster'])
deg.crc[, 'cluster'] <- as.factor(deg.crc[, 'cluster'])


#get upregulated DEGs 0.25 as avg_logFC threshold
deg.bc <- deg.bc[deg.bc$avg_logFC > 0.25, ]
deg.crc <- deg.crc[deg.crc$avg_logFC > 0.25, ]
deg.lc <- deg.lc[deg.lc$avg_logFC > 0.25, ]
deg.ovc <- deg.ovc[deg.ovc$avg_logFC > 0.25, ]


#get Icm DEGs of all celltypes Net for each cancer type
deg.crc.ici <- deg.crc[deg.crc$gene %in% ici.all,]
deg.bc.ici <- deg.bc[deg.bc$gene %in% ici.all,]
deg.lc.ici <- deg.lc[deg.lc$gene %in% ici.all,]
deg.ovc.ici <- deg.ovc[deg.ovc$gene %in% ici.all,]


data5 <- as.data.frame(table(deg.crc.ici$cluster))
p5 <- ggbarplot(data5, x = 'Var1', y = 'Freq', ylab ='Num of upregulated Icm',xlab = '', title = 'Colorectal Cancer') + 
  rotate_x_text(-45) + 
  theme(text = element_text(size = 16))  + 
  scale_y_continuous(expand = c(0, 0)) 


data6 <- as.data.frame(table(deg.bc.ici$cluster))
p6 <- ggbarplot(data6, x = 'Var1', y = 'Freq', ylab ='Num of upregulated Icm',xlab = '', title = 'Breast Cancer') + 
  rotate_x_text(-45) + 
  theme(text = element_text(size = 16))  + 
  scale_y_continuous(expand = c(0, 0)) 

data7 <- as.data.frame(table(deg.lc.ici$cluster))
p7 <- ggbarplot(data7, x = 'Var1', y = 'Freq', ylab ='Num of upregulated Icm',xlab = '', title = 'Lung Cancer') + 
  rotate_x_text(-45) + 
  theme(text = element_text(size = 16))  + 
  scale_y_continuous(expand = c(0, 0)) 

data8 <- as.data.frame(table(deg.ovc.ici$cluster))
p8 <- ggbarplot(data8, x = 'Var1', y = 'Freq', ylab ='Num of upregulated Icm',xlab = '', title = 'Ovarian Cancer') + 
  rotate_x_text(-45) + 
  theme(text = element_text(size = 16))  + 
  scale_y_continuous(expand = c(0, 0)) 

library(patchwork)
pdf('~/HumanNetv3/Graphs/barplot_ICm_upDEGs_filtered_celltypes.pdf', 6,6)
(p5+p6) / (p7 + p8)
dev.off()




#draw grouped bar plot for how many ICMs are detected for DEG or scNET
#these are for scNET
ici.count.crc <- unlist(lapply(net.list.crc, function(net){ 
  all.genes <- unique(c(net[,1],net[,2]))
  print(all.genes[all.genes %in% ici.all])
  return(length(all.genes[all.genes %in% ici.all]))
}))

ici.count.bc <- unlist(lapply(net.list.bc, function(net){ 
  all.genes <- unique(c(net[,1],net[,2]))
  print(all.genes[all.genes %in% ici.all])
  return(length(all.genes[all.genes %in% ici.all]))
}))

ici.count.lc <- unlist(lapply(net.list.lc, function(net){ 
  all.genes <- unique(c(net[,1],net[,2]))
  print(all.genes[all.genes %in% ici.all])
  return(length(all.genes[all.genes %in% ici.all]))
}))

ici.count.ovc <- unlist(lapply(net.list.ovc, function(net){ 
  all.genes <- unique(c(net[,1],net[,2]))
  print(all.genes[all.genes %in% ici.all])
  return(length(all.genes[all.genes %in% ici.all]))
}))

#filter to 5 common celltype
desired.order <- c('Fibroblast','EC','B_cell','Myeloid','T_cell')
ici.net.crc <- ici.count.crc[desired.order]
ici.net.bc <- ici.count.bc[desired.order]
ici.net.lc <- ici.count.lc[desired.order]
ici.net.ovc <- ici.count.ovc[desired.order]

#merge the DEG and LLS data
library(dplyr)
data5 <- data5 %>% slice(match(desired.order, Var1))#crc
data6 <- data6 %>% slice(match(desired.order, Var1))#BC
data7 <- data7 %>% slice(match(desired.order, Var1))#LC
data8 <- data8 %>% slice(match(desired.order, Var1))#OVC

#add to data df
column.name <- c('celltype', 'DEG_count', 'NET_count')
data5$NET <- ici.net.crc
colnames(data5) <- column.name
data6$NET <- ici.net.bc
colnames(data6) <- column.name
data7$NET <- ici.net.lc
colnames(data7) <- column.name
data8$NET <- ici.net.ovc
colnames(data8) <- column.name

#melt it to plot grouped bar chart
data5 <- reshape2::melt(data5[,c('celltype','DEG_count','NET_count')],id.vars = 1)
data6 <- reshape2::melt(data6[,c('celltype','DEG_count','NET_count')],id.vars = 1)
data7 <- reshape2::melt(data7[,c('celltype','DEG_count','NET_count')],id.vars = 1)
data8 <- reshape2::melt(data8[,c('celltype','DEG_count','NET_count')],id.vars = 1)


#draw plot
p9 <- ggplot(data5,aes(x = celltype,y = value)) + 
  geom_bar(aes(fill = variable, alpha=0.5),stat = "identity",position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c('yellow','blue'))+
  #theme(text = element_text(size = 16))  + 
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle('Colorectal Cancer') +
  xlab('') +
  ylab('Num of Icm genes detected') +
  rotate_x_text(-45)



p10 <- ggplot(data6,aes(x = celltype,y = value)) + 
  geom_bar(aes(fill = variable, alpha=0.5),stat = "identity",position = "dodge") +
  scale_fill_manual(values = c('yellow','blue')) +
  theme_classic() +
  #theme(text = element_text(size = 16))  + 
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle('Breast Cancer') +
  xlab('') +
  ylab('Num of Icm genes detected') +
  rotate_x_text(-45)


p11 <- ggplot(data7,aes(x = celltype,y = value)) + 
  geom_bar(aes(fill = variable, alpha=0.5),stat = "identity",position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c('yellow','blue')) +
  #theme(text = element_text(size = 16))  + 
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle('Lung Cancer') +
  xlab('') +
  ylab('Num of Icm genes detected') +
  rotate_x_text(-45)



p12 <- ggplot(data8,aes(x = celltype,y = value)) + 
  geom_bar(aes(fill = variable, alpha=0.5),stat = "identity",position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c('yellow','blue')) +
  #theme(text = element_text(size = 16))  + 
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle('Ovarian Cancer') +
  xlab('') +
  ylab('Num of Icm genes detected') +
  rotate_x_text(-45)


library(patchwork)
pdf('~/HumanNetv3/Graphs/barplot_ICm_DEGvsNET_filtered_celltypes.pdf', 10,10)
(p9+p10) / (p11 + p12)
dev.off()




#make a table for detected ICMs for each cancer each celltype
setwd('~/HumanNetv3/table/')
col.names <- c('Fibroblast_scNET','Fibroblast_DEG','EC_scNET','EC_DEG','B_cell_scNET','B_cell_DEG',
               'Myeloid_scNET','Myeloid_DEG','T_cell_scNET','T_cell_DEG')

#CRC
ici.genes.crc <- lapply(net.list.crc, function(net){ 
  all.genes <- unique(c(net[,1],net[,2]))
  print(all.genes[all.genes %in% ici.all])
  return((all.genes[all.genes %in% ici.all]))
})


#combine with DEG detected ICMs
deg.bcell <- deg.crc.ici$gene[deg.crc.ici$cluster == 'B_cell']
deg.tcell <- deg.crc.ici$gene[deg.crc.ici$cluster == 'T_cell']
deg.caf <- deg.crc.ici$gene[deg.crc.ici$cluster == 'Fibroblast']
deg.ec <- deg.crc.ici$gene[deg.crc.ici$cluster == 'EC']
deg.mye <- deg.crc.ici$gene[deg.crc.ici$cluster == 'Myeloid']

#combine them all...this is a messy code...
crc.df <-  qpcR:::cbind.na(ici.genes.crc[['Fibroblast']], deg.caf,
                           ici.genes.crc[['EC']], deg.ec,
                           ici.genes.crc[['B_cell']], deg.bcell,
                           ici.genes.crc[['Myeloid']], deg.mye,
                           ici.genes.crc[['T_cell']], deg.tcell)

colnames(crc.df) <- col.names
write.table(crc.df, sep='\t',quote=F, row.names = F, col.names = T, file='./CRC_detectedICM_DEGvssNET.tsv')




#BC
ici.genes.bc <- lapply(net.list.bc, function(net){ 
  all.genes <- unique(c(net[,1],net[,2]))
  print(all.genes[all.genes %in% ici.all])
  return((all.genes[all.genes %in% ici.all]))
})


#combine with DEG detected ICMs
bc.deg.bcell <- deg.bc.ici$gene[deg.bc.ici$cluster == 'B_cell']
bc.deg.tcell <- deg.bc.ici$gene[deg.bc.ici$cluster == 'T_cell']
bc.deg.caf <- deg.bc.ici$gene[deg.bc.ici$cluster == 'Fibroblast']
bc.deg.ec <- deg.bc.ici$gene[deg.bc.ici$cluster == 'EC']
bc.deg.mye <- deg.bc.ici$gene[deg.bc.ici$cluster == 'Myeloid']

#combine them all...this is a messy code...
bc.df <-  qpcR:::cbind.na(ici.genes.bc[['Fibroblast']], bc.deg.caf,
                           ici.genes.bc[['EC']], bc.deg.ec,
                           ici.genes.bc[['B_cell']], bc.deg.bcell,
                           ici.genes.bc[['Myeloid']], bc.deg.mye,
                           ici.genes.bc[['T_cell']], bc.deg.tcell)

colnames(bc.df) <- col.names
write.table(bc.df, sep='\t',quote=F, row.names = F, col.names = T, file='./BC_detectedICM_DEGvssNET.tsv')


#LC
ici.genes.lc <- lapply(net.list.lc, function(net){ 
  all.genes <- unique(c(net[,1],net[,2]))
  print(all.genes[all.genes %in% ici.all])
  return((all.genes[all.genes %in% ici.all]))
})


#combine with DEG detected ICMs
lc.deg.bcell <- deg.lc.ici$gene[deg.lc.ici$cluster == 'B_cell']
lc.deg.tcell <- deg.lc.ici$gene[deg.lc.ici$cluster == 'T_cell']
lc.deg.caf <- deg.lc.ici$gene[deg.lc.ici$cluster == 'Fibroblast']
lc.deg.ec <- deg.lc.ici$gene[deg.lc.ici$cluster == 'EC']
lc.deg.mye <- deg.lc.ici$gene[deg.lc.ici$cluster == 'Myeloid']

#combine them all...this is a messy code...
lc.df <-  qpcR:::cbind.na(ici.genes.lc[['Fibroblast']], lc.deg.caf,
                          ici.genes.lc[['EC']], lc.deg.ec,
                          ici.genes.lc[['B_cell']], lc.deg.bcell,
                          ici.genes.lc[['Myeloid']], lc.deg.mye,
                          ici.genes.lc[['T_cell']], lc.deg.tcell)

colnames(lc.df) <- col.names
write.table(lc.df, sep='\t',quote=F, row.names = F, col.names = T, file='./LC_detectedICM_DEGvssNET.tsv')




#OVC
ici.genes.ovc <- lapply(net.list.ovc, function(net){ 
  all.genes <- unique(c(net[,1],net[,2]))
  print(all.genes[all.genes %in% ici.all])
  return((all.genes[all.genes %in% ici.all]))
})


#combine with DEG detected ICMs
ovc.deg.bcell <- deg.ovc.ici$gene[deg.ovc.ici$cluster == 'B_cell']
ovc.deg.tcell <- deg.ovc.ici$gene[deg.ovc.ici$cluster == 'T_cell']
ovc.deg.caf <- deg.ovc.ici$gene[deg.ovc.ici$cluster == 'Fibroblast']
ovc.deg.ec <- deg.ovc.ici$gene[deg.ovc.ici$cluster == 'EC']
ovc.deg.mye <- deg.ovc.ici$gene[deg.ovc.ici$cluster == 'Myeloid']

#combine them all...this is a messy code...
ovc.df <-  qpcR:::cbind.na(ici.genes.ovc[['Fibroblast']], ovc.deg.caf,
                          ici.genes.ovc[['EC']], ovc.deg.ec,
                          ici.genes.ovc[['B_cell']], ovc.deg.bcell,
                          ici.genes.ovc[['Myeloid']], ovc.deg.mye,
                          ici.genes.ovc[['T_cell']], ovc.deg.tcell)

colnames(ovc.df) <- col.names
write.table(ovc.df, sep='\t',quote=F, row.names = F, col.names = T, file='./OVC_detectedICM_DEGvssNET.tsv')

