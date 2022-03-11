#this code reads HNv3, and celltype specific net, and make edge and node attribute table
library(qpcR)

hnet <- read.table('~/HumanNetv3/HumanNetv3_networks/HumanNetv3-XC_symbol_LLS.tsv', stringsAsFactors = F)

#sort hnet
#sort edge alphabetically
net <- hnet
net <- as.data.frame(t(apply(net,1,sort)))
net <- net[,c(2,3,1)]
hnet <- net

#read BC nets
bnet <- read.table('~/HumanNetv3/XC_nobatch/BC/net_B_cell_sort.tsv', stringsAsFactors = F)
tnet <- read.table('~/HumanNetv3/XC_nobatch/BC/net_T_cell_sort.tsv', stringsAsFactors = F)
myenet <- read.table('~/HumanNetv3/XC_nobatch/BC/net_Myeloid_sort.tsv',stringsAsFactors = F)
cafnet <- read.table('~/HumanNetv3/XC_nobatch/BC/net_Fibroblast_sort.tsv', stringsAsFactors = F)
ecnet <- read.table('~/HumanNetv3/XC_nobatch/BC/net_EC_sort.tsv', stringsAsFactors = F)


#make node attribute
all.genes <- union(hnet[,1], hnet[,2])
tnet.genes <- union(tnet[,1], tnet[,2])
bnet.genes <- union(bnet[,1], bnet[,2])
myenet.genes <- union(myenet[,1], myenet[,2])
cafnet.genes <- union(cafnet[,1], cafnet[,2])
ecnet.genes <- union(ecnet[,1], ecnet[,2])

#for later use, not really needed here
#tmp <- qpcR:::cbind.na(all.genes, tnet.genes, bnet.genes, myenet.genes, cafnet.genes, ecnet.genes)
#node.table <- as.data.frame(tmp)

#write 1 0 values for coloring in cytoscape
node.table <- data.frame(all.genes = all.genes)
node.table$tnet <- ifelse(all.genes %in% tnet.genes, 1, 0)
node.table$bnet <- ifelse(all.genes %in% bnet.genes, 1, 0)
node.table$myenet <- ifelse(all.genes %in% myenet.genes, 1, 0)
node.table$cafnet <- ifelse(all.genes %in% cafnet.genes, 1, 0)
node.table$ecnet <- ifelse(all.genes %in% ecnet.genes, 1, 0)

#save output
write.table(node.table, '~/HumanNetv3/XC_nobatch/BC/HNv3_BC_nodetable.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

#make edge attribute
hnet_links = paste(hnet[,1], hnet[,2], sep = '|')
tnet_links = paste(tnet[,1], tnet[,2], sep = '|')
bnet_links = paste(bnet[,1], bnet[,2], sep = '|')
myenet_links = paste(myenet[,1], myenet[,2], sep = '|')
cafnet_links = paste(cafnet[,1], cafnet[,2], sep = '|')
ecnet_links = paste(ecnet[,1], ecnet[,2], sep = '|')


hnet$tnet_links <- ifelse(hnet_links %in% tnet_links, 1, 0)
hnet$bnet_links <- ifelse(hnet_links %in% bnet_links, 1, 0)
hnet$myenet_links <- ifelse(hnet_links %in% myenet_links, 1, 0)
hnet$cafnet_links <- ifelse(hnet_links %in% cafnet_links, 1, 0)
hnet$ecnet_links <- ifelse(hnet_links %in% ecnet_links, 1, 0)


write.table(hnet, '~/HumanNetv3/XC_nobatch/BC/HNv3_BC_edgetable.tsv', quote = F, col.names = T, row.names = F, sep = '\t')


