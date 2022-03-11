#this code calcuates raw PCC coexpressiox network and write to file 
library(Seurat)
library(parallel)
library(doParallel)

setwd('/home3/junhacha/HumanNetv3/Benchmark/rawPCC/')
data <- read.table('/home3/junhacha/HumanNetv3/Benchmark/BC_countdata_EC.tsv', sep = '\t', header = T)
cell.type <- 'EC'

seurat <- CreateAssayObject(counts=data)
#normalize data
seurat <- NormalizeData(seurat)
seurat.data = GetAssayData(seurat)

#filter by coverage
coverage <- readRDS('/home3/junhacha/HumanNetv3/data/coverage_genes_all.rds')
seurat.data <- seurat.data[rownames(seurat.data) %in% coverage, ]


cores <- makeCluster(detectCores(), type='PSOCK') # grabs max available
cores <- 6  
options('mc.cores' = cores)
registerDoParallel(cores)

mat <- t(as.matrix(seurat.data))
res <- foreach(i = seq_len(ncol(mat)),
               .combine = rbind,
               .multicombine = TRUE,
               .inorder = FALSE,
               .packages = c('data.table', 'doParallel')) %dopar% {
                 cor(mat[,i], mat, method = 'pearson')
               }

#because most is zero...we will get a lot of NAs...
rownames(res) <- colnames(res)


res[!lower.tri(res)] <- NA
corr.net <- reshape2::melt(res, na.rm = T)

#take links above 0.8 for raw PCC..we will not perfrom LLS thresholding
corr.net <- corr.net[corr.net[,3] > 0.8,] 

#set name for network
output <- paste0(cell.type, '_rawPCCnet0.8')

#write network
write.table(corr.net, file = paste0('/home3/junhacha/HumanNetv3/Benchmark/rawPCC/',output), quote = F, row.names = F, sep = '\t', col.names = F)


#sort network
sort.file.name <- paste0(output, '_possorted')
system(paste('sort -nrk 3,3', output, '>', sort.file.name))


#remove unsorted network
cat("removing unsorted network\n")
system(paste('rm', output))


