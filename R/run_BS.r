#!/usr/bin/env Rscript
#this code takes an expression data as an input, runs bigSCale correlation calculaion, and produce a unsorted PCC graph
#rds file must be a seurat file
#./Rscript.R [-i FULL path to exprs UMI data] [-s sort] [-o 20K_monocytes]
#./Rscript.R --help for detail
#RunbigSCale -> coverage filter -> cut correlation -> LLS.py -> draw benchmark graph
#Junha Cha
#07-12-2019

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-s", "--sort", default='sort', help="run absort to get clearer top bins [default %(default)s]")
parser$add_argument("-n", "--nCore", type="integer", default=4, help="number of threads for sorting [default %(default)s]")
parser$add_argument("-c", "--cutoff", type="integer", default=0.95, help="cutoff threshold for compute.network [default %(default)s]")
parser$add_argument("-i", "--input", help="FULL PATH to exprs matrix or seurat rds")
parser$add_argument("-r", "--reuse", default='F', help="type T or F. if -r T,  calcuated bs object will be used [default %(default)s]")
parser$add_argument("-o", "--output", help="filename, current directory will be added to prefix") #no dir prefix
parser$add_argument("-d", "--outdir", help="saved directory")

args <- parser$parse_args()

library(bigSCale)
library(reshape2)
library(float)
library(Seurat)
library(ggplot2)
library(data.table)


ReadCoverage <- function() {
  cat("Reading Coverage\n")
  #read coverage and map data
  coverage <- read.table(file = paste0(pubdata.path, '/coverage_symbol'), colClasses = 'character')[,1]
  #REF <- read.table(file = '/home3/junhacha/public_Data/id_symbol_map', sep = '\t')
  #ref <- as.list(as.character(REF[,1]))
  #names(ref) <- REF[,2] #ref is the map from symbol to id (to compare with gold standard genes)
  #rm(REF)
  return(coverage)
}

ReadData <- function(expr.file) {
  cat("Reading Data\n")
  datatype <- tail(unlist(strsplit(expr.file, "\\.")), n=1)
  if (datatype == 'csv'){
    seperate = ','
    data <- read.table(file = expr.file, header =T, sep = seperate)
  }
  else if (datatype == 'tsv'){
    seperate = '\t'
    data <- read.table(file = expr.file, header =T, sep = seperate)
  }
  else if(datatype == 'rds'){
    seurat <- readRDS(expr.file)
    data <- as.matrix(seurat@assays$RNA@counts)
    rm(seurat)
  }
  else{
    print('unknown file type: check input file')
    quit()
  }
  #read.table and make appropriate adjustment for bigscale input
  return(data)
}


path.to.exprs <- args$input
cutoff <- args$cutoff
sort.type <- args$sort
ncores <- args$nCore
celltype.name <- args$output
out_dir <- args$outdir
setwd(out_dir)
if (args$reuse == 'T'){
  reuse.bs <- TRUE
}else if (args$reuse == 'F'){
  reuse.bs <- FALSE
}

#test arguments
args = commandArgs(trailingOnly=TRUE)
if (sort.type != 'absort' & sort.type != 'sort'){
  cat('Wrong sorting method')
  quit()
} 

if (reuse.bs != F & reuse.bs != T){
  cat('wrong input of bigscale object usage')
  quit()
}


# pubic_data path
pubdata.path <- '/home2/bsb0613/Research/RawData/Network_DB'
#name of the data directory that will be used to name bigscale PCCnet
# directory.name <- tail(unlist(strsplit(path.to.exprs, '/')),  n=3)[1] 
# directory.name <- tail(unlist(strsplit(out_dir,"/")),1)
# output.file <- paste0(directory.name,"_",celltype.name)
output.file <- paste0(celltype.name)

coverage <- ReadCoverage()
data <- ReadData(path.to.exprs)
data.prefiltered <- data[rownames(data) %in% coverage, ]


#run Bigscale and saveRDS bigscale object
#assgin each bigscale object(all and prefiltered) from saved to inprogress file
#dont return anything from mclapply, just save the bigscale objects
if (reuse.bs == F){
  cat("Running BigScale...\n")
  
  cat(paste('start...\n'))
  bigscale.prefiltered <- compute.network(expr.data = data.prefiltered, gene.names = rownames(data.prefiltered), clustering = 'recursive')
  saveRDS(bigscale.prefiltered, file = paste0('./', output.file, '_bs.rds'))
} else if (reuse.bs == T){
  #read from the saved objects
  cat("Calling bigscale objects...\n")
  bigscale.prefiltered <- readRDS(file = paste0('./', output.file, '_bs.rds'))
}



#process corr
zscore <-dbl(bigscale.prefiltered$tot.scores)
colnames(zscore) <- colnames(bigscale.prefiltered$correlations)
#write genespace to know what gene we are making network with (they are filtered)
write.table(colnames(zscore), file = paste0(output.file,'_Genespace_BS.txt'), quote=F, row.names=F, col.names=F)

cat('calculating PCC...')
corr.prefiltered <- cor(zscore, method = 'pearson')


#make PCC unsorted prefiltered using for each 
cat('Melting correlation matrix\n')
corr.prefiltered[!lower.tri(corr.prefiltered)] <- NA
net <- reshape2::melt(corr.prefiltered, na.rm = T)

if (sort.type == 'absort'){
  net <- net[abs(net$value) > quantile(abs(net$value), cutoff), ] #absolute value??
}
if (sort.type == 'sort'){
  net <- net[net$value > quantile(net$value, cutoff), ]
  #net <- net[net$value > 0,]
}



output1 <- paste0('./', output.file, '_BS_PCCnet')
print(paste('network name:', output1))

#writeoutput in inprogress folder
write.table(net, file = output1, quote = F, row.names = F, sep = '\t', col.names = F)


#system absort commaned for net
cat(paste('sorting with...', sort.type, '\n'))
cat("Sorting network...\n")

if (sort.type == 'absort'){
  sort.file.name <- paste0(output1, '_absorted')
  system(paste('~/bin/absort', output1, ncores))
  
}else if(sort.type == 'sort'){
  sort.file.name <- paste0(output1, '_possorted')
  system(paste('sort -nrk 3,3', output1, '>', sort.file.name))
}



#check if the input for regression.py was properly made
if (!file.exists(sort.file.name)){
  cat("Error: Problem with sorting. Read bigscale.rds from file and try again")
  quit()
}


#run regression analysis for net1 and net2
cat("Running LLS.py ...\n")
system(paste('python3 /home3/junhacha/bin/LLS.py', sort.file.name))
# system(paste('python3.6 ~/bin/condLLS.py', sort.file.name, paste0('Genespace_BS_', output.file,'.txt')))

#remove unsorted network
cat("removing unsorted network\n")
system(paste('rm', output1))



#check if the output for regression was properly made
if(!file.exists(paste0(sort.file.name, '.binlls'))){
  cat("Error: Problem with LLS.py analysis prefiltered. Read bigscale.rds from inprogress file and try again\n")
  quit()
}

#read LLS benchmark output and save regression plots
LLS.prefiltered <- read.table(file = paste0(sort.file.name, '.binlls'), sep='\t',header = T)


#save plots
cat('draw to pdf...done!\n')
setwd('./')
pdf(paste0(output.file, '_BS_benchmark.pdf'), width = 14)
par(mfrow = c(1,2))
plot(LLS.prefiltered$GeneCoverage / 18802 *100, LLS.prefiltered$cumLLS, main = 'Benchmark', xlab = 'coverage', ylab = 'cumLLS', type = 'l')
plot(LLS.prefiltered$MeanBinStatistics,LLS.prefiltered$BinLLS, main = paste0('BS Regression 1000bin (',sort.type,')'), xlab = 'avg PCC', ylab = "LLS", pch=19)
dev.off()

folder.path <- getwd()
#get bins and cut to final network
tryCatch({
  bm_res <- fread(file = paste0(sort.file.name, '.binlls'))
  # Get optimal PCC threshold
  prev_lls <- bm_res[1, BinLLS]
  pcc_thres <- NA
  for (row in seq_len(nrow(bm_res))) {
    cur_lls <- bm_res[row, BinLLS]
    if (row > 1) {
      if (abs(prev_lls - cur_lls) >= 2 || cur_lls <= -0.5) {
        last_pos_row <- row - 1
        while (bm_res[last_pos_row, BinLLS] < 0) {
          last_pos_row <- last_pos_row - 1
        }
        pcc_thres <- bm_res[row - 1, MeanBinStatistics]
        bin <- last_pos_row
        cat(paste0("Detected valid PCC threshold: ", pcc_thres, "\n"))
        cat(paste0("Detected valid Bin: ", bin, "\n"))
        break
      }
    }
    if (row > 2) {
      test_df <- bm_res[c(row - 2, row - 1, row)]
      test_df <- test_df[BinLLS < 0]
      if (nrow(test_df) > 1) {
        last_pos_row <- row - 2
        while(bm_res[last_pos_row, BinLLS] < 0) {
          last_pos_row <- last_pos_row - 1
        }
        pcc_thres <- bm_res[last_pos_row, MeanBinStatistics]
        bin <- last_pos_row
        cat(paste0("Detected valid PCC threshold: ", pcc_thres, "\n"))
        cat(paste0("Detected valid Bin: ", bin, "\n"))
        break
      }
    }
    prev_lls <- cur_lls
  }
  
  bin_suggest <- last_pos_row
  
  p1 <- ggplot(data = bm_res, aes(x = MeanBinStatistics, y = BinLLS)) +
    geom_point(show.legend = F) + theme_bw() +
    labs(title = paste0('Regression_auto_threshold (Bin = ',bin_suggest,')'), x = "avg PCC", y = "LLS")
  if (!is.na(pcc_thres)) {
    pcc_lab <- as.character(round(pcc_thres, 6))
    y_pos <- min(bm_res[, BinLLS]) + 1
    p1 <- p1 + geom_vline(xintercept = pcc_thres, color = "red", linetype = "dashed")
    p1 <- p1 + annotate(geom = "vline", xintercept = pcc_thres, linetype = "dashed") +
      annotate(geom = "text", label = pcc_lab, x = pcc_thres, y = y_pos, color = "red", angle = 90, vjust = 1.2)
  }
}, error=function(e) {
  p1 <- ggplot() + theme_void() +labs(title = "error")
})


ggsave(paste0(sort.file.name,'_Bin_automation.pdf'), plot = p1, width = 7, height = 7)
write.table(bin_suggest,file=paste0(sort.file.name,"_bin_tmp"),sep='\t',row.names=F,col.names=F)

network <- read.table(sort.file.name, sep = '\t', header = F)
network.f <- network[1:(1000*bin_suggest),]
write.table(network.f, file = paste0(sort.file.name, '_LLScut.tsv'), quote=F, row.names = F, col.names = F, sep='\t')

