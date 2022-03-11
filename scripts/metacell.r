#!/usr/bin/env Rscript
#need pandas pip install in the envrionment if you don't have it
#usage ./metacell.R [outputfolder] [data.name] [path to count data]
#example ./metacell.R /home3/junhacha/HumanNetv3/Benchmark/MetaCell/Bcell/ BC_metacell_Bcell /home3/junhacha/HumanNetv3/Benchmark/BC_countdata_f_B_cell.tsv
#Junha Cha 

library(ggpubr)
library(metacell)
library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)

#must set up a directory to save Robjects==========================================================
args = commandArgs(trailingOnly=TRUE)

folder.path <- args[1] #'/home3/junhacha/HumanNetv3/Benchmark/MetaCell/Bcell/' #add the slash in the last..its important! the programs makes the last folder if non-existant
data.name <- args[2] #'BC_metacell_Bcell'
data.path <- args[3] #'/home3/junhacha/HumanNetv3/Benchmark/BC_countdata_f_B_cell.tsv'
#==================================================================================================



if(!dir.exists(folder.path)) dir.create(folder.path)
scdb_init(folder.path, force_reinit=T)

mcell_import_scmat_tsv(mat_nm = data.name, 
                       fn = data.path,
                       meta_fn = NULL,
                       dset_nm = data.name)
mat = scdb_mat(data.name)
#===========================================================================================================
print(dim(mat@mat))

#set directory for figures
if(!dir.exists(paste0(folder.path,'/figs/'))) dir.create(paste0(folder.path,"/figs/"))
scfigs_init(paste0(folder.path,"/figs/"))


## Exploring and controlling for quality of raw UMI matrix 
#distribution of UMI plot
mcell_plot_umis_per_cell(data.name)

#list of mito genes for dying cells and IG genes that represent strong clonal signatures, rather than cellular identity
mat = scdb_mat(data.name)
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^IGJ", nms, v=T), 
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T), 
             grep("^IGL", nms, v=T))
bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))

#ignore these genes
mcell_mat_ignore_genes(new_mat_id=data.name, mat_id=data.name, bad_genes, reverse=F) 

#filter cells with low count cells (here less than 800), this value should be modified for each dataset...but here lets just use 800 for all

mcell_mat_ignore_small_cells(data.name, data.name, 500)

## Selecting feature genes based on varaince and 
#calcuate sacaled variance(varaiance divided by mean)
mcell_add_gene_stat(gstat_id=data.name, mat_id=data.name, force=T)

#select informatic features based on threshold
mcell_gset_filter_varmean(gset_id=paste0(data.name, "_feats"), gstat_id=data.name, T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id=paste0(data.name, "_feats"), gstat_id=data.name, T_tot=100, T_top3=2)


#plot the selected genes to refine the parameters
mcell_plot_gstats(gstat_id=data.name, gset_id=paste0(data.name,"_feats"))
## Building the balanced cell graph
#main part of metacell
mcell_add_cgraph_from_mat_bknn(mat_id=data.name, 
                               gset_id = paste0(data.name,"_feats"),
                               graph_id= paste0(data.name,"_graph"),
                               K=100,  #this is the most important parameters
                               dsamp=T)


## Resampling and generating the co-clustering graph
# this process ensures that the metacells are indeed clustered from similar cells
mcell_coclust_from_graph_resamp(
  coc_id=paste0(data.name,"_coc500"), 
  graph_id=paste0(data.name,"_graph"),
  min_mc_size=20, 
  p_resamp=0.75, 
  n_resamp=500) #number of resampling will determine the computation time


mcell_mc_from_coclust_balanced(
  coc_id=paste0(data.name, "_coc500"), 
  mat_id= data.name,
  mc_id= paste0(data.name,"_mc"), 
  K=30, #number of neighbors we wish to minimally associate with each cell (determines the size of metacells) 
  min_mc_size=10, 
  alpha=2) #smaller alpha results in harsh filtering


## Removing outlier cells that deiviate strongly from their metacell expression profile

# detect the outliers with heatmap (don't do this for large datasets)
#mcell_plot_outlier_heatmap(mc_id=paste0(data.name,"_mc"), mat_id = data.name, T_lfc=3)

#filter the detected outlier cells
mcell_mc_split_filt(new_mc_id = paste0(data.name, "_mc_f"), 
            mc_id = paste0(data.name, "_mc"), 
            mat_id=data.name,
            T_lfc=3, #increase when there is too much outliers error 
            plot_mats=F)
            

## Selecting markers and coloring metacells
#test_filtered_mc_f_ is the filtered metacell expression profile
data.metacell <- scdb_mc(paste0(data.name, '_mc_f'))
expr.mc <- data.metacell@mc_fp # this is the final metacell expression matrix

#get PCC correlation and save matrix
data <- expr.mc
data <-as.data.frame(t(data))

#filter data to coding genes
#hnv3 <- read.table('~/HumanNetv3/HumanNetv3_networks/HumanNetv3-XC_symbol_linkOnly.tsv', sep = '\t')
#coverage.genes <- union(as.character(hnv3[,1]), as.character(hnv3[,2]))
#saveRDS(coverage.genes, '/home3/junhacha/HumanNetv3/Benchmark/HNv3_coverage_genes.rds')

hnv3.genes <- readRDS('/home3/junhacha/HumanNetv3/Benchmark/HNv3_coverage_genes.rds')
data.f <- data[,colnames(data) %in% hnv3.genes]

corr.mat <- cor(data, method='pearson')
corr.mat[!lower.tri(corr.mat)] <- NA

corr.net <- reshape2::melt(corr.mat, na.rm = T)

#take only positive values
corr.net <- corr.net[corr.net[,3] > quantile(corr.net[,3], 0.90), ]

#set name for network
output1 <- paste0(folder.path, data.name, '_PCCnet')
cat(paste('network name:', output1,'\n'))

#write network
write.table(corr.net, file = output1, quote = F, row.names = F, sep = '\t', col.names = F)

sort.file.name <- paste0(output1, '_possorted')
system(paste('sort -nrk 3,3', output1, '>', sort.file.name))


#run regression analysis for net1 and net2
cat("Running LLS.py ...\n")
system(paste('python3 /home3/junhacha/bin/LLS.py', sort.file.name))


#remove unsorted network
cat("removing unsorted network\n")
system(paste('rm', output1))


#read LLS benchmark output and save regression plots
LLS.prefiltered <- read.table(file = paste0(sort.file.name, '.binlls'), sep='\t',header = T)


#save plots
cat('draw to pdf...done!\n')
pdf(paste0(folder.path, data.name, '_benchmark.pdf'), width = 14)
par(mfrow = c(1,2))
plot(LLS.prefiltered$GeneCoverage / 18802 *100, LLS.prefiltered$cumLLS, main = 'Benchmark', xlab = 'coverage', ylab = 'cumLLS', type = 'l')
plot(LLS.prefiltered$MeanBinStatistics,LLS.prefiltered$BinLLS, main = 'Regression 1000bin', xlab = 'avg PCC', ylab = "binLLS", pch=19)
dev.off()


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
    labs(title = paste0('Regression_auto_threshold (Metascell Bin = ',bin_suggest,')'), x = "avg PCC", y = "LLS")
  if (!is.na(pcc_thres)) {
    pcc_lab <- as.character(round(pcc_thres, 6))
    y_pos <- min(bm_res[, BinLLS]) + 1
    p1 <- p1 + geom_vline(xintercept = pcc_thres, color = "red", linetype = "dashed")
    p1 <- p1 + annotate(geom = "vline", xintercept = pcc_thres, linetype = "dashed") +
      annotate(geom = "text", label = pcc_lab, x = pcc_thres, y = y_pos, color = "red", angle = 90, vjust = 1.2)
  }
}, error=function(e) {
  p1 <- ggplot() + themo_void() +labs(title = "error")
})


ggsave(paste0(folder.path,data.name,'_Bin_automation.pdf'), plot = p1, width = 7, height = 7)
write.table(bin_suggest,file=paste0(folder.path,data.name,"_bin_tmp"),sep='\t',row.names=F,col.names=F)

network <- read.table(sort.file.name, sep = '\t', header = F)
network.f <- network[1:(1000*bin_suggest),]
write.table(network.f, file = paste0(sort.file.name, '_LLScut.tsv'), quote=F, row.names = F, col.names = F, sep = '\t')





