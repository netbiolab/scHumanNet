#get survival data for Breast Cancer TCGA do this in R4 env
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(ggfortify)
library(survival)
library(EDASeq)
library(survminer)
library(gridExtra)
library(TCGAbiolinks) #download this using devtools::install_github()

#set working directory
setwd('~/HumanNetv3/TCGA/')

query <- GDCquery(project="TCGA-BRCA", 
                  data.type = "Gene Expression Quantification", 
                  data.category = "Transcriptome Profiling", 
                  workflow.type = "HTSeq - Counts", 
                  experimental.strategy = "RNA-Seq")


samplesDown <- getResults(query,cols=c("cases"))
#get primary solid tumors(TP) and TN (e.g.TN are normal solid tissue)
#for normalization purposes with DESeq2
#500 TP
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")

#theres 44 normal
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")


queryDown <- GDCquery(project = "TCGA-BRCA", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP, dataSmNT))

GDCdownload(queryDown)
se.data <- GDCprepare(queryDown) #the data is SummerizedExperiment class
saveRDS(se.data, './BRCA_TCGA_TP_NT_samples_SummarizedExperiment.rds')
#se.data <- GDCprepare(query = queryDown,directory = '~/GDCdata/')

#no samples are dropped in the cor.cut threshold 0.6, data are harmonious
data.preprocessed <- TCGAanalyze_Preprocessing(object = se.data, 
                                               cor.cut = 0.6,
                                               datatype = "HTSeq - Counts")

dataNorm <- TCGAanalyze_Normalization(tabDF = data.preprocessed,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent") 



#map ENSG to GENE symbol
library(EnsDb.Hsapiens.v86) #coordinates for GRch38
edb <- EnsDb.Hsapiens.v86
#read data
data.brca <- as.data.frame(dataNorm)

#map ensg to genename
gene.ids <- mapIds(edb, keys = rownames(data.brca), column = 'GENENAME', keytype = 'GENEID')

#make a column in the data
data.brca$geneid <- gene.ids

#filter data to mapped genes
data.brca <- data.brca[!is.na(data.brca$geneid),]

#remove duplicates
data.brca <- data.brca[!duplicated(data.brca$geneid),]

#make gene ids rownames
rownames(data.brca) <- data.brca$geneid
data.brca$geneid <- NULL

saveRDS(data.brca, './TCGA_BRCA_TP_NT_samples_GCNormCounts.rds')


#get clinical data and save
#data from counts(coldata) and clin.data can be mapped via unique IDs
clin.brca <- GDCquery_clinic("TCGA-BRCA", "clinical")

#filter samples to those in our RNA data
sample.id <- se.data@colData$patient

#there are multipe RNA samples with same patients..that's why it doesn't match
#we will normalize with all samples TP&NT included then filter for TP 
dim(clin.brca[clin.brca$submitter_id %in% sample.id,])
dim(clin.brca)

#normalize count data with DESeq2 pipeline
library(DESeq2)

metadata <- se.data@colData
dds <- DESeqDataSetFromMatrix(countData = data.brca, colData = metadata, design = ~ sample_type) #NP and TP of design does not affect normalization...
dds <- estimateSizeFactors(dds)
data.brca.n <- counts(dds, normalized=TRUE)


#filter se.data by clinical info, and with that to data and save
#normal Solid Tissue came from cancer patients as well so there sample patients with both NT and TP samples, # use TP samples!
se.data.f <- se.data[,se.data@colData$sample_type == 'Primary Tumor']
data.brca.n.f <- data.brca.n[,colnames(assay(se.data.f))] #filter normalized data to TP samples

head(se.data.f@colData$patient)
clin.brca.f <- clin.brca[clin.brca$submitter_id %in% se.data.f@colData$patient, ] #filter metadata to 500 samples

dim(clin.brca.f)
dim(data.brca.n.f)
table(se.data.f@colData$patient %in% clin.brca.f$submitter_id) #check that every RNA seq data column is represented in the metadata

lnames.rna <- unlist(lapply(strsplit(colnames(data.brca.n.f), '-'), function(x){paste(x[1], x[2], x[3], sep='-')}))

#order colnames so it matches with clin.brca.f
data.brca.n.f <- data.brca.n.f[, clin.brca.f$submitter_id] 

#filter out 12 male samples and NA samples...they might be unwanted sources of variation
clin.brca.f <- clin.brca.f[clin.brca.f$gender == 'female',]
clin.brca.f <- clin.brca.f[!is.na(clin.brca.f$gender), ]
data.brca.n.f <- data.brca.n.f[,colnames(data.brca.n.f) %in% clin.brca.f$submitter_id]

#make sure they match
identical(colnames(data.brca.n.f), clin.brca.f$submitter_id)

data.brca.n.f <- as.data.frame(data.brca.n.f)
#save final se, metadata, and normalized rna data
saveRDS(data.brca.n.f, './BRCA_TP_DESeq2_norm_counts_Female_samples.rds')
saveRDS(clin.brca.f, './BRCA_TP_clin_df_Female_samples.rds')
saveRDS(se.data.f, './BRCA_TP_SummarizedExperiment_raw.rds')


#read data
data.brca.n.f <- readRDS('./BRCA_TP_DESeq2_norm_counts_Female_samples.rds')
clin.brca.f <- readRDS('./BRCA_TP_clin_df_Female_samples.rds')
#get list on genes that have significant impact on survival based on expression quantile
ggi97.hub.df <- readRDS('~/HumanNetv3/data/GGI97_neighbor_hubgenes_DF.rds')

#normalize for tcell...since we are looking for Tcell hub genes
tcell.sig.mean <- apply(data.brca.n.f, 2, FUN = function(x){(x['CD3D'] + x['CD3E'] + x['CD3G']) / 3})

#divde each row by tcell.sig.mean, we must transpose the final result from gene * cell
data.tcell.n <- as.data.frame(t(apply(data.brca.n.f, 1, function(x) x/tcell.sig.mean)))
saveRDS(data.tcell.n, '~/HumanNetv3/TCGA/RNAseq_data_Tcell_Normalized.rds')
#data.tcell.n <- readRDS('~/HumanNetv3/TCGA/RNAseq_data_Tcell_Normalized.rds')

#double check that order match of rnaseq data and metadata
identical(colnames(data.tcell.n), clin.brca.f$submitter_id)



#cehck survival for top genes...this did not work
for (i in 1:5){
  gene <- ggi97.hub.df$names[1:(10*i)] #from 10 to 50
  #gene <- 'CD28'
  print(i)
  
  #get samples by gene expression of this gene with rna seq data
  #gene.data <- apply(data.brca.n.f[rownames(data.brca.n.f) %in% gene,],2,median)
  gene.data <- apply(data.brca.n.f[rownames(data.brca.n.f) %in% gene,],2,median)
  gene.num <- length(gene)
  

  #this works becuase rna seq and clin data are in the same order!
  #0;low 1;medium 2;high
  #BRCA has generally high survival rates...let's cut it to top 30 bottom 30 and remove middle
  clin.brca.f$gene.group <- cut(gene.data, breaks = quantile(gene.data, c(0,0.3,0.7, 1)), labels=c('Low','Med','High'), include.lowest=TRUE)
  
  #need to add numeric variables to follow-up and vital status
  data <- clin.brca.f
  rownames(data) <- clin.brca.f$submitter_id
  data$s <- grepl("dead", data$vital_status, ignore.case = TRUE)
  
  #if not dead, days_to_last_follow_up becomes days to death (incorporating censored data)
  notDead <- is.na(data$days_to_death)
  if (any(notDead == TRUE)) {
    data[notDead, "days_to_death"] <-
      data[notDead, "days_to_last_follow_up"]
  }
  
  # Column with groups
  data$type <- as.factor(data[, 'gene.group'])
  data <- data[, c("days_to_death", "s", "type")]
  
  #To exclude med group
  data <- data[data$type != "Med",]
  
  
  # create the formula for survival analysis
  fit <- survfit(Surv(days_to_death, s) ~ type, data=data)
  #get summary for 5 year
  p <- survminer::ggsurvplot(fit,
                             data=data,
                             risk.table = TRUE,
                             pval = TRUE,
                             break.time.by = 500,
                             xlim=c(0,3000),
                             #legend.title = paste0("Top_",gene.num),
                             legend.title = gene,
                             risk.table.y.text.col=T,
                             risk.table.y.text=T,
                             risk.table.col = 'strata',
                             legend.labs = c('Low', 'High'),
                             palette = c('Blue', 'Red')
  )
  p.l[[i]] <- p
  print(gene)
  
}
p


res <- arrange_ggsurvplots(p.l, print=F, ncol = 5, nrow=2)
ggsave(file = './survival_10to50.pdf', res, width = 25, height=10)


#unbiased approach shows that median of collective genes are clinically bad...we should try to interpret selected hubs based on manual curation
top.ggi97.hub <- head(ggi97.hub.df$names, 50)
#okay...to find what top LLS gene might be associated with good clincal outcome, lets find what is associated with good outcome first
#Do this for just nor counts and Tnormalized dataset
tokenStop <- 1
tabSurvKMcomplete <- NULL
for( i in 1: round(nrow(data.brca.n.f)/100)){
  message( paste( i, "of ", round(nrow(data.brca.n.f)/100)))
  tokenStart <- tokenStop
  tokenStop <-100*i
  
  # > 0.7;high, 0.3 <; Low (by qunatile)
  tabSurvKM<-TCGAanalyze_SurvivalKM(clin.brca.f,
                                    data.brca.n.f,
                                    Genelist = rownames(data.brca.n.f)[tokenStart:tokenStop],
                                    Survresult = F,
                                    ThreshTop=0.7,
                                    ThreshDown=0.3)
  
  tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
}

#tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.05,]
tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]

tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
  rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
]

prognostic.genes <- tabSurvKMcompleteDEGs
prognostic.good <- prognostic.genes[prognostic.genes[,'Group2 Deaths with Top'] < prognostic.genes[,'Group2 Deaths with Down'],]
prognostic.bad <- prognostic.genes[prognostic.genes[,'Group2 Deaths with Top'] > prognostic.genes[,'Group2 Deaths with Down'],]


adjusted_p_val_of_good_KM <- p.adjust(prognostic.good$pvalue,method="hochberg",n=length(prognostic.good$pvalue))
adjusted_p_val_of_bad_KM <- p.adjust(prognostic.bad$pvalue,method="hochberg",n=length(prognostic.bad$pvalue))

prognostic.good <- cbind(prognostic.good, adj_p_val=adjusted_p_val_of_good_KM)
prognostic.bad <- cbind(prognostic.bad, adj_p_val=adjusted_p_val_of_bad_KM)

prognostic.good <- prognostic.good[prognostic.good$adj_p_val < 0.05,]
prognostic.bad <- prognostic.bad[prognostic.bad$adj_p_val < 0.05,]

good_KM <- rownames(prognostic.good)
bad_KM <- rownames(prognostic.bad)

setwd('/home3/junhacha/HumanNetv3/TCGA/')
saveRDS(good_KM, 'Good_KM_normCounts.rds')

# for norm counts we get 1...TNFRSF18
top.ggi97.hub[top.ggi97.hub %in% good_KM]


tokenStop <- 1
tabSurvKMcomplete <- NULL
for( i in 1: round(nrow(data.tcell.n)/100)){
  message( paste( i, "of ", round(nrow(data.tcell.n)/100)))
  tokenStart <- tokenStop
  tokenStop <-100*i
  
  # > 0.7;high, 0.3 <; Low (by qunatile)
  tabSurvKM<-TCGAanalyze_SurvivalKM(clin.brca.f,
                                    data.tcell.n,
                                    Genelist = rownames(data.tcell.n)[tokenStart:tokenStop],
                                    Survresult = F,
                                    ThreshTop=0.7,
                                    ThreshDown=0.3)
  
  tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
}

#tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.05,]
tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]

tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
  rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
]

prognostic.genes <- tabSurvKMcompleteDEGs
prognostic.good <- prognostic.genes[prognostic.genes[,'Group2 Deaths with Top'] < prognostic.genes[,'Group2 Deaths with Down'],]
prognostic.bad <- prognostic.genes[prognostic.genes[,'Group2 Deaths with Top'] > prognostic.genes[,'Group2 Deaths with Down'],]


adjusted_p_val_of_good_KM <- p.adjust(prognostic.good$pvalue,method="hochberg",n=length(prognostic.good$pvalue))
adjusted_p_val_of_bad_KM <- p.adjust(prognostic.bad$pvalue,method="hochberg",n=length(prognostic.bad$pvalue))

prognostic.good <- cbind(prognostic.good, adj_p_val=adjusted_p_val_of_good_KM)
prognostic.bad <- cbind(prognostic.bad, adj_p_val=adjusted_p_val_of_bad_KM)

prognostic.good <- prognostic.good[prognostic.good$adj_p_val < 0.05,]
prognostic.bad <- prognostic.bad[prognostic.bad$adj_p_val < 0.05,]

good_KM <- rownames(prognostic.good)
bad_KM <- rownames(prognostic.bad)

setwd('/home3/junhacha/HumanNetv3/TCGA/')
saveRDS(good_KM, 'Good_KM_Tcellnorm.rds')


# for Tcell_norm counts we get none....so lets just do normCount
top.ggi97.hub[top.ggi97.hub %in% good_KM]


#get survival plot for overlapped genes
#get samples by gene expression of this gene with rna seq data

gene <- 'TNFRSF18'
#for single genes we should do this!!
gene.data <- as.numeric(data.brca.n.f[gene,])

#this works becuase rna seq and clin data are in the same order!
#0;low 1;medium 2;high
#BRCA has generally high survival rates...let's cut it to top 30 bottom 30 and remove middle
#clin.brca.f$gene.group <- cut(gene.data, breaks = quantile(gene.data, c(0,0.3,0.7, 1)), labels=c('Low','Med','High'), include.lowest=TRUE)

#lets cut it 50 high to match with GGI97 criteria
clin.brca.f$gene.group <- cut(gene.data, breaks = quantile(gene.data, c(0,0.5, 1)), labels=c('Low','High'), include.lowest=TRUE)

#need to add numeric variables to follow-up and vital status
data <- clin.brca.f
rownames(data) <- clin.brca.f$submitter_id
data$s <- grepl("dead", data$vital_status, ignore.case = TRUE)

#if not dead, days_to_last_follow_up becomes days to death (incorporating censored data)
notDead <- is.na(data$days_to_death)
if (any(notDead == TRUE)) {
  data[notDead, "days_to_death"] <-
    data[notDead, "days_to_last_follow_up"]
}

# Column with groups
data$type <- as.factor(data[, 'gene.group'])
data <- data[, c("days_to_death", "s", "type")]

#To exclude med group
data <- data[data$type != "Med",]


# create the formula for survival analysis
fit <- survfit(Surv(days_to_death, s) ~ type, data=data)
#get summary for 5 year
p <- survminer::ggsurvplot(fit,
                           data=data,
                           risk.table = TRUE,
                           pval = TRUE,
                           break.time.by = 500,
                           xlim=c(0,3000),
                           legend.title = gene,
                           risk.table.y.text.col=T,
                           risk.table.y.text=T,
                           risk.table.col = 'strata',
                           legend.labs = c('Low', 'High'),
                           palette = c('Blue', 'Red')
)
p.l <- list()
p.l[[1]] <- p
res <- arrange_ggsurvplots(p.l, print=F, ncol = 1, nrow=1,)
ggsave(file = '/home3/junhacha/HumanNetv3/Graphs/Survival_TCGA_TNFRSF18_50high.pdf', res, width = 5, height=5)



#draw network that visualize TNFRSF18 as cetner for tcell net
net.list <- readRDS('~/HumanNetv3/XC_nobatch/BC/BC_celltypeNet_list_LLS.rds')
tnet <- net.list[['T_cell']]
tnet.genes <- unique(c(tnet[,1],tnet[,2]))

library(igraph)
#get subnet that contatins TNFRSF18
aaa <- tnet[tnet[,1] %in% 'TNFRSF18' | tnet[,2] %in% 'TNFRSF18',]

#draw those that are connetected to GGI97
gipr.genes <- unique(c(aaa[,1],aaa[,2]))

tnet.sub <- graph.data.frame(aaa, directed=F)
#visualize Tnet.sub
plot.igraph(tnet.sub)
p2 <- plot.igraph(tnet.sub, asp = 0, main = "BC Tnet with GIPR",
            ## colors =======================================
            vertex.color = rgb(0.8,0.4,0.3,0.8),      ## color
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
            vertex.label.font = 2 ## 1 = plain, 2 = bold, 3 = italic, 4 = bold/italic
            
)

gipr.genes



#get graph for MKI67
bbb <- tnet[tnet[,1] %in% 'MKI67' | tnet[,2] %in% 'MKI67',]
mki67.genes <- unique(c(bbb[,1],bbb[,2]))
mki67.sub <- graph.data.frame(bbb, directed=F)

#how many of them are GGI97genes?
bc.sig <- read.table('~/HumanNetv3/BC_signaturegenes_33union.csv', header = T, sep = ',')
ggi97.genes <- bc.sig$gene[bc.sig$GGI97 == 'yes']
ggi97.genes.tnet <- as.character(ggi97.genes[ggi97.genes %in% tnet.genes])

table(ggi97.genes.tnet %in% mki67.genes)

#of 24 18 genes are GGI97 genes!!
#make group and color


#draw sub graph 
p3 <- plot.igraph(mki67.sub, asp = 0, main = "BC Tnet with MKI67",
                  ## colors =======================================
                  vertex.color = rgb(0.8,0.4,0.3,0.8),      ## color
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
                  vertex.label.font = 2 ## 1 = plain, 2 = bold, 3 = italic, 4 = bold/italic
                  
)


#lets see if we can draw them together in one network 
#get subnet that contatins TNFRSF18
ccc <- tnet[tnet[,1] %in% c('TNFRSF18','MKI67') | tnet[,2] %in% c('TNFRSF18','MKI67'),]
tnetsub2.genes <- unique(c(ccc[,1],ccc[,2]))

#draw those that are connetected to GGI97
tnet.sub2 <- graph.data.frame(ccc, directed=F)

#make GGI97gene group in the tnet
all.genes <- data.frame(genes = tnetsub2.genes)
all.genes$ggi97.label <- ifelse(all.genes$genes %in% ggi97.genes.tnet, 1, 0)

tnet.sub2 <- set_vertex_attr(tnet.sub2, "GGI97_tnet", index = V(tnet.sub2), as.character(all.genes$ggi97.label))
list.vertex.attributes(tnet.sub2)

#color
V(tnet.sub2)$color = rgb(0.8,0.4,0.3,0.8)
V(tnet.sub2)$color = ifelse(V(tnet.sub2)$GGI97_tnet == '1', "#013220", V(tnet.sub2)$color)


#make GGI97 genes bigger in size
V(tnet.sub2)$size <- V(tnet.sub2)$GGI97_tnet
normalize_01 <- function(x) (x - min(x)) / (max(x) - min(x)) + 0.8
V(tnet.sub2)$size <- normalize_01(as.numeric(V(tnet.sub2)$GGI97_tnet)) * 5

list.vertex.attributes(tnet.sub2)
#set legend
legend_cats <- data.frame(attr = unique(vertex_attr(tnet.sub2, "GGI97_tnet")),
                          color = unique(V(tnet.sub2)$color))
legend_cats <- legend_cats[order(legend_cats$attr), c(1, 2)]



#visualize Tnet.sub2
plot.igraph(tnet.sub2, asp = 0, main = "Breast Cancer Tnet with GIPR & GGI97",
                  ## colors =======================================
                  vertex.frame.color = "white",             ## border color
                  ## shapes =======================================
                  vertex.shape = "circle",                  ## none, circle, square, csquare, 
                  ## vrectangle, pie, raster, sphere
                  ## rectangle, crectangle
                  ## sizes =======================================
                  vertex.size = 3,                         ## size, default = 15
                  vertex.size2 = NA,                      ## second dimension size (for parallelograms)
                  ## color =======================================
                  vertex.label.color = "black",
                  ## font family =======================================
                  vertex.label.family = "Times",
                  ## font face =======================================
                  vertex.label.font = 1, ## 1 = plain, 2 = bold, 3 = italic, 4 = bold/italic
                  vertex.label.cex = normalize_01(as.numeric(V(tnet.sub2)$GGI97_tnet)) * 0.5
                  
)
legend(x = "bottomleft",      ## position, also takes x,y coordinates
       legend = legend_cats$attr,
       pch = 19,              ## legend symbols see ?points
       col = legend_cats$color,
       bty = "n",
       title = "MKI97 Genes in Tnet")

                    


#draw a survival plot of mean GGI97 signature expression
df <- read.table('/home3/junhacha/HumanNetv3/BC_signaturegenes_33union.csv', header = T, sep = ',')
ggi97.genes <- as.character(df$gene[df$GGI97 == 'yes'])


#how many of them are in the RNA seq data...76! probably others are deprecated
table(ggi97.genes %in% rownames(data.brca.n.f))


#use all 76 as one would in bulk, clinical practice
gene <- ggi97.genes[ggi97.genes %in% rownames(data.brca.n.f)]

#get samples by gene expression of this gene with rna seq data
gene.data <- apply(data.brca.n.f[rownames(data.brca.n.f) %in% gene,],2,median)
gene.num <- length(gene)

#double check that order match of rnaseq data and metadata
identical(colnames(data.brca.n.f), clin.brca.f$submitter_id)

#this works becuase rna seq and clin data are in the same order!
#0;low 1;medium 2;high
#divde samples top50 and bottome 50...same as GITR
clin.brca.f$gene.group <- cut(gene.data, breaks = quantile(gene.data, c(0,0.5, 1)), labels=c('Low','High'), include.lowest=TRUE)

#need to add numeric variables to followup and vital status
data <- clin.brca.f
rownames(data) <- clin.brca.f$submitter_id
data$s <- grepl("dead", data$vital_status, ignore.case = TRUE)

#if not dead, days to death becomes days_to_last_follow_up
notDead <- is.na(data$days_to_death)
if (any(notDead == TRUE)) {
  data[notDead, "days_to_death"] <-
    data[notDead, "days_to_last_follow_up"]
}

# Column with groups
data$type <- as.factor(data[, 'gene.group'])
data <- data[, c("days_to_death", "s", "type")]

#To exclude med group
data <- data[data$type != "Med",]

#exclude data with over 2000day survival...roughly a 5year survival
#data <- data[data$days_to_death <= 2000, ]


# create the formula for survival analysis
fit <- survfit(Surv(days_to_death, s)~type, data=data)
p <- survminer::ggsurvplot(fit,
                           data=data,
                           risk.table = TRUE,
                           pval = TRUE,
                           break.time.by = 500,
                           xlim=c(0,3000),
                           legend.title = "GGI97 Sig",
                           risk.table.y.text.col=T,
                           risk.table.y.text=T,
                           risk.table.col = 'strata',
                           legend.labs = c('Low', 'High'),
                           palette = c('Blue', 'Red')
)

plot.list<- list()
plot.list[[1]] <- p
res <- arrange_ggsurvplots(plot.list, print=F, ncol = 1, nrow=1)
ggsave(file = './GGI97 signature Survival TCGA.pdf', res, width = 5, height=5)




#get scatter plot to see correlation between TNFRSF18 and CD3, CD8, FOXP3
library(ggpubr)

rna.data <- readRDS('/home3/junhacha/HumanNetv3/TCGA/BRCA_TP_DESeq2_norm_counts_Female_samples.rds')
tcell.sig.mean <- apply(rna.data, 2, FUN = function(x){(x['CD3D'] + x['CD3E'] + x['CD3G']) / 3})
aaa <- as.data.frame(t(rna.data))

test <- aaa[,c('TNFRSF18', 'CD3D', 'CD8A', 'FOXP3')]
test$Tcell.markers <- tcell.sig.mean

p1 <- ggscatter(test, x = "TNFRSF18", y = "Tcell.markers", 
                add = 'reg.line', conf.int = T,
                add.params = list(color='blue', fill = 'lightgray')) +
  stat_cor(method = "pearson", label.x = 5000, label.y = 5000)

p2 <- ggscatter(test, x = "TNFRSF18", y = "CD8A",
                add = 'reg.line', conf.int = T,
                add.params = list(color='blue', fill = 'lightgray')) +
  stat_cor(method = "pearson", label.x = 5000, label.y = 5000)

p3 <- ggscatter(test, x = "TNFRSF18", y = "FOXP3",
                add = 'reg.line', conf.int = T,
                add.params = list(color='blue', fill = 'lightgray')) +
  stat_cor(method = "pearson", label.x = 5000, label.y = 2000)

library(patchwork)
pdf('/home3/junhacha/HumanNetv3/TCGA/GITR_Tcells_correlations_TCGA.pdf', width = 12, height=4)
p1+p2+p3
dev.off()







        