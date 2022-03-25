#' SortAddLLS
#'
#' @description Sort SCINET output edges alphabetically and add LLS weight derived from HumanNetv3 reference interactome
#'
#' @param Celltype.specific.networks SCINET output from function run.SCINET.clusters().
#' @param reference.network reference network input, default is HumanNetv3.
#'
#' @return List of dataframe edgelist, with gene interaction and weights from SCINET and HumanNetv3.
#' @export
#'
#' @examples
#' SortAddLLS(Celltype.specific.networks)
SortAddLLS <- function(Celltype.specific.networks = NULL, reference.network = graph.hnv3){
  sorted.net.list <- lapply(Celltype.specific.networks, function(cell.net){
    int.net <- intersection(reference.network, cell.net, keep.all.vertices=F)
    names(edge_attr(int.net)) <- c('LLS', 'scinet_weight')
    tmp <- t(apply(get.data.frame(int.net)[,1:2],1,sort))
    x <- cbind(tmp, get.data.frame(int.net)[,3:4])
  })
  return(sorted.net.list)
}


#' GetCentrality
#'
#' @description Get centrality values for each nodes of scHumanNet list. Ribosomal genes are excluded when calculating centralities
#'
#' @param method Centrality measure to calcuate, supports degree(strength), betweenness, closeness, and eigenvector
#' @param net.list Output of SortAddLLS()
#'
#' @return List of named vector, each value corresponding to node's centrality value
#' @export
#'
#' @examples
#' strenght.list <- GetCentrality(method = 'degree', net.list = Celltype.specific.networks)
GetCentrality <- function(method = method, net.list = NULL){

  cell.list <- list()
  for (cell.net in seq_along(net.list)){
    LLS.net <- net.list[[cell.net]][,c(1:3)]
    colnames(LLS.net) <- c('gene1', 'gene2', 'weight')
    net.graph <- graph_from_data_frame(LLS.net, directed=FALSE)

    if (method == 'degree'){
      genes.cent <- strength(net.graph)[order(strength(net.graph), decreasing = T)]
    }
    else if (method == 'betweenness'){
      E(net.graph)$weight <- 1
      genes.cent <- betweenness(net.graph)[order(betweenness(net.graph), decreasing = T)]
    }
    else if (method == 'closeness'){
      E(net.graph)$weight <- 1
      genes.cent <- closeness(net.graph)[order(closeness(net.graph), decreasing = T)]
    }
    else if (method == 'eigenvector'){
      eigen.cent <- eigen_centrality(net.graph, weights = E(net.graph)$weight)$vector
      genes.cent <- eigen.cent[order(eigen.cent, decreasing = T)] # change this measure
    }
    else{
      print('wrong method input')
    }

    #remove ribo.genes
    genes.all <- V(net.graph)$name
    ribo.genes <-  genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
    mito.genes <- genes.all[grep('^MRPS[0-9]*|^MRPL[0-9]*', genes.all)]
    bad.genes <- c(ribo.genes,mito.genes)
    genes.cent.f <- genes.cent[!(names(genes.cent) %in% bad.genes)]

    #to normalize for network size, we use percentile rank
    final.rank <- dplyr::percent_rank(genes.cent.f)
    cell.list[[cell.net]] <- final.rank
  }
  names(cell.list) <- names(net.list)


  #remove NAs from each cancer type
  cent.list.final = list()
  for (cell in names(cell.list)){
    data.cent <- cell.list[[cell]]
    na.slot <- is.na(names(data.cent))
    cent.list.final[[cell]] <- data.cent[!na.slot]
  }

  return(cent.list.final)
}


#' CombinePercRank
#'
#' @description Calcuate percentile rank of each centrality measure. For those not in the network but in the nodes of HNv3, 0 value is assigned.
#'
#' @param perc.rank.list Output of GetCentrality()
#'
#' @return Dataframe of celltypes and their centrality values
#' @export
#'
#' @examples
#' sorted.net.list <- SortAddLLS(Celltype.specific.networks)
#' strength.list <- GetCentrality(net.list = sorted.net.list)
#' rank.df.final <- CombinePercRank(perc.rank.list = strength.list)
CombinePercRank <- function(perc.rank.list = perc.rank.list){
  coverage.genes <- V(graph.hn3)$name
  df <- lapply(perc.rank.list, data.frame)
  df <- df %>%
    imap(~setNames(.x, .y)) %>%
    map(tibble::rownames_to_column) %>%
    reduce(full_join, by = "rowname") %>%
    mutate_all(~replace(., is.na(.), 0))

  add.genes <- coverage.genes[!(coverage.genes %in% df$rowname)]

  add.df <- data.frame(matrix(0, nrow = length(add.genes), ncol = ncol(df)))
  add.df[,1] <- add.genes
  colnames(add.df) <- colnames(df)
  rank.df.final <- rbind(df, add.df)

  rownames(rank.df.final) <- rank.df.final$rowname
  rank.df.final$rowname <- NULL

  return(rank.df.final)
}

#' TopHub
#'
#' @description Get top n genes in terms of centrality for each scHumanNet
#'
#' @param perc.rank.list Output of GetCentrality()
#'
#' @return Dataframe of celltypes and their top n centrality genes
#' @export
#'
#' @examples
#' sorted.net.list <- SortAddLLS(Celltype.specific.networks)
#' strength.list <- GetCentrality(net.list = sorted.net.list)
#' rank.df.final <- CombinePercRank(perc.rank.list = strength.list)
#' top.df <- TopHub(rank.df.final, top.n = 50)
TopHub <- function(rank.df.final = NULL, top.n = NULL){
  top.genes <- list()
  for (cell in colnames(rank.df.final)){
    perc.rank <- unlist(rank.df.final[cell])
    names(perc.rank) <- rownames(rank.df.final)
    perc.rank.top <- head(names(perc.rank[order(perc.rank, decreasing = T)]),50)
    top.genes[[cell]] <- perc.rank.top
  }

  top.df <- data.frame(matrix(ncol = 1, nrow = 50))
  for (genes in top.genes){
    top.df <- cbind(top.df,genes)
  }
  top.df[,1] <- NULL
  colnames(top.df) <- names(rank.df.final)
  return(top.df)
}



#' DiffPR
#'
#' @description Get differnence of normalized centrality values from the output of GetCentrality()
#'
#' @param rank.df.final Output of CombinePercRank
#' @param meta meta.data Dataframe that contatins annotated celltypes and condition columns
#' @param celltypes String character of column name that stores annoated celltypes
#' @param condition String character of column name that stores disease vs control
#' @param control Character string that states which of the two condition in the condition column name that will be the control
#'
#' @return Dataframe of celltypes and their diffPR values of each genes in the scHumanNet
#' @export
#'
#' @examples
#' sorted.net.list <- SortAddLLS(Celltype.specific.networks)
#' strength.list <- GetCentrality(net.list = sorted.net.list)
#' rank.df.final <- CombinePercRank(perc.rank.list = strength.list)
#' diffPR.df <- DiffPR(rank.df.final, celltypes = 'celltypes_merged', condition = 'diagnosis', control = 'Control', meta = meta)
DiffPR <- function(rank.df.final = NULL,
                   celltypes = NULL ,
                   condition = NULL,
                   control = NULL,
                   meta = NULL){

  rank.list <- list()

  for (celltype in names(table(meta[,celltypes]))){
    #print(celltype)
    conditions <- as.character(unique(meta[,condition]))
    control = conditions[conditions == control]
    disease = conditions[conditions != control]
    #subset brain percentile rank
    df <- rank.df.final[,c(paste(control, celltype, sep = '_'), paste(disease, celltype,sep = '_'))]
    #its percentile rank...ASD-CTL lets higher value have higher centrality in ASD
    df$diff.rank <- df[,2] - df[,1]

    #make a list for each celltype get genes that have the most rank differential
    rank.list[[celltype]] <- df$diff.rank
    names(rank.list[[celltype]]) <- rownames(df)
    rank.list[[celltype]] <- rank.list[[celltype]][order(abs(rank.list[[celltype]]), decreasing = T)]
  }

  #check that all element of list contain same genes(all genes)
  df.final <- data.frame(matrix(nrow=length(rank.list[[1]]), ncol=2*length(rank.list)))

  for (i in seq(1,length(rank.list))){
    t <- 2*i - 1
    df.final[,t] <- names(rank.list[[i]])
    df.final[,t+1] <- rank.list[[i]]
  }

  column.name <- names(rank.list)
  final.column.names <- vector()
  for (i in seq_along(column.name)){
    celltype <- column.name[i]
    column.name2 <- paste(column.name[[i]], paste0(disease,'-',control),sep = '_')
    index = 2*i -1
    final.column.names[index] <- celltype
    final.column.names[index+1] <- column.name2
  }
  colnames(df.final) <- final.column.names

  return(df.final)

}


#' FindDiffHub
#'
#' @description Calcuate statistical signficance of diffPR values calculated with DiffPR()
#'
#' @param rank.df.final Output from CombinePercRank
#' @param meta metadata data.frame used for this analayiss
#' @param celltypes column name where celltypes annotaiton are stored
#' @param condition column name where disase and control annotaition is stored
#' @param control character string that specific which value is the control
#' @param net.list Output of SortAddLLS, network list to perfrom the analysis
#' @param centrality the centrality method used in the getCentrality. both must match! this is used to find the null distribution
#' @param q.method Method to perfrom multiple testing correction. accepted values are: c("BH","holm","hochberg","bonferroni","BY","fdr","none")
#'
#' @return Dataframe Percentile rank of centrality of each networks, their difference, and their statistical significance
#' @export
#'
#' @examples
#' sorted.net.list <- SortAddLLS(Celltype.specific.networks)
#' strength.list <- GetCentrality(net.list = sorted.net.list, method = 'degree')
#' rank.df.final <- CombinePercRank(perc.rank.list = strength.list)
#' diffPR.df <- DiffPR(rank.df.final, celltypes = 'celltypes_merged', condition = 'diagnosis', control = 'Control', meta = meta)
#' diffPR.df.sig <- FindDiffHub(rank.df.final = rank.df.final, celltypes = 'celltypes_merged', condition = 'diagnosis', control = 'Control', meta = meta, net.list=sorted.net.list, q.method='BH', centrality="degree")
FindDiffHub <- function(rank.df.final = NULL,
                        celltypes = NULL ,
                        condition = NULL,
                        control = NULL,
                        meta = NULL, 
                        net.list = NULL,
                        centrality = "degree",
                        q.method = "BH"){
  
  final.df.list <- list()
  for (celltype in names(table(meta[,celltypes]))){
    #progress bar
    print(paste0("Finding DiffHubs in ",celltype,"..."))

    conditions <- as.character(unique(meta[,condition]))
    control = conditions[conditions == control]
    disease = conditions[conditions != control]
    
    #get diffPR.df of ctrl disase for each celltype
    df <- rank.df.final[,c(paste(control, celltype, sep = '_'), paste(disease, celltype,sep = '_'))]
  
    #get disease net
    disease.net.name <- paste(disease, celltype,sep = '_')
    disease.net <- igraph::graph_from_data_frame(net.list[[disease.net.name]], directed = F)
    shuffled.weight1 <- sample(E(disease.net)$LLS) #rewire does not suporte weight..so we are going to shuffle both node and edges, while preserving topology
    
    #get control net
    control.net.name <- paste(control, celltype,sep = '_')
    control.net <- igraph::graph_from_data_frame(net.list[[control.net.name]], directed = F)
    shuffled.weight2 <- sample(E(control.net)$LLS) #rewire does not suporte weight..so we are going to shuffle both node and edges, while preserving topology
    
    
    #get df where at least one column has non-zero value. we don't need to evaluate diffPR zero..
    df.f <- df[!(df[,1] ==0 & df[,2] == 0),]
    #write.table(df.f, paste0('~/test/df_f_', celltype,'.tsv'), sep = '\t', quote = F, row.names = T, col.names = T)
    
    #its percentile rank...ASD-CTL lets higher value have higher centrality in ASD
    df.f$gene <- rownames(df.f)
    df.f$DiffPR <- df.f[,2] - df.f[,1]
    
    
    #make this once for all gene comparision, we assume that null diffPR distribution is same for all genes
    #to get significance of our diffPR value, we make n random diffPR and compare
    null.distribution <- vector()
    #make two random network and append diffPR values to null.distribution
    while (length(null.distribution) < 1000000){
      #random disease network
      random.disease <- rewire(disease.net, with = each_edge(1, loops = F))
      E(random.disease)$LLS <- shuffled.weight1
      random.net.disease <- igraph::as_data_frame(random.disease)
      random.cent.disease <- GetCentrality(method = centrality, net.list =list(disease.net.name = random.net.disease))
      random.cent.disease <- as.data.frame(random.cent.disease[[1]])
      random.cent.disease$gene <- rownames(random.cent.disease)
      
      #random control network
      random.control <- rewire(disease.net, with = each_edge(1, loops = F))
      E(random.control)$LLS <- shuffled.weight2
      random.net.control <- igraph::as_data_frame(random.control)
      random.cent.control <- GetCentrality(method = centrality, net.list =list(control.net.name = random.net.control))
      random.cent.control <- as.data.frame(random.cent.control[[1]])
      random.cent.control$gene <- rownames(random.cent.control)
      
      #bind two centrality
      diffPR.random <- dplyr::full_join(random.cent.disease, random.cent.control, by='gene')
      diffPR.random[is.na(diffPR.random)] <- 0 #replace NA with 0
      diffPR.values <- as.numeric(diffPR.random[,1]) - as.numeric(diffPR.random[,3])
      null.distribution <- c(null.distribution, diffPR.values)
    }


    
    #remove RPS RPL MRPS MRPL from union geneset of df.f to save time..these genes will be disregarded
    genes.all <- df.f$gene
    ribo.genes <-  genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
    mito.genes <- genes.all[grep('^MRPS[0-9]*|^MRPL[0-9]*', genes.all)]
    nduf.genes <- genes.all[grep('^NDUF', genes.all)]
    bad.genes <- c(ribo.genes,mito.genes, nduf.genes)
    
    df.f.f <- df.f[!(df.f$gene %in% bad.genes),]
    
    # iterate over union set of genes and get pvalue
    pb = txtProgressBar(min = 0, max = nrow(df.f), style = 3)
    pvalue.all <- vector()
    for (g in 1:nrow(df.f.f)){
      Sys.sleep(0.05)
      setTxtProgressBar(pb,g)
      gene <- rownames(df.f)[g]
      diffPR.gene <- df.f$DiffPR[g]

      #get pvalue, $frank is 10times faster
      distribution.all <- c(null.distribution, diffPR.gene)
      pvalue <- data.table::frank(-abs(distribution.all), ties.method = "min")[length(distribution.all)] / length(distribution.all) # ties are averaged..highly unlikely here
      pvalue.all[g] <- pvalue
    }
    close(pb)
    #make pvalue fo reach celltype
    df.f.f$pvalue <- pvalue.all
    df.f.f$qvalue <- p.adjust(df.f.f$pvalue, method = q.method, n = nrow(df.f.f))
    df.f.f$celltype <- rep(celltype, nrow(df.f.f))
  
  }
  
  final.df.list[[celltype]] <- df.f
  diffPR.df.result <- as.data.frame(do.call("rbind", final.df.list))
  return(dffPR.df.result)
}
  



