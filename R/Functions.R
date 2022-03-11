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
    LLS.net <- net.list[[cell.net]][,c('1','2','LLS')]
    colnames(LLS.net) <- c('gene1', 'gene2', 'weight')
    net.graph <- graph_from_data_frame(LLS.net, directed=FALSE)

    if (method == 'degree'){
      genes.cent <- strength(net.graph)[order(strength(net.graph), decreasing = T)]
    }
    else if (method == 'betweenness'){
      genes.cent <- betweenness(net.graph)[order(betweenness(net.graph), decreasing = T)]
    }
    else if (method == 'closeness'){
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
    genes.cent.ribo.r <- genes.cent[!(names(genes.cent) %in% ribo.genes)]

    #to normalize for network size, we use percentile rank
    final.rank <- dplyr::percent_rank(genes.cent.ribo.r)
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

  for (i in seq(1,8)){
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
#' @param diffPR.df Output form diffPR()
#' @param p.value Thresholded pvalue, e.g. 0.05 or 0.01
#'
#' @return Dataframe of celltypes and their centrality values
#' @export
#'
#' @examples
#' sorted.net.list <- SortAddLLS(Celltype.specific.networks)
#' strength.list <- GetCentrality(net.list = sorted.net.list)
#' rank.df.final <- CombinePercRank(perc.rank.list = strength.list)
#' diffPR.df <- DiffPR(rank.df.final, celltypes = 'celltypes_merged', condition = 'diagnosis', control = 'Control', meta = meta)
#' diffPR.df.sig <- DiffPR(diffPR.df = diffPR.df, p.value = 0.05)
FindDiffHub <- function(diffPR.df = NULL, p.value = NULL){
  diffPR.df.list <- list()
  for (i in seq(1,16,2)){
    index <- (i+1) / 2
    celltype <- colnames(df.final)[i]
    #make diffPR value with gene names
    diffPR <- as.vector(df.final[,i+1])
    names(diffPR) <- df.final[,i]

    diffPR.nonzero <- diffPR[diffPR != 0]
    rank <- seq(1,length(diffPR.nonzero),1)
    pvalue <- rank / length(diffPR.nonzero) #non paramatric pvalue calculation

    diffPR.df <- cbind(names(diffPR.nonzero),diffPR.nonzero,pvalue, rep(celltype, length(diffPR.nonzero)))

    #sort by diffPRvalue
    diffPR.df <- diffPR.df[order(diffPR.df[,2]),]

    diffPR.df.list[[index]] <- diffPR.df

    print(paste(celltype, 'network:', length(diffPR.nonzero), 'nodes'))

    #calculate p-value based on absolute rank (we will consider up to 0.01)

  }
  diffPR.df.final <- as.data.frame(do.call("rbind", diffPR.df.list))
  colnames(diffPR.df.final) <- c('gene','diffPR','pvalue','celltype')

  #change numeric values to numeric class
  diffPR.df.final[c('diffPR','pvalue')] <- sapply(diffPR.df.final[c('diffPR','pvalue')], function(x){as.numeric(as.character(x))})

  #get genes that have less then 0.01
  pvalue <- p.value
  diffPR.df.final.pvalue <- diffPR.df.final[diffPR.df.final$pvalue < pvalue,]

  return(diffPR.df.final.pvalue)
}



