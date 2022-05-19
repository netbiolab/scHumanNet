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
    nduf.genes <- genes.all[grep('^NDUF', genes.all)]
    bad.genes <- c(ribo.genes, mito.genes, nduf.genes)
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

  #convert all factor column to character
  i <- sapply(meta, is.factor)
  meta[i] <- lapply(meta[i], as.character)

  rank.list <- list()

  for (celltype in names(table(meta[,celltypes]))){
    #print(celltype)
    conditions <- as.character(unique(meta[,condition]))
    control = conditions[conditions == control]
    disease = conditions[conditions != control]

    #subset df
    #either celltype_condition or condition_celltype
    celltype_condition_cols <- c(paste(control, celltype, sep = '_'), paste(disease, celltype,sep = '_'), paste(celltype, control, sep = '_'), paste(celltype, disease,sep = '_'))
    colnames.in <- colnames(rank.df.final)[colnames(rank.df.final) %in% celltype_condition_cols]
    df <- rank.df.final[,colnames.in]

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
#' @param meta metadata data.frame used for this analayis
#' @param celltypes column name where celltypes annotaiton are stored
#' @param condition column name where disase and control annotaition is stored
#' @param control character string that specific which value is the control
#' @param net.list Output of SortAddLLS, network list to perfrom the analysis
#' @param centrality the centrality method used in the getCentrality. both must match! this is used to find the null distribution
#' @param q.method Method to perfrom multiple testing correction. accepted values are: c("BH","holm","hochberg","bonferroni","BY","fdr","none")
#' @param min.cells minimum cells required to perform diff network analysis default value is 500
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
                        q.method = "BH",
                        min.cells = 500){


  #check each parameter input
  if (!(celltypes %in% names(meta))){
    print(paste('celltypes column', celltypes, 'does not exist in metadata'))
    stop('this column should only have celltype names, not conditions...')
  }

  if (!(condition %in% names(meta))){
    print(paste('condition column', condition, 'does not exist in metadata'))
    stop('this column should only have conditions')
  }

  if (!(control %in% names(meta))){
    print(paste('condition', control, 'does not exist in condition column'))
    stop('this column should only have conditions')
  }

  if (!(q.method %in% c("BH","holm","hochberg","bonferroni","BY","fdr","none"))){
    stop('wrong input for Q-value method')
  }


  #set 1 core for data.table frank
  data.table::setDTthreads(threads = 1)

  #convert all factor column to character
  i <- sapply(meta, is.factor)
  meta[i] <- lapply(meta[i], as.character)

  final.df.list <- list()

  #check that each celltyppe pass minimum threshold
  celltypes.analyze <- vector()
  for (celltype in names(table(meta[,celltypes]))){
    control.cells <- meta[(meta[,celltypes] == celltype & meta[,condition] == control),]
    disease.cells <- meta[(meta[,celltypes] == celltype & meta[,condition] == disease),]

    if (nrow(control.cells) >= min.cells & nrow(disease.cells) >= min.cells){
      celltypes.analyze <- c(celltypes.analyze, celltype)

    }
    else{
      print(paste('not enough cells in', celltype,', excluding diffHub analysis...'))

    }
  }



  for (celltype in celltypes.analyze){
    #progress bar
    print(paste0("Finding DiffHubs in ",celltype,"..."))

    #get disease and control variables
    conditions <- as.character(unique(meta[,condition]))
    control = conditions[conditions == control]
    disease = conditions[conditions != control]

    #get diffPR.df of ctrl disase for each celltype
    celltype_condition_cols <- c(paste(control, celltype, sep = '_'), paste(disease, celltype,sep = '_'), paste(celltype, control, sep = '_'), paste(celltype, disease,sep = '_'))
    colnames.in <- colnames(rank.df.final)[colnames(rank.df.final) %in% celltype_condition_cols]
    df <- rank.df.final[,colnames.in]

    #get disease net # this is not needed in the null distribution generation
    disease.net.name <- colnames.in[grep(disease, colnames.in)]
    disease.net <- igraph::graph_from_data_frame(net.list[[disease.net.name]], directed = F)
    shuffled.weight1 <- sample(E(disease.net)$LLS) #rewire does not suporte weight..so we are going to shuffle both node and edges, while preserving topology

    #get control net
    control.net.name <- colnames.in[!(colnames.in %in% disease.net.name)]
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
      # we make two random network from control..and find null diffPR values

      #random control network 1
      random.control1 <- rewire(control.net, with = each_edge(0.9, loops = F))
      #random.disease <- rewire(control.net, with = keeping_degseq(niter = vcount(control.net * 50)))
      E(random.control1)$LLS <- shuffled.weight2
      random.net.control1 <- igraph::as_data_frame(random.control1)
      random.cent.control1 <- GetCentrality(method = centrality, net.list =list(disease.net.name = random.net.control1))
      random.cent.control1 <- as.data.frame(random.cent.control1[[1]])
      random.cent.control1$gene <- rownames(random.cent.control1)

      #random control network
      random.control2 <- rewire(control.net, with = each_edge(0.9, loops = F))
      #random.control <- rewire(control.net, with = keeping_degseq(niter = vcount(control.net * 50)))
      E(random.control2)$LLS <- shuffled.weight2
      random.net.control2 <- igraph::as_data_frame(random.control2)
      random.cent.control2 <- GetCentrality(method = centrality, net.list =list(control.net.name = random.net.control2))
      random.cent.control2 <- as.data.frame(random.cent.control2[[1]])
      random.cent.control2$gene <- rownames(random.cent.control2)

      #bind two centrality
      diffPR.random <- dplyr::full_join(random.cent.control1, random.cent.control2, by='gene')
      diffPR.random[is.na(diffPR.random)] <- 0 #replace NA with 0
      diffPR.values <- as.numeric(diffPR.random[,1]) - as.numeric(diffPR.random[,3])
      null.distribution <- c(null.distribution, diffPR.values)
    }

    #hist(null.distribution, col='grey')

    #remove RPS RPL MRPS MRPL from union geneset of df.f to save time..these genes will be disregarded
    genes.all <- df.f$gene
    ribo.genes <-  genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
    mito.genes <- genes.all[grep('^MRPS[0-9]*|^MRPL[0-9]*', genes.all)]
    nduf.genes <- genes.all[grep('^NDUF', genes.all)]
    bad.genes <- c(ribo.genes, mito.genes, nduf.genes)

    df.f.f <- df.f[!(df.f$gene %in% bad.genes),]

    # iterate over union set of genes and get pvalue
    pb = txtProgressBar(min = 0, max = nrow(df.f.f), style = 3)
    pvalue.all <- vector()
    for (g in 1:nrow(df.f.f)){
      Sys.sleep(0.05)
      setTxtProgressBar(pb,g)
      gene <- rownames(df.f)[g]
      diffPR.gene <- df.f.f$DiffPR[g]

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
    colnames(df.f.f) <- c("Control_scHumanNet", "Disease_scHumanNet", "gene", "diffPR", "pvalue", "qvalue", "celltype")

    final.df.list[[celltype]] <- df.f.f
  }

  names(final.df.list) <- NULL
  diffPR.df.result <- do.call("rbind", final.df.list)

  return(diffPR.df.result)
}




