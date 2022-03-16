#' DeconvoluteNet
#'
#' @description Deconvolution of geneset connectivity for each scHumanNet constructed
#'
#' @param network.list network list of edgelist, this is the output of SortAddLLS
#' @param geneset character vector of geneset it can also be a list of multiple geneset
#'
#' @return dataframe of that stores absolute connectivity, connectivty normalized for number of detected geneset in each scHumanNet, and the celltype of scHumanNet. Note: Normalized connectivity should only be used for multiple genesets comparison, and is the value of connectivty devided by the detected number of geneset within each scHumanNet
#' @export
#'
#' @examples
#' data("ICMs")
#' DeconvoluteNet(network.list = sorted.net.list, geneset = icm.genes)
DeconvoluteNet <- function(network.list = NULL, geneset = NULL){

  #get node list for each celltype
  node.list <- lapply(network.list, function(net){
    genes.all <- unique(c(as.character(net[,1]),as.character(net[,2])))
  })
  names(node.list) <- names(network.list)


  #check if geneset is a list
  if (is.list(geneset)){
    #if the list does not have a specific name provide them as default
    if (is.null(names(geneset))){
      names(geneset) <- paste0('Geneset_',seq(length(geneset)))
    }

    sig.list <- list()
    for (i in seq_along(geneset)){
      sig.genes <- geneset[[i]]
      connectivity.sig <- lapply(network.list, function(net){
        nrow(net[(net[,1] %in% sig.genes & net[,2] %in% sig.genes), ])
      })
      names(connectivity.sig) <- names(network.list)
      connectivity.sig <- dplyr::bind_rows(connectivity.sig)
      sig.list[[i]] <- connectivity.sig
    }

    connectivity.sig.all <- as.data.frame(dplyr::bind_rows(sig.list))
    rownames(connectivity.sig.all) <- names(geneset)

    df <- t(connectivity.sig.all)

    # Gathering data, rearragne datraframe
    data <- reshape::melt(df)
    colnames(data) <- c('scHumanNet', 'signature_name', 'connectivity')

    #add gene sig length
    data$signature_gene_num <- lengths(bc.sig.list[as.character(data$signature_name)])

    #add geneset length detected in each celltype net
    detected.sig <- list()
    for (i in seq_along(table(data$signature_name))){
      signature <- geneset[[i]]
      sig.f <- lapply(node.list, function(node){length(signature[signature %in% node])})
      names(sig.f) <- names(node.list)
      detected.sig[[i]] <- sig.f
    }

    names(detected.sig) <- names(table(data$signature_name))
    detected.sig <- unlist(detected.sig)
    data$detected.sig.num <- detected.sig

    #normalize connectivity count values by number of signature genes
    #data$connectivity.n <- data$connectivity / data$signature_gene_num
    data$connectivity.n. <- data$connectivity / data$detected.sig.num

  }
  else{
    print('Only one geneset detected')
    sig.list <- list()
    sig.genes <- geneset
    connectivity.sig <- lapply(network.list, function(net){
      nrow(net[(net[,1] %in% sig.genes & net[,2] %in% sig.genes), ])
      })
    connectivity.sig <- dplyr::bind_rows(connectivity.sig)
    sig.list[[1]] <- connectivity.sig

    connectivity.sig.all <- as.data.frame(dplyr::bind_rows(sig.list))
    df <- t(connectivity.sig.all)

    # Gathering data, rearragne datraframe
    data <- as.data.frame(df)
    colnames(data) <- 'connectivity'

    #add gene sig length
    data$signature_gene_num <- rep(length(geneset), nrow(data))

    #add genes detected within the geneset in each celltype net
    sig.f <- lapply(node.list, function(node){length(geneset[geneset %in% node])})
    detected.sig <- unlist(sig.f)
    data$detected.sig.num <- detected.sig #user can know how many of their genesets were actually in the scHumanNets, but normalization does not happen here becuase detected num also is part of the connectivity measurement

    #normalize connectivity count values by number of signature genes....
    #data$connectivity.n <- data$connectivity / data$signature_gene_num
    #data$connectivity.n.detected <- data$connectivity / data$detected.sig.num ...this is not informative in single geneset undection should also be regarded

  }

  return(data)
}


#' Connectivity
#'
#' @description Deconvolution of geneset connectivity for each scHumanNet constructed
#'
#' @param network dataframe edgelist of scHumanNet. Must input only one network
#' @param geneset character vector of geneset
#' @param simulate.num number to simulate the for the null distribution, default is 10,000 times
#'
#' @return A list. First element the null distribution of connectivity from rejection sampling, second element non-parametric pvalue of your geneset, third element the detected genes used for analysis
#' @export
#'
#' @examples
#' data("ICMs")
#' Connectivity(network.list = sorted.net.list, geneset = icm.genes)
Connectivity <- function(network = NULL, geneset = NULL, simulate.num = 10000){

  # conserve similar degree for randomly selected
  network.graph <- igraph::graph_from_data_frame(network)
  degree.centrality <- igraph::degree(network.graph)

  #only use detected number of input geneset in the network
  if (isTRUE(geneset %in% names(degree.centrality))){
    print('all genes in the network')
  }
  else {
    detected.genes <- geneset[geneset %in% names(degree.centrality)]
    print(paste(length(detected.genes), 'genes detected out of', length(geneset), 'input genes'))
    print(paste('connectivity of', length(detected.genes), 'will be assessed...'))
  }


  #for progress bar
  pb = txtProgressBar(min = 0, max = simulate.num, style = 3)

  connectivity.random <- vector()
  for (n in 1:simulate.num){
    Sys.sleep(0.05)
    setTxtProgressBar(pb,n)


    #pick random genes meeting criteria
    random.geneset <- vector()
    selected <- vector()
    for (i in seq_along(detected.genes)){
      gene <- detected.genes[i]
      gene.degree <- degree.centrality[gene]

      #this is to ensure that we don't pick the same genes randomly
      degree.centrality.f <- degree.centrality[!(names(degree.centrality) %in% selected)]

      #pick a pool of nodes that have within +-10 percent of nodes
      degree.range <- c(floor(0.9 * gene.degree), ceiling(1.1 * gene.degree))
      node.pool <- degree.centrality.f[degree.range[1] <= degree.centrality.f | degree.centrality.f >= degree.range[2]]

      #pick a random node from this pool
      gene.random <- sample(names(node.pool), 1)

      #store variable
      random.geneset[i] <- gene.random
      selected[i] <- gene.random
    }

    connectivity.random.n <- nrow(network[(network[,1] %in% random.geneset & network[,2] %in% random.geneset), ])
    connectivity.random[n] <- connectivity.random.n
  }
  close(pb)


  #get pvalue of user geneset connectivty and the connectiviity.random
  connectivity <- nrow(network[(network[,1] %in% geneset & network[,2] %in% geneset), ])
  connectivity.final <- c(connectivity.random, connectivity)
  pvalue <- rank(-(connectivity.final),ties.method = 'last')[simulate.num + 1] / simulate.num
  #distribution
  output.list <- list(null.distribution = connectivity.random, p.value = pvalue, detected.geneset = detected.genes)

  return(output.list)
}



