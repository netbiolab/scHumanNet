#' Find statistically significant hub genes in each scHumanNets
#'
#' @param net.list output of SortAddLLS. list of network dataframe
#' @param centrality default is degree, which is calculated by sum of edge weights
#' @param q.method default BH, input one of the following c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
#' @param threshold a threshold value to cut the significant hub genes. default 0.05
#'
#' @return a dataframe with Percentile Rank Centrality, gene, pvalue, qvalue and the celltype
#' @export
#'
#' @examples
#' sorted.net.list <- SortAddLLS(Celltype.specific.networks, reference.network = graph.hn3)
#' sig.hub.df <- FindAllHub(net.list = sorted.net.list)
FindAllHub <- function(net.list = NULL,
                        centrality = "degree",
                        q.method = "BH",
                        threshold = 0.05){

  # run this for all celltypes in the net.list
  final.df.list <- list()
  centrality.list <- GetCentrality(net.list, method='degree')
  for (celltype in names(net.list)){
    #progress bar
    print(paste0("Finding Hubs in ",celltype,"..."))

    #get df of absolute centrality value for each network
    #get union geneset dataframe and fill with 0
    celltype.net <- igraph::graph_from_data_frame(net.list[[celltype]], directed=F)
    celltype.cent <- igraph::strength(celltype.net, weights = E(celltype.net)$LLS)
    celltype.pr <- centrality.list[[celltype]]

    #make this once for all gene comparision, we assume that null diffPR distribution is same for all genes
    #to get significance of our PR value, we make n random PR and compare
    #make ~100K random distribution of
    null.distribution.celltype <- vector()
    while (length(null.distribution.celltype) < 10000){

      #randomly shuffle non.celltype network for null model
      celltype.net <- igraph::graph_from_data_frame(net.list[[celltype]], directed = F)
      shuffled.weight <- sample(E(celltype.net)$LLS) #rewire does not suport weight..so we are going to shuffle both node and edges, while preserving topology

      #lower probabilities will adhere a longer tail and a more strict null
      random.net <- rewire(celltype.net, with = each_edge(0.5, loops = F))
      #random.net <- rewire(celltype.net, with = keeping_degseq(niter = vcount(celltype.net * 10)))
      E(random.net)$LLS <- shuffled.weight
      random.cent <- igraph::strength(random.net, weights = E(random.net)$LLS)

      null.distribution.celltype <- c(null.distribution.celltype, random.cent)
    }

    #hist(null.distribution.celltype, col='grey')


    # iterate over the output of Getcentrality only
    pb = txtProgressBar(min = 0, max = length(celltype.pr), style = 3)
    pvalue.all <- vector()
    for (g in 1:length(celltype.pr)){
      Sys.sleep(0.05)
      setTxtProgressBar(pb,g)
      gene <- names(celltype.pr)[g]
      gene.cent <- celltype.cent[gene]

      #get pvalue, $frank is 10times faster
      distribution.all <- c(null.distribution.celltype, gene.cent)
      pvalue <- data.table::frank(-distribution.all, ties.method = "min")[length(distribution.all)] / length(distribution.all) #-value to inverse rank
      pvalue.all[g] <- pvalue
    }
    close(pb)
    #make pvalue fo reach celltype
    df.f <- as.data.frame(celltype.pr)
    df.f$gene <- rownames(df.f)
    df.f$pvalue <- pvalue.all
    df.f$qvalue <- p.adjust(df.f$pvalue, method = q.method, n = nrow(df.f))
    df.f$celltype <- rep(celltype, nrow(df.f))
    #df.f$qvalue2 <- p.adjust(df.f$pvalue, method = 'bonferroni', n = nrow(df.f))


    #filter to signficant genes by qvalue
    df.f.f <- df.f[df.f$qvalue < threshold,]

    #unify the column name for rbind
    colnames(df.f.f)[1] <- 'Centrality_PR'
    dim(df.f.f)

    final.df.list[[celltype]] <- df.f.f
  }

  names(final.df.list) <- names(net.list)
  bind.df.result <- do.call("rbind", final.df.list)

  return(bind.df.result)
}
