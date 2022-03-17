#' GenesetHnv3
#'
#' @param geneset geneset of interest. provided multipe genelist as a named list or as a character vector
#' @param reference.network default HumanNetv3
#'
#' @return a dataframe of genelist with connectivity, detected genes in the HumanNetv3 and the normalized connectivity
#' @export
#'
#' @examples
#' data("BC_signatures")
#' hnv3.connectivity.bcsig <- GenesetHnv3(geneset = bc.sig.list)
GenesetHnv3 <- function(geneset = NULL, reference.network = graph.hn3){
  humannet <- igraph::as_data_frame(graph.hn3)
  hnv3.genes <- unique(c(as.character(humannet[,1]), as.character(humannet[,2])))

  if (is.list(geneset)){
    #if the list does not have a specific name provide them as default
    if (is.null(names(geneset))){
      names(geneset) <- paste0('Geneset_',seq(length(geneset)))
      print('genesets does not have names, will be given Geneset_n')
    }


    hnv3.conn.list <- list()
    for (i in seq_along(geneset)){
      sig.genes <- geneset[[i]]
      hnv3.connectivity.sig <- nrow(humannet[(humannet[,1] %in% sig.genes & humannet[,2] %in% sig.genes), ])
      hnv3.conn.list[[i]] <- hnv3.connectivity.sig
    }


    sig.detect <- lapply(geneset, function(sig){
      length(sig[sig %in% hnv3.genes])
    })
    names(sig.detect) <- names(geneset)
    sig.detect <- unlist(sig.detect)

    hnv3.connectivity.sig <- as.data.frame(unlist(hnv3.conn.list))
    rownames(hnv3.connectivity.sig) <- names(bc.sig.list)
    hnv3.connectivity.sig$siggene_num <- lengths(bc.sig.list)
    hnv3.connectivity.sig$siggene_num_detected <- sig.detect
    hnv3.connectivity.sig$name <- rownames(hnv3.connectivity.sig)
    colnames(hnv3.connectivity.sig)[1] <- 'connectivity'

  }
  else if (is.vector(geneset)){
    hnv3.connectivity.sig <- nrow(humannet[(humannet[,1] %in% geneset & humannet[,2] %in% geneset), ])
    sig.detect <- length(geneset[geneset %in% hnv3.genes])
    hnv3.connectivity.sig <- data.frame("connectivity" = hnv3.connectivity.sig,
                                        "siggene_num" = length(geneset),
                                        "siggene_num_detected" = length(sig.detect),
                                        "name" = "User_geneset1")


  }

  return(hnv3.connectivity.sig)

}
