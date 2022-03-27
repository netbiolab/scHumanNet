#' Title Get top percent of diffHubs default is 5percent
#'
#' @param diffPR.df output dataframe of diffPR(),
#' @param top.percent top percentage of diffPR hub genes to retrieve, default is 0.05
#'
#' @return dataframe of top percentage diffPR genes for each celltype
#' @export
#'
#' @examples
#' sorted.net.list <- SortAddLLS(Celltype.specific.networks)
#' strength.list <- GetCentrality(net.list = sorted.net.list, method = 'degree')
#' rank.df.final <- CombinePercRank(perc.rank.list = strength.list)
#' diffPR.df <- DiffPR(rank.df.final, celltypes = 'celltypes_merged', condition = 'diagnosis', control = 'Control', meta = meta)
#' top.diffPR.df <- TopDiffHub(diffPR.df)
TopDiffHub <- function(diffPR.df = NULL, top.percent = 0.05){
  diffPR.df.list <- list()
  for (i in seq(1,ncol(diffPR.df),2)){
    index <- (i+1) / 2
    celltype <- colnames(diffPR.df)[i]
    #make diffPR value with gene names
    diffPR <- as.vector(diffPR.df[,(i+1)])
    names(diffPR) <- diffPR.df[,i]

    diffPR.nonzero <- diffPR[diffPR != 0]
    rank <- seq(1,length(diffPR.nonzero),1)
    top <- rank / length(diffPR.nonzero) #top ranked diffPR genes

    diffPR.df <- cbind(names(diffPR.nonzero),diffPR.nonzero, top, rep(celltype, length(diffPR.nonzero)))

    #sort by diffPRvalue
    diffPR.df <- diffPR.df[order(diffPR.df[,2]),]

    diffPR.df.list[[index]] <- diffPR.df

    print(paste(celltype, 'network:', length(diffPR.nonzero), 'nodes'))

    #calculate p-value based on absolute rank (we will consider up to 0.01)

  }
  diffPR.df.final <- as.data.frame(do.call("rbind", diffPR.df.list))
  colnames(diffPR.df.final) <- c('gene','diffPR','top_percentage','celltype')

  #change numeric values to numeric class
  diffPR.df.final[c('diffPR','top.percentage')] <- sapply(diffPR.df.final[c('diffPR','top_percentage')], function(x){as.numeric(as.character(x))})

  #get genes that have less then 0.05
  threshold <- top.percent
  diffPR.df.final.pvalue <- diffPR.df.final[diffPR.df.final$top_percentage < threshold,]

  return(diffPR.df.final.pvalue)
}
