#' Subsets rows/columns of an ACTIONetExperiment (ACE) object
#'
#' @param i,j rows and columns
#'
#' @export
setMethod("[", c("ACTIONetExperiment", "ANY", "ANY"), function(x, i, j, ..., drop = TRUE) {
    rn <- x@rowNets
    cn <- x@colNets
    rf <- x@rowMaps
    cf <- x@colMaps
    
    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i, rownames(x), 
                fmt)
        }
        i <- as.vector(i)
        
        if (length(rn) > 0) {
            for (k in 1:length(rn)) {
                tmp = rn[[k]]
                rn[[k]] = tmp[i, i]
            }
        }
        if (length(rf) > 0) {
            for (k in 1:length(rf)) {
                tmp = rf[[k]]
                rf[[k]] = tmp[i, , drop = FALSE]
            }
        }
    }
    
    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(j, colnames(x), 
                fmt)
        }
        j <- as.vector(j)
        
        if (length(cn) > 0) {
            for (k in 1:length(cn)) {
                tmp = cn[[k]]
                cn[[k]] = tmp[j, j]
            }
        }
        if (length(cf) > 0) {
            for (k in 1:length(cf)) {
                tmp = cf[[k]]
                cf[[k]] = tmp[j, , drop = FALSE]
            }
        }
    }
    
    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, rowNets = rn, colNets = cn, rowMaps = rf, colMaps = cf, 
        check = FALSE)
})
