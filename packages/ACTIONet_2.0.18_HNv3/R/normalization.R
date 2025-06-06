#' @export
normalize.scran <- function(
  ace,
  batch_attr = NULL,
  assay_name = "counts",
  BPPARAM = SerialParam()
) {

    .check_and_load_package(c("scran", "scater"))
    batch_attr = .get_attr_or_split_idx(ace, batch_attr, return_vec = TRUE)

    ace = scran::computeSumFactors(
      ace,
      clusters = batch_attr,
      assay.type = assay_name,
      BPPARAM = BPPARAM
    )

    ace = scater::logNormCounts(ace, exprs_values = assay_name)

    return(ace)
}

#' @export
normalize.multiBatchNorm <- function(
  ace,
  batch_attr,
  assay_name = "counts",
  BPPARAM = SerialParam()
) {

    .check_and_load_package(c("scran", "batchelor"))
    batch_attr = .get_attr_or_split_idx(ace, batch_attr, return_vec = TRUE)
    sce_temp = as.SingleCellExperiment(ace)

    sce_temp = batchelor::multiBatchNorm(
      sce_temp,
      batch = batch_attr,
      assay.type = assay_name,
      BPPARAM = BPPARAM
    )

    SummarizedExperiment::assays(ace)[["logcounts"]] = SummarizedExperiment::assays(sce_temp)[["logcounts"]]

    return(ace)
}

#' @export
normalize.Linnorm <- function(
  ace,
  assay_name = "counts"
) {

    .check_and_load_package("Linnorm")
    SummarizedExperiment::assays(ace)[["logcounts"]] = Linnorm(SummarizedExperiment::assays(ace)[[assay_name]])
    return(ace)
}

#' @export
normalize.default <- function(
  ace,
  assay_name = "counts",
  log_scale = TRUE
) {

    S = SummarizedExperiment::assays(ace)[[assay_name]]
    B = rescale.matrix(S, log_scale)
    rownames(B) = rownames(ace)
    colnames(B) = colnames(ace)
    SummarizedExperiment::assays(ace)[["logcounts"]] = B
    return(ace)
}

rescale.matrix <- function(
  S,
  log_scale = FALSE
) {

    if (is.matrix(S)) {
        cs = fastColSums(S)
        cs[cs == 0] = 1
        B = median(cs) * Matrix::t(Matrix::t(S) / cs)
        if (log_scale == TRUE) {
            B = log1p(B)
        }
    } else {
        A = as(S, "dgTMatrix")
        cs = fastColSums(S)
        cs[cs == 0] = 1
        x = median(cs) * (A@x/cs[A@j + 1])
        if (log_scale == TRUE) {
            x = log1p(x)
        }
        B = Matrix::sparseMatrix(
          i = A@i + 1,
          j = A@j + 1,
          x = x,
          dims = dim(A)
        )
    }

    return(B)
}

#' @export
normalize.ace <- function(
  ace,
  norm_method = "default",
  batch_attr = NULL,
  assay_name = "counts",
  BPPARAM = SerialParam()
) {

    if (norm_method == "scran") {

        ace = normalize.scran(
          ace = ace,
          batch_attr = batch_attr,
          assay_name = assay_name,
          BPPARAM = BPPARAM
        )

    } else if (norm_method == "multiBatchNorm") {

        ace = normalize.multiBatchNorm(
          ace = ace,
          batch_attr = batch_attr,
          assay_name = assay_name,
          BPPARAM = BPPARAM
        )

    } else if (norm_method == "linnorm") {

        ace = normalize.Linnorm(
          ace = ace,
          assay_name = assay_name
        )

    } else {

        ace = normalize.default(
          ace = ace,
          assay_name = assay_name,
          log_scale = TRUE
        )
        norm_method = "default"
    }

    metadata(ace)$normalization.method = norm_method
    
    return(ace)
}
