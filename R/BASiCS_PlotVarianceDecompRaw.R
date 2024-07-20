#' Plot variance decomposition results.
#' @param Decomp The output of \code{\link{BASiCS_VarianceDecomp}}.
#' @param beside If \code{TRUE}, bars are placed beside each other.
#' If \code{FALSE}, bars are stacked.
#' @param nBatch Number of batches.
#' @param main Plot title.
#' @param xlabs x-axis labels. Defaults to "Batch 1", "Batch 2", etc.
#' @param ylab y axis label.
#' @return A ggplot object.
BASiCS_PlotVarianceDecompRaw <- function(
    Decomp,
    beside = FALSE,
    nBatch = ((ncol(Decomp) - 2) / 3) - 1,
    main = "Overall variance decomposition",
    xlabs = if (nBatch == 1) "Overall"
      else c(
        "Overall",
        paste("Batch", seq_len(nBatch))
      ),
    ylab = "Raw variance"
  ) {
  outmat <- matrix(
    apply(Decomp[, grep("Raw", colnames(Decomp))], 2, mean),
    nrow = 3, byrow = FALSE
  )
  rownames(outmat) <- c("Raw Total", "Raw Technical", "Raw Biological")
  colnames(outmat) <- xlabs
  mdf <- reshape2::melt(outmat)
  ggplot(mdf,
      aes(x = .data$Var2, y = .data$value, fill = .data$Var1)
    ) +
    geom_col(position = if (!beside) "dodge" else "stack") +
    scale_fill_brewer(palette = "Set1", name = NULL) +
    labs(title = main, x = NULL, y = ylab)
}
