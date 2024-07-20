#' @name BASiCS_VarianceDecompRaw
#' @aliases BASiCS_VarianceDecompRaw
#'
#' @title Decomposition of gene expression variability according to BASiCS
#'
#' @description Function to decompose total variability of gene
#' expression into biological and technical components.
#'
#' @param Chain an object of class \code{\linkS4class{BASiCS_Chain}}
#' @param OrderVariable Ordering variable for output.
#' Possible values: \code{'GeneName'}, \code{'BioVarGlobal'},
#'  \code{'TechVarGlobal'} and \code{'ShotNoiseGlobal'}.
#' Default: \code{OrderVariable = "BioVarGlobal"}.
#' @param Plot If \code{TRUE}, a barplot of the variance decomposition
#' (global and by batches, if any) is generated. Default: \code{Plot = TRUE}.
#' @param main Plot title.
#' @param ylab y axis label.
#' @param beside If \code{TRUE}, bars are placed beside each other.
#' If \code{FALSE}, bars are stacked.
#' @param palette Palette to be passed to \code{\link{scale_fill_brewer}}
#' to create a discrete colour mapping.
#' @param legend Labels for variance components.
#' @param names.arg X axis labels.
#'  
#'
#' @return A \code{\link[base]{data.frame}} whose first 4 columns correspond to
#' \describe{
#' \item{\code{GeneName}}{Gene name (as indicated by user)}
#' \item{\code{BioVarGlobal}}{Percentage of variance explained by a biological
#'                            component (overall across all cells)}
#' \item{\code{TechVarGlobal}}{Percentage of variance explained by the technical
#'                             component (overall across all cells)}
#' \item{\code{ShotNoiseGlobal}}{Percentage of variance explained by the shot
#'                               noise component (baseline Poisson noise,
#'                               overall across all cells)}
#' }
#' If more than 1 batch of cells are being analysed, the remaining columns
#' contain the corresponding variance decomposition calculated within each
#' batch.
#'
#' @examples
#'
#' # For illustration purposes we load a built-in 'BASiCS_Chain' object
#' # (obtained using the 'BASiCS_MCMC' function)
#' data(ChainSC)
#'
#' VD <- BASiCS_VarianceDecompRaw(ChainSC)
#'
#' @details See vignette
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references
#'
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology.
#'
#' @rdname BASiCS_VarianceDecompRaw
#' @export
BASiCS_VarianceDecompRaw <- function(Chain,
                                  OrderVariable = c(
                                    "BioVarGlobal", 
                                    "GeneName",
                                    "TechVarGlobal",
                                    "ShotNoiseGlobal"
                                  ),
                                  Plot = TRUE,
                                  main = "Overall variance decomposition",
                                  ylab = "% of variance",
                                  beside = FALSE,
                                  palette = "Set1",
                                  legend = c(
                                    "Biological",
                                    "Technical",
                                    "Shot noise"
                                  ),
                                  names.arg = if (nBatch == 1) "Overall"
                                    else c(
                                      "Overall",
                                      paste("Batch", seq_len(nBatch))
                                    )
                                  ) {
  
  OrderVariable <- match.arg(OrderVariable)
  if (!is(Chain, "BASiCS_Chain")) {
    stop("'Chain' is not a BASiCS_Chain class object.")
  }

  q.bio <- ncol(Chain@parameters$delta)
  UniqueBatch <- colnames(Chain@parameters$theta)
  nBatch <- length(UniqueBatch)

  # Calculating variance decomposition
  VarDecomp <- HiddenVarDecompRaw(Chain)

  # Global values
  BioVarGlobal <- matrixStats::colMedians(VarDecomp$BioVarGlobal)
  TechVarGlobal <- matrixStats::colMedians(VarDecomp$TechVarGlobal)
  ShotNoiseGlobal <- 1 - BioVarGlobal - TechVarGlobal

  RawBioVarGlobal <- matrixStats::colMedians(VarDecomp$RawBioVarGlobal)
  RawTechVarGlobal <- matrixStats::colMedians(VarDecomp$RawTechVarGlobal)
  RawTotalVarGlobal <- matrixStats::colMedians(VarDecomp$RawTotalVarGlobal)

  Genes <- seq_len(q.bio)
  GeneName <- colnames(Chain@parameters$mu)

  if (nBatch > 1) {
    VarDecompBatch <- NULL
    RawVarDecompBatch <- NULL

    for (Batch in seq_len(nBatch)) {
      BioVarAux <- matrixStats::colMedians(VarDecomp$BioVarBatch[, , Batch])
      TechVarAux <- matrixStats::colMedians(VarDecomp$TechVarBatch[, , Batch])
      RawBioVarAux <- matrixStats::colMedians(VarDecomp$RawBioVarBatch[, , Batch])
      RawTechVarAux <- matrixStats::colMedians(VarDecomp$RawTechVarBatch[, , Batch])
      RawTotalVarAux <- matrixStats::colMedians(VarDecomp$RawTotalVarBatch[, , Batch])

      VarDecompBatch <- cbind(
        VarDecompBatch,
        BioVarAux,
        TechVarAux,
        1 - BioVarAux - TechVarAux
      )

      RawVarDecompBatch <- cbind(
        RawVarDecompBatch,
        RawBioVarAux,
        RawTechVarAux,
        RawTotalVarAux
      )
    }

    colnames(VarDecompBatch) <- paste0(
      rep(c("BioVarBatch", "TechBatch", "ShotNoiseBatch"), nBatch),
      rep(seq_len(nBatch), each = 3)
    )

    colnames(RawVarDecompBatch) <- paste0(
      rep(c("RawBioVarBatch", "RawTechBatch", "RawTotalVarBatch"), nBatch),
      rep(seq_len(nBatch), each = 3)
    )

    out <- cbind.data.frame(
      GeneIndex = Genes,
      GeneName = GeneName,
      BioVarGlobal = BioVarGlobal,
      TechVarGlobal = TechVarGlobal,
      ShotNoiseGlobal = ShotNoiseGlobal,
      RawBioVarGlobal = RawBioVarGlobal,
      RawTechVarGlobal = RawTechVarGlobal,
      RawTotalVarGlobal = RawTotalVarGlobal,
      VarDecompBatch,
      RawVarDecompBatch,
      stringsAsFactors = FALSE
    )
  } else {
    out <- cbind.data.frame(
      GeneIndex = Genes,
      GeneName = GeneName,
      BioVarGlobal = BioVarGlobal,
      TechVarGlobal = TechVarGlobal,
      ShotNoiseGlobal = ShotNoiseGlobal,
      RawBioVarGlobal = RawBioVarGlobal,
      RawTechVarGlobal = RawTechVarGlobal,
      RawTotalVarGlobal = RawTotalVarGlobal,
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- Genes

  # Re-order before output
  orderVar <- switch(OrderVariable,
    "GeneName" = GeneName,
    "BioVarGlobal" = BioVarGlobal,
    "TechVarGlobal" = TechVarGlobal,
    "ShotNoiseGlobal" = ShotNoiseGlobal
  )

  out <- out[order(orderVar, decreasing = TRUE), ]

  if (Plot) {
    p <- BASiCS_PlotVarianceDecompRaw(
      out, xlabs = names.arg, beside = beside, main = main, ylab = ylab
    )
    print(p)
  }
  return(out)
}
