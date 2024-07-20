HiddenVarDecompRaw <- function(Chain) {
  if (!is(Chain, "BASiCS_Chain")) {
    stop("'Chain' is not a BASiCS_Chain class object.")
  }

  N <- nrow(Chain@parameters$delta)
  q.bio <- ncol(Chain@parameters$delta)
  UniqueBatch <- colnames(Chain@parameters$theta)
  nBatch <- length(UniqueBatch)
  CellName <- colnames(Chain@parameters$s)

  if (nBatch > 1) {
    Theta <- matrixStats::rowMedians(Chain@parameters$theta)
  } else {
    Theta <- as.vector(Chain@parameters$theta)
  }

  if("phi" %in% names(Chain@parameters)) {
    PhiS <- matrixStats::rowMedians(Chain@parameters$phi * Chain@parameters$s)
  } else {
    PhiS <- matrixStats::rowMedians(Chain@parameters$s)
  }

  Aux <- (1 / (PhiS * Chain@parameters$mu)) + Chain@parameters$delta * (Theta + 1)
  
  RawTechVarGlobal <- Theta
  RawBioVarGlobal <- Chain@parameters$delta * (Theta + 1)
  RawTotalVarGlobal <- Aux + Theta

  TechVarGlobal <- Theta / (Aux + Theta)
  BioVarGlobal <- (Chain@parameters$delta * (Theta + 1)) / (Aux + Theta)

  TechVarBatch <- array(0, dim = c(N, q.bio, nBatch))
  BioVarBatch <- array(0, dim = c(N, q.bio, nBatch))
  RawTechVarBatch <- array(0, dim = c(N, q.bio, nBatch))
  RawBioVarBatch <- array(0, dim = c(N, q.bio, nBatch))
  RawTotalVarBatch <- array(0, dim = c(N, q.bio, nBatch))

  if (nBatch > 1) {
    for (Batch in seq_len(nBatch)) {
      SBatch <- Chain@parameters$s[, grep(UniqueBatch[Batch], CellName)]
      if("phi" %in% names(Chain@parameters)) {
        PhiBatch <- Chain@parameters$phi[, grep(UniqueBatch[Batch], CellName)]
        PhiSBatch <- matrixStats::rowMedians(PhiBatch * SBatch)
      } else {
        PhiSBatch <- matrixStats::rowMedians(SBatch)
      }

      AuxBatch <- (1 / (PhiSBatch * Chain@parameters$mu)) + Chain@parameters$delta * (Chain@parameters$theta[, Batch] + 1)
      
      RawTechVarBatch[, , Batch] <- Chain@parameters$theta[, Batch]
      RawBioVarBatch[, , Batch] <- Chain@parameters$delta * (Chain@parameters$theta[, Batch] + 1)
      RawTotalVarBatch[, , Batch] <- AuxBatch + Chain@parameters$theta[, Batch]

      TechVarBatch[, , Batch] <- Chain@parameters$theta[, Batch] / (AuxBatch + Chain@parameters$theta[, Batch])
      BioVarBatch[, , Batch] <- (Chain@parameters$delta * (Chain@parameters$theta[, Batch] + 1)) / (AuxBatch + Chain@parameters$theta[, Batch])
    }
  }

  if (nBatch > 1) {
    list(
      TechVarGlobal = TechVarGlobal,
      BioVarGlobal = BioVarGlobal,
      RawTechVarGlobal = RawTechVarGlobal,
      RawBioVarGlobal = RawBioVarGlobal,
      RawTotalVarGlobal = RawTotalVarGlobal,
      TechVarBatch = TechVarBatch,
      BioVarBatch = BioVarBatch,
      RawTechVarBatch = RawTechVarBatch,
      RawBioVarBatch = RawBioVarBatch,
      RawTotalVarBatch = RawTotalVarBatch
    )
  } else {
    list(
      TechVarGlobal = TechVarGlobal,
      BioVarGlobal = BioVarGlobal,
      RawTechVarGlobal = RawTechVarGlobal,
      RawBioVarGlobal = RawBioVarGlobal,
      RawTotalVarGlobal = RawTotalVarGlobal
    )
  }
}
