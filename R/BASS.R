BASS <- function(R, maxP = 5, Print = TRUE) {
  # Program to compute Goldberg's Bass Ackwards Procedure
  # from a correlation matrix (R). PC with Varimax Rotation
  # Niels Waller, May 10, 2006
  #
  # Arguments:
  # R: Input correlation matrix
  # maxP: Maximum number of components to rotate
  # Print: Boolean flag to print summarized findings to screen
  
  varNames <- rownames(R)
  ULU <- eigen(R)
  U <- ULU$vectors
  L <- ULU$values
  
  key <- sign(apply(U, 2, sum))
  key[key < 0] <- -1
  U <- U %*% diag(key)
  P <- U %*% diag(sqrt(L))
  
  p <- ncol(R)
  CrossLevelCors <- vector("list", maxP - 1)
  T <- vector("list", maxP - 1)
  PCloadings <- vector("list", maxP - 1)
  
  for (i in 2:maxP) {
    vout <- varimax(P[, 1:i], normalize = TRUE, eps = 1e-15)
    T[[i - 1]] <- vout$rotmat
    PCloadings[[i - 1]] <- vout$loadings[1:p, ]
    rownames(PCloadings[[i - 1]]) <- varNames
  }
  
  Z <- paste("Z", 1:maxP, sep = "")
  V <- paste("V", 1:maxP, sep = "")
  
  if (Print) {
    for (i in 1:(maxP - 1)) {
      cat("\nCorrelation of", Z[i], "with", V[i + 1], "\n\n")
      print(round(T[[i]][1, ], 3))
    }
  }
  
  for (i in 1:(maxP - 2)) {
    S <- cbind(diag(i), matrix(0, i, 1))
    out <- t(T[[i]]) %*% S %*% T[[i + 1]]
    rownames(out) <- paste(V[i], ".", 1:i, sep = "")
    colnames(out) <- paste(V[i + 1], ".", 1:(i + 1), sep = "")
    CrossLevelCors[[i]] <- out
  }
  
  list(T = T, cors = CrossLevelCors, loadings = PCloadings)
}
