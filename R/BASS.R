BASS<-function(R, maxP=5, Print = "ON"){
#---------------------------------------------------------

# Program to compute Goldberg's Bass Ackwards Procedure
# from a correlation matrix (R). PC with Varimax Rotation
# Niels Waller, May 10, 2006
#
# Program arguments:
# R = input correlation matrix
# maxP = maximum number of components to rotate
# Print = ON/OFF to print summarzed findings to screen
#----------------------------------------------------------
varNames <- rownames(R, do.NULL = FALSE, prefix = "var")

ULU <- eigen(R)
U <- ULU$vectors
L <- ULU$values
key <- sign(apply(U, 2, sum))
key[key < 0] <- -1
U <- U %*% diag(key)
P <- U %*% diag(sqrt(L))
p <- ncol(R)
CrossLevelCors <- list(rep(0, p - 1))
T <- list(rep(0, p - 1))
PCloadings <- list(rep(0, p - 1))
for (i in 2:maxP) {
vout <- varimax(P[, 1:i], normalize = TRUE, eps = 1e-15)

T[[i - 1]] <- vout$rotmat
PCloadings[[i - 1]] <- vout$loadings[1:p, ]
rownames(PCloadings[[i - 1]]) <- varNames
}
Z <- paste("Z", 1:maxP, sep = "")
V <- paste("V", 1:maxP, sep = "")
if (Print == "ON") {
cat("nCorrelation of", Z[1], " with ", V[2], "n")

}
out <- T[[1]][1, ]
dim(out) <- c(1, 2)
rownames(out) <- Z[1]
colnames(out) <- paste(V[2], ".", 1:2, sep = "")
CrossLevelCors[[1]] <- out
if (Print == "ON") {
print(round(out, 3))

}
for (i in 2:(maxP - 1)) {
if (Print == "ON") {
cat("nnnCorrelation of", V[i], " with ", V[i + 1], "nn")

}
S <- cbind(diag(i), matrix(0, i, 1))
out <- t(T[[i - 1]]) %*% S %*% T[[i]]
rownames(out) <- paste(V[i], ".", 1:i, sep = "")
colnames(out) <- paste(V[i + 1], ".", 1:(i + 1), sep = "")
CrossLevelCors[[i]] <- out
if (Print == "ON") {
print(round(out, 3))

}
}
invisible(list(T = T, cors = CrossLevelCors, loadings = PCloadings))
}