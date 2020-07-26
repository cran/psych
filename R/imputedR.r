# A simple function to impute missing data correlations
#we find the mean r for each row and column and replace the missing value with geometric mean of the row and column for the missing value
#we can either find the colmean correlation or the median
#May 7, 2020

imputeR <- function(r,mean=TRUE,z=TRUE){
diag(r) <- NA
if(z){r <- fisherz(r)}
if(mean) {
rowR <- rowMeans(r,na.rm=TRUE)
colR <- colMeans(r,na.rm=TRUE)
} else {
rowR <- apply(r,1,median,na.rm=TRUE)
colR <- apply(r,2,median,na.rm=TRUE)}

if(z) {rowR <- fisherz2r(rowR)
       colR <- fisherz2r(colR)}
       
       
imputed <- (rowR %o% colR) #the outer product of the mean correlations
diag(imputed) <- NA
r [is.na(r)] <- sqrt(imputed[is.na(r)])
diag(r) <- 1
r [is.na(r) ] <- -sqrt(abs(imputed[is.na(r)]))
return(r)
}
