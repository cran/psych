    "count.pairwise" <-
function (x, y=NULL) 
{
    sizex <- dim(x)[2]
    if (length(y)>0) 
        {sizey <- dim(y)[2]}  else  {sizey <- dim(x)[2]}
    result <- matrix(1, nrow = sizey, ncol = sizex)
    xnames <- names(x)
    colnames(result) <- names(x)
    if (((is.data.frame(y)) | (is.matrix(y)))) { 
        rownames(result) <- names(y) }  else {rownames(result) <- names(x)}
    if (length(y) ==0) {
        for (i in 1:sizex) {
            for (j in 1:(i)) {
                result[j, i] <- sum((!is.na(x[, j])) & (!is.na(x[, 
                  i])))
                result[i, j] <- result[j, i]
            }
        }
    } else {
        for (i in 1:sizex) {
            for (j in 1:sizey) {
                result[j, i] <- sum((!is.na(x[, i])) & (!is.na(y[, 
                  j])))
            }
        }
    }
    return(result) }
    
    
    
