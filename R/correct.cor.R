"correct.cor" <-
function(x,y) { n=dim(x)[1]   
        { diag(x) <- y 
        if (n> 1)  {
        for (i in 2:n) {
           k=i-1
           for (j in 1:k) {
              x[j,i] <- x[j,i]/sqrt(y[i]*y[j])  }   #fix the upper triangular part of the matrix
             }}
           return(x)  }}
           
"rangeCorrection" <- function(r,sdu,sdr) {
return(  r * sdu/(sdr* sqrt(1-r^2 + r^2*(sdu^2/sdr^2))))
}
