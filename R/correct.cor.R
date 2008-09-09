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
