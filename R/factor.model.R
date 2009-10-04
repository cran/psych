"factor.model" <-
function(f) { 
    result<- f %*% t(f)
    return (result)}

