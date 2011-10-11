"cor2dist" <- function(x) {
   if(dim(x)[1] != dim(x)[2]) {x <- cor(x,use="pairwise")}
   dist <-  sqrt(2*(1-x))
   return(dist)
   }