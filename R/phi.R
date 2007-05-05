# slight changes to combine phi and phi1 from W. Revelle
# Leo GŸrtler 07-09-06
# from http://www.personality-project.org/r/useful.r
"phi" <-
 function(t,digits=2)
{  # expects: t is a 2 x 2 matrix or a vector of length(4)
   stopifnot(prod(dim(t)) == 4 || length(t) == 4)
   if(is.vector(t)) t <- matrix(t, 2)
   r.sum <- rowSums(t)
   c.sum <- colSums(t)
   total <- sum(r.sum)
   r.sum <- r.sum/total
   c.sum <- c.sum/total
   v <- prod(r.sum, c.sum)
   phi <- (t[1,1]/total - c.sum[1]*r.sum[1]) /sqrt(v)
return(round(phi,2))  }