
test.all <- function(p) {
 library(p,character.only=TRUE)
  ob <- paste("package",p,sep=":")
  ol <- objects(ob)
  nf <- length(ol)
  for(i in 1:nf) {
    fn <- as.character(ol[[i]])
    example(topic=fn,package=p,character.only=TRUE)
    }
 detach(ob,character.only=TRUE)
}


