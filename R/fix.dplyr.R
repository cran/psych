"fix.dplyr" <- function (object) {
   if(is.data.frame(object)) {   
   cn <- class(object)
df <- which(cn=="data.frame")
cn.not <- cn[-df]
cn <- c("data.frame",cn.not)
class(object) <- cn
      } 
invisible(object) 
}