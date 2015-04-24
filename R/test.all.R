
test.all <- function(pl,package="psych",dependencies = c("Depends", "Imports", "LinkingTo"),find=FALSE,skip=NULL) {
 if (find) {
     pl <-tools::dependsOnPkgs(package,dependencies=dependencies) 
     if(!is.null(skip) && skip %in% pl) {pl <- pl[-which(skip ==pl)]}
     }
 np <- length(pl)
 if(np > 0 ) {
 for(i in 1:np) {
   p <- pl[i]
 test <- require(p,character.only=TRUE)
 if(!test) {cat("\nCould not find package ",p, "\n")
      next  
      }
 cat(paste("\nNow testing package " ,p ))
  ob <- paste("package",p,sep=":")
  ol <- objects(ob)
  nf <- length(ol)
  options("example.ask"=FALSE)
  for(i in 1:nf) {
    fn <- as.character(ol[[i]])
    example(topic=fn,package=p,character.only=TRUE)
    }
 detach(ob,character.only=TRUE)
} } else {cat("\nNo dependencies for package ", package) }
}


#tools::package_dependencies(reverse = TRUE)   #lists all the reverse dependencies
#tools::check_packages_in_dir(dir,reverse = list())   #might check them,  unclear
#rd <-reverse_dependencies_with_maintainers("psych")