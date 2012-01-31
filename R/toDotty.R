"toDotty" <- 
function(graph,filename,...) {
if(!require(Rgraphviz) && !require(graph)) {stop("You have called a function requiring Rgraphviz and it is not installed.  Install it and try again.")} else {toDot(graph,filename,...)}
}