"toDotty" <- 
function(graph,filename,...) {
if(!require(Rgraphviz)) {stop("You have called a function requiring Rgraphviz and it is not installed.  Install it and try again.")} else {toDot(graph,filename,...)}
}