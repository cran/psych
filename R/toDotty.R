"toDotty" <- 
function(graph,filename,...) {
if(!require(Rgraphviz) && !require(graph)) {stop("You have called a function requiring Rgraphviz and it is not installed.  Install it and try again.")} else {
require(Rgraphviz)
require(graph)
toDot(graph,filename,...)}
}