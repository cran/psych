"toDotty" <- 
function(graph,filename,...) {
if(!requireNamespace('Rgraphviz') && !requireNamespace(graph)) {stop("You have called a function requiring Rgraphviz and it is not installed.  Install it and try again.")
toDot <- function() {}   #dummy function} else {
Rgraphviz::toDot(graph,filename,...)}
}