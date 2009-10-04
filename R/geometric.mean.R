"geometric.mean" <- 
function(x,na.rm=TRUE)
{ 
exp(mean(log(x),na.rm=na.rm)) }

