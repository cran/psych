"geometric.mean" <- 
function(x,na.rm=TRUE){ if (is.null(nrow(x))) {exp(mean(log(x),na.rm=na.rm)) } else {
exp(apply(log(x),2,mean,na.rm=na.rm))} }

