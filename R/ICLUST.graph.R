#created August 4, 2006
#last revised September 5, 2006
#create dot file for graphviz from ICLUST output
#inspired by sem.path.diagram by J. Fox

"ICLUST.graph"  <- 
function(ic.results, out.file,min.size=1,short=FALSE,labels=NULL,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 12),  rank.direction=c("RL","TB","LR","BT"), digits=2,title="ICLUST", ...){
    
      if(!missing(out.file)){
        out <- file(out.file, "w")
        on.exit(close(out))
        }
        else out <- stdout()
        
    results <- ic.results$results
    if (length(labels)==0) {
    var.labels <- rownames(ic.results$loadings)} else {var.labels=labels}
   # if(dim(var.labels)[1] < dim(var.labels)[2]) {var.labels <- t(var.labels)}
    clusters <- as.matrix(ic.results$clusters)
   #if(length(clusters)==length(var.labels) )  {clusters <- as.matrix(clusters)}
	num <- nrow(results)
    if (short) {var.labels <- paste("V",1:nrow,(var.labels),sep="")}   
   
  rank.direction <- match.arg(rank.direction)
  #first some basic setup parameters 
  cat( file=out,paste('digraph ICLUST', ' {\n', sep=""))
  cat(file=out, paste('  rankdir=', rank.direction, ';\n', sep=""))
  cat(file=out, paste('  size="',size[1],',',size[2],'";\n', sep=""))
  cat(file=out, paste('  node [fontname="', node.font[1], 
        '" fontsize=', node.font[2], ' shape=box, width=2];\n', sep=""))
  cat(file=out, paste('  edge [fontname="', edge.font[1],
        '" fontsize=', edge.font[2], '];\n', sep=""))
  cat(file=out, paste(' label = "' ,title,'";
	fontsize=20;\n', sep=""))
 
   #create the items as boxes  
   #add the sign from the clusters 
  num.var <- nrow(results)+1   #how many variables?
  if (num.var > dim(clusters)[1]) {num.var <-  dim(clusters)[1]}
     for (i in 1:num.var) { if (max(clusters[i,]) > 0 ) {
      	cat(file=out,paste('V',i,'  [label = "',var.labels[i], '"];\n', sep="")) } else {
      	cat(file=out,paste('V',i,'  [label = "-',var.labels[i], '"];\n', sep="")) }
     	 }
   
  #show the cluster structure with ellipses
  
  cat(file=out,paste('node [shape=ellipse, width ="1"];\n', sep=""))
  #draw the edges
  for (i in 1:num) {if(results[i,1]>0) { #avoid printing null results
     cat(file=out,paste(row.names(results)[i],  '-> ', results[i,1], ' [ label = ',round(results[i,"r1"],digits),' ];\n', sep=""))
     cat(file=out,paste(row.names(results)[i],  '-> ', results[i,2], ' [ label = ',round(results[i,"r2"],digits),' ];\n', sep=""))
     }}
    
   #label the clusters with alpha and beta
   for (i in 1:num) {if(results[i,1]>0) { #don't print blank results
   if (results[i,"size"] > min.size) {
     cat(file=out,paste(row.names(results)[i],  '  [label =   "',row.names(results)[i],'\\n  alpha= ',round(results[i,"alpha"],digits),'\\n beta=  ' ,round(results[i,"beta"],digits),'\\nN= ',results[i,"size"], '"] ;\n', sep=""))
     } else {cat(file=out,paste(row.names(results)[i],' ;\n', sep="")) } #short names for small clusters
     }}
  
  #keep the boxes all at the same rank (presumably the left side)
  cat(file=out, paste('{ rank=same;\n', sep=""))
  for (i in 1:num.var) { cat(file=out,paste('V',i,';', sep=""))
  }   
   cat(file=out, paste('}}', sep=""))   # we are finished
 } # end of ICLUST.graph 
