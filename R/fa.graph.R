#take the output from a factor analysis and create the dot code
#Created April 20, 2012 to allow making dot code without using rgraphviz

"fa.graph" <-
function(fa.results,out.file=NULL,labels=NULL,cut=.3,simple=TRUE,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction=c("RL","TB","LR","BT"), digits=1,main="Factor Analysis", ...){
    	 Phi <- NULL  #the default case
  
   if(!missing(out.file)){
        out <- file(out.file, "w")  #set it to write permission
        on.exit(close(out))
        }
        else out <- stdout()   #use the normal console for output
        
 if((!is.matrix(fa.results)) && (!is.data.frame(fa.results)))  {factors <- as.matrix(fa.results$loadings)
                if(!is.null(fa.results$Phi)) Phi <- fa.results$Phi} else {factors <- as.matrix(fa.results)}
   rank.direction <- match.arg(rank.direction)
   
  #first some basic setup parameters 
   if (length(labels)==0) { labels <- rownames(fa.results$loadings)} else {labels=labels}
   num.var <- dim(factors)[1]   #how many variables?
   if (is.null(num.var) ){num.var <- length(factors)
       num.factors <- 1} else {
   num.factors <- dim(factors)[2]}
   vars <- paste("V",1:num.var,sep="")   
   fact <- paste("F",1:num.factors,sep="")
   
    #first some basic setup parameters 
  cat( file=out,paste('digraph Factor ', ' {\n', sep=""))
  cat(file=out, paste('  rankdir=', rank.direction, ';\n', sep=""))
  cat(file=out, paste('  size="',size[1],',',size[2],'";\n', sep=""))
  cat(file=out, paste('  node [fontname="', node.font[1], 
        '" fontsize=', node.font[2], ' shape=box, width=2];\n', sep=""))
  cat(file=out, paste('  edge [fontname="', edge.font[1],
        '" fontsize=', edge.font[2], '];\n', sep=""))
  #cat(file=out, paste(' label = "' ,main,'";fontsize=20;\n', sep=""))

  if(simple) {  # do simple structure if requested     	
 rowmax <- apply(abs(factors),1,max)
 factors[abs(factors) < rowmax] <- 0
     }
 
     #create the items as boxes 
 for (i in 1:num.var) {
      	cat(file=out,paste('V',i,'  [label = "',labels[i], '"];\n', sep="")) 
      	} 
 
  #show the factors as ellipses
   cat(file=out,paste('node [shape=ellipse, width ="1"];\n', sep=""))
     
     #draw the loadings    
for( nf in 1:num.factors) { 
  for (i in 1:num.var) {if(abs(factors[i,nf]) > cut ) { #avoid printing null results
    cat(file=out,paste(colnames(factors)[nf],  '-> V', i, ' [ label = ',round(factors[i,nf],digits),' ];\n', sep=""))
     }
     }
    }
 
#draw the interfactor correlations
 if(!is.null(Phi)) {
for (f in 2:num.factors)  {
     for (f1 in 1:(f-1)) { if(abs(Phi[f,f1]) > cut) {
        cat(file=out,paste(colnames(factors)[f],  ' -> ', colnames(factors)[f1], ' [ label = ',round(Phi[f,f1],digits),' , dir="both" ];\n', sep=""))
                         }
                               } 
 }   
}

  #keep the boxes all at the same rank (presumably the left side)
  cat(file=out, paste('{ rank=same;\n', sep=""))
  for (i in 1:num.var) { cat(file=out,paste('V',i,';', sep="")) } 
  cat(file=out, paste('}', sep="")) 
  cat(file=out, paste('{ rank=same;\n', sep=""))
  for (nf in 1:num.factors) { cat(file=out,paste(paste(colnames(factors)[nf],';', sep="")))
  }   
   cat(file=out, paste('}}', sep=""))   # we are finished
   }
 
