#factor analysis and sem diagrams
#based upon fa.graph with some ideas taken from the diagram and shape packages of Karline Soetaert 
#version of September 20, 2009
#developed to replace Rgraphviz which is too much of a pain to install
#Rgraphviz uses a NEL  (Node, Edge) represenation while diagram uses a complete linking matrix
#thus, I am trying to combine these two approaches
#modified 31/5/14 to allow for drawing factor extension derived omegas
#and again 7/2/18 to include omegaDirect output

"omega.diagram" <-
  function(om.results,sl=TRUE,sort=TRUE,labels=NULL,flabels=NULL,cut=.2,gcut=.2,simple=TRUE,errors=FALSE,
    digits=1,e.size=.1,rsize=.15,side=3,main=NULL,cex=NULL,color.lines=TRUE
    ,marg=c(.5,.5,1.5,.5),adj=2, ...) {
     if(color.lines) { colors <- c("black","red")} else {colors <- c("black","black") }
  Phi <- NULL  #the default case
  if(is.null(cex)) cex <- 1
   old.par<- par(mar=marg)  #give the window some narrower margins
    on.exit(par(old.par))  #set them back
 
 #Figure out what type of input we are drawing   (done July 2, 2018)
 #fixed for R 4.0.0  December 11, 2019
  extend<- FALSE
  
#which kind of input are we drawing?
    if(length(class(om.results)) > 1)  {
      omegaSem <- omegaDirect <- omega <- NULL  #strange fix to R 4,0,0 compiler
    names <- cs(matrix, extend, omegaSem, omegaDirect, omega)
    value <- inherits(om.results,names,which=TRUE)  # value <- class(x)[2]
    if(any(value > 1) ) { result <- names[which(value > 0)]} else {result <- "other"}
    
    } else {result <- "extend"}    #just a raw matrix
   if(result == "matrix") result <- "extend"

    
#dispatch to the right option
    switch(result,
extend = {extend <- TRUE
               #class(om.results)[2] <- "omega"
               factors <- om.results
               nvar <- num.var <- nrow(factors)
               num.factors <- ncol(factors) -1
                if(sort) {temp  <- fa.sort(factors[,-1])
                   temp2 <- factors[,1,drop=FALSE]     #added the drop option November 4, 2018
          factors <- cbind(g=temp2[rownames(temp),1],temp)}  #added column 1
            }, 
omegaSem = {
    #did we do an omegaSem or just an omegaFromSem?
    if(is.null(om.results$omega.efa$cfa.loads)) {cfa.loads <- om.results$cfa.loads} else {cfa.loads <- om.results$omega.efa$cfa.loads}
      #  factors <- as.matrix(om.results$omega.efa$cfa.loads[,2:ncol(om.results$omega.efa$cfa.loads)])
      class(cfa.loads) <- c("psych","omegaSem")
      if(sort) cfa.loads <- fa.sort(cfa.loads)
      factors <- as.matrix(cfa.loads)
        gloading <- cfa.loads[,1,drop=FALSE]
        nvar <- num.var <- nrow(gloading)
        num.factors <- ncol(factors) -1
        sl=TRUE
        main <- "Omega from SEM" 
    },
    
 omegaDirect = {factors <- om.results$loadings
     nvar <- num.var <- nrow(factors)
     num.factors <- ncol(factors) -1
                  if(sort) {temp  <- fa.sort(factors[,-1])
                   temp2 <- factors[,1]
          factors <- cbind(g=temp2[rownames(temp)],temp)} },
    
omega ={
        if(extend) class(om.results) <- c("psych","omega")
 if(sort) om.results <- fa.sort(om.results)   #usually sort, but sometimes it is better not to do so
          
 if (sl) {factors <- as.matrix(om.results$schmid$sl) 
         
         if(is.null(main)) {main <- "Omega with Schmid Leiman Transformation" }
         } else {factors <- as.matrix(om.results$schmid$oblique)
          if(is.null(main)) {main <- "Hierarchical (multilevel) Structure" }
         }
       gloading <- om.results$schmid$gloading
        nvar <- num.var <- dim(factors)[1]   #how many variables?
   if (sl ) {num.factors <- dim(factors)[2] - 1 - (!extend) *3 } else {num.factors <- dim(factors)[2]
       }
 },
 other ={warning("I am sorry, I don't know how to diagram this input")}
      
    )
if(result !="other") {#skip to the end if we don' know what we are doing
 
#now draw the figure
   
   e.size <- e.size * 10/ nvar   #this is an arbitrary setting that seems to work
#first some basic setup parameters 
  
   
   
    vars <- paste("V",1:num.var,sep="")  
   if (!is.null(labels)) {vars <- paste(labels)} else{vars <- rownames(factors) }
   if(!is.null(flabels)) {fact <- flabels} else { if(sl)  {fact <- c("g",paste("F",1:num.factors,"*",sep="")) } else {fact <- c(paste("F",1:num.factors,sep="")) }  } # e.g.  "g"  "F'1" "F2" "F3"
 
   colnames(factors)[1:length(fact)] <- fact
   var.rect <- list()
   arrows.list <- list()
   fact.rect <- list()
    max.len <- max(nchar(rownames(factors)))*rsize
   cex <-  min(cex,40/nvar)
   xleft <-  -max.len/2 #0#
   xright <- nvar + 1  # or hard code to 3?
  plot(0,type="n",xlim=c(xleft-max.len,xright+1),ylim=c(1,nvar+1),frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main)
   if (sl) {vloc <- (xright)/2
            gloc <- xleft
            grouploc <-xright
            start <- 1
            end <- num.factors+1} else {
            vloc <- xleft
            gloc <- xright
            grouploc <- (xright)/2
            start <- 0
            end <- num.factors
            }
  for (v in 1:nvar) { 
 	 var.rect[[v]] <- dia.rect(vloc,nvar-v+1,rownames(factors)[v],xlim=c(0,nvar),ylim=c(0,nvar),cex=cex,draw=FALSE,...)
     }
   f.scale <- (nvar+ 1)/(num.factors+1)
   f.shift <- nvar/num.factors
   for (f in 1:num.factors) {
   		fact.rect[[f]] <- dia.ellipse(grouploc,(num.factors+1-f)*f.scale,colnames(factors)[f+start],xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,draw=TRUE,...)
     		for (v in 1:nvar)  {
    			if (abs(factors[v,f+start]) > cut) {d.arrow <- dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,col=colors[((sign(factors[v,f+start])<0) +1)],lty=((sign(factors[v,f+start])<0)+1),labels=round(factors[v,f+start],digits),adj=f %% adj +1,draw=FALSE)
                   if(!is.null(d.arrow) )arrows.list <- c(arrows.list,d.arrow)            }
                             }
                     }
   
  g.ellipse <-  dia.ellipse(gloc,(num.var+1)/2,"g",xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)
   if(!sl) { 
   for (f in 1:num.factors) {
            d.arrow <-   dia.arrow(from=g.ellipse,to=fact.rect[[f]],col=colors[((sign(gloading[f])<0) +1)],lty=((sign(gloading[f])<0) +1),labels=round(gloading[f],digits),adj=f %% adj +1,draw=FALSE)                
            arrows.list <- c(arrows.list,d.arrow)                }
             } else {
              for (i in 1:nvar) {
              if(abs(factors[i,1]) > gcut) {
             d.arrow <- dia.arrow(from=g.ellipse,to=var.rect[[i]]$left,col=colors[((sign(factors[i,1])<0) +1)],lty=((sign(factors[i,1])<0)+1),labels=round(factors[i,1],digits),adj=1,draw=FALSE)}
               if(!is.null(d.arrow) ) arrows.list <- c(arrows.list,d.arrow)         
                         }
 
	  }
  if (errors) {for (v in 1:nvar) {
       dia.self(location=var.rect[[v]],scale=.5,side=side)  }
       }
    } #end of normal case
    
 tv <- matrix(unlist(var.rect),nrow=nvar,byrow=TRUE)
      all.rects.x <- tv[,5] #the center of the figure
      all.rects.y <- tv[,2]
      all.rects.names <- rownames(factors)
      dia.rect(all.rects.x, all.rects.y,all.rects.names) 
 #   multi.rect(var.rect,...)

#draw the factors  not implemented 
   fv <- matrix(unlist(fact.rect),nrow=num.factors,byrow=TRUE)
      all.rects.x <- fv[,5] #the center of the figure
      all.rects.y <- fv[,6]
      all.rects.names <- colnames(factors)[2:(num.factors+1)]
    #  dia.ellipse1(all.rects.x, all.rects.y,all.rects.names) 
      

 
#    now show all the loadings as arrows
#    tv <- matrix(unlist(arrows.list),byrow=TRUE,ncol=20)
#    cname<- colnames(tv)
#    tv  <- data.frame(tv)
#    tv[,c(1:18,20)] <- nchar2numeric(tv[,c(1:18,20)])
#    colnames(tv) <- cname
#         text(tv[,1],tv[,2],tv[,3],cex=tv[,4])
#         len1 <- as.numeric(tv[1,9])
#         len2 <- as.numeric(tv[1,16]) 
#         arrows(x0=tv[,5],y0=tv[,6],x1=tv[,7],y1=tv[,8],length=len1,angle=30,code=1,col=tv[,19],lty=tv[,20],...)
#         arrows(x0=tv[,12],y0=tv[,13],x1=tv[,14],y1=tv[,15],length=len2,angle=30,code=2,col=tv[,19],lty=tv[,29],...)
# 
   multi.arrow(arrows.list,...)   
}                  