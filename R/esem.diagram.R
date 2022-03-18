#adapted from my structure.diagram function to draw the output from esem
#August 4, 2016
"esem.diagram" <-
function(esem=NULL,labels=NULL,cut=.3,errors=FALSE,simple=TRUE,regression=FALSE,lr=TRUE,
  digits=1,e.size=.1,adj=2,main="Exploratory Structural Model", ...){

#a helper function
 sort.f <- function(x) {
   nvar <- ncol(x)
   if(is.null(nvar)) {return(x)} else {
   nitem <- nrow(x)
   cluster <- data.frame(item <- seq(1:nitem),clust = rep(0,nitem))
   cluster$clust <- apply(abs(x),1,which.max)
   ord <- sort(cluster$clust,index.return=TRUE)
   x[1:nitem,] <- x[ord$ix,]
   rownames(x) <- rownames(x)[ord$ix]
   return(x)}}


#first some default values

 xmodel <- sort.f(esem$loadsX)
 ymodel <- sort.f(esem$loadsY)
 Phi       <-  esem$Phi
 num.y  <- num.x <- 0   #we assume there is nothing there
 vars<- NULL

 

 num.xvar <- dim(xmodel)[1]   #how many x variables?
 num.yvar <- dim(ymodel)[1]
 num.xfactors <- dim(xmodel)[2]
 num.yfactors <- dim(ymodel)[2]
 if(is.null(num.xvar)) num.xvar <- length(xmodel)
 if(is.null(num.yvar)) num.yvar <- length(ymodel)
 if(is.null(num.xfactors)) num.xfactors <- 1
 if(is.null(num.yfactors)) num.yfactors <- 1
 
# if(max(num.xvar,num.yvar) < 10) e.size <- e.size* 2  #make the ellipses bigger for small problems
  e.size <- e.size * 10/ max(num.xvar,num.yvar)	
  if(is.null(labels)) { xvars <-  rownames(xmodel)} else { xvars <-  vars <- labels}
  if(is.null(ncol(xmodel))) xvars <- names(xmodel)
  
#  if(is.null(vars) ) {xvars <- paste0("x",1:num.xvar)  }
  fact <- colnames(xmodel)
  if (is.null(fact)) { fact <- paste0("X",1:num.xfactors) }
  
  
  
  
     if(is.null(ncol(ymodel))) {yvars <-names(ymodel) } else {yvars <- rownames(ymodel)}
     if(is.null(yvars)) {yvars <- paste0("y",1:num.y)  }
     if(is.null(labels)) {vars <- c(xvars,yvars)} else {yvars <- labels[(num.xvar+1):(num.xvar+num.y)]} 
     
      yfact <- colnames(ymodel)
      if(is.null(yfact)) {yfact <- paste0("Y",1:num.yfactors) }
      fact <- c(fact,yfact)
     
                        
   
   
    num.var <- num.xvar + num.y
    num.factors <- num.xfactors + num.yfactors
    
     sem <- matrix(rep(NA),6*(num.var*num.factors + num.factors),ncol=3)    #this creates an output model for sem analysis
     lavaan <- vector("list",num.xfactors + num.yfactors) #create a list for lavaan
     colnames(sem) <- c("Path","Parameter","Value")
    var.rect <- list()
     fact.rect <- list()
#     
    
                           


length.labels <-  0    # a filler for now
#plot.new() is necessary if we have not plotted before
#strwd <- try(strwidth(xvars),silent=TRUE)
strwd <- try(length.labels <- max(strwidth(xvars),strwidth(yvars),strwidth("abc"))/1.8,silent=TRUE)  #although this throws an error if the window is not already open, we don't show it
#if (class(strwd) == "try-error" ) {plot.new() } 
if (inherits(strwd, "try-error" )) {length.labels = max(nchar(xvars),3)/1.8 } 
#length.labels <- max(strwidth(xvars),strwidth("abc"))/1.8

if(lr) {limx <- c(-(length.labels),max(num.xvar,num.yvar)+3 )  
        limy <-  c(0,max(num.xvar,num.yvar)+1) } else {
        limy <- c(-(length.labels),max(num.xvar,num.yvar) +3 )  
        limx <-  c(0,max(num.xvar,num.yvar)+1)
       if( errors) limy <- c(-1,max(num.xvar,num.yvar)+2)}
scale.xaxis <- 3


if(lr) {plot(0,type="n",xlim=limx,ylim=limy,frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main)} else  {plot(0,type="n",xlim=limx,ylim=limy,frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main) }
    
   #now draw the x part

   k <- num.factors  
   x.adjust <- num.xvar/ max(num.xvar,num.yvar)

    for (v in 1:num.xvar) {   
 	if(lr) { var.rect[[v]] <- dia.rect(0,(num.xvar-v+1)/x.adjust,xvars[v],xlim=limx,ylim=limy,...) } else { var.rect[[v]] <- dia.rect(v,0,xvars[v],xlim=limy,ylim=limx,...) }
     }
  nvar <- num.xvar
  f.scale <- (num.xvar+ 1)/(num.xfactors+1)/x.adjust
  
  factors <- round(xmodel,digits)
  if(is.null(ncol(factors))) factors <- matrix(factors)
  if (num.xfactors >0) { 
 for (f in 1:num.xfactors) {
   	if(!regression) {if(lr) {fact.rect[[f]] <- dia.ellipse(limx[2]/scale.xaxis,(num.xfactors+1-f)*f.scale,fact[f],xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)} else {fact.rect[[f]] <- dia.ellipse(f*f.scale,limy[2]/scale.xaxis,fact[f],ylim=c(0,nvar),xlim=c(0,nvar),e.size=e.size,...)
   	}
   	} else {if(lr) {fact.rect[[f]] <- dia.rect(limx[2]/scale.xaxis,(num.xfactors+1-f)*f.scale,fact[f],xlim=c(0,nvar),ylim=c(0,nvar),...)} else {
   	                fact.rect[[f]] <- dia.rect(f*f.scale,limy[2]/scale.xaxis,fact[f],xlim=c(0,nvar),ylim=c(0,nvar),...)}
   	    }
   	
     		for (v in 1:num.xvar)  {
     		    if(is.numeric(factors[v,f])) {
     		    if(simple && (abs(factors[v,f]) == max(abs(factors[v,])) )  && (abs(factors[v,f]) > cut) | (!simple && (abs(factors[v,f]) > cut))) { if (!regression) {if(lr){dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,labels =factors[v,f],col=((sign(factors[v,f])<0) +1),lty=((sign(factors[v,f])<0) +1),adj=f %% adj +1) 
    			} else {dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$top,labels =factors[v,f],col=((sign(factors[v,f])<0) +1),lty=((sign(factors[v,f])<0) +1)) 
    			}
    		} else {dia.arrow(to=fact.rect[[f]]$left,from=var.rect[[v]]$right,labels =factors[v,f],col=((sign(factors[v,f])<0) +1))} }
                                    }  else {
                if (factors[v,f] !="0") {
               if (!regression) { if(lr) {dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,labels =factors[v,f]) } else {dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$top,labels =factors[v,f])} 
    			} else {if(lr) {dia.arrow(to=fact.rect[[f]],from=var.rect[[v]]$right,labels =factors[v,f])} else {dia.arrow(to=fact.rect[[f]],from=var.rect[[v]]$top,labels =factors[v,f])}
    			   } }
                                    } }
                              } 
      if (num.xfactors ==1) { 
                   lavaan[[1]] <- paste(fact[1],"=~ ")
                      for(i in 1:num.xvar) {
                       sem[i,1] <- paste(fact[1],"->",vars[i],sep="")
                         lavaan[[1]] <-  paste0(lavaan[[1]], ' + ', vars[i]) 
                          if(is.numeric(factors[i])) {sem[i,2] <- vars[i]} else {sem[i,2] <- factors[i] }
                        }}  #end of if num.xfactors ==1 
       k <- num.xvar+1 
  
                   k <- 1

                   for (f in 1:num.xfactors) { #if (!is.numeric(factors[i,f]) ||  (abs(factors[i,f]) > cut))
                    lavaan[[f]] <- paste0(fact[f] ," =~ ")
                     for (i in 1:num.xvar) {
                   
                   if((!is.numeric(factors[i,f] ) && (factors[i,f] !="0"))||  ((is.numeric(factors[i,f]) && abs(factors[i,f]) > cut ))) {
                              
                              sem[k,1] <- paste(fact[f],"->",vars[i],sep="")
                              lavaan[[f]] <-  paste0(lavaan[[f]], ' + ', vars[i]) 
                             if(is.numeric(factors[i,f])) {sem[k,2] <- paste("F",f,vars[i],sep="")} else {sem[k,2] <- factors[i,f]}
                              k <- k+1 }   #end of if 
                         }  
                        }      
   }  #end of if num.xfactors >0     
  if(errors) {  for (i in 1:num.xvar) {if(lr) { dia.self(var.rect[[i]],side=2) } else { dia.self(var.rect[[i]],side=1)}
                                          sem[k,1] <- paste(vars[i],"<->",vars[i],sep="")
                                          sem[k,2] <- paste("x",i,"e",sep="")
                                    k <- k+1 }
                    }  
     
    

 #do it for y model 
  if(!is.null(ymodel)) { 
  	if(lr) {   y.adj <- num.yvar/ max(num.xvar,num.yvar)
  	#y.adj <-  num.yvar/2 - num.xvar/2 
  	        f.yscale <- limy[2]/(num.yfactors+1)
   	        y.fadj <- 0} else {
   	       #  y.adj <-  num.xvar/2 - num.yvar/2
   	         y.adj <- num.yvar/2
   	         f.yscale <- limx[2]/(num.yfactors+1)
   	          y.fadj <- 0} 
   	 
        for (v in 1:num.yvar) { if(lr){   var.rect[[v+num.xvar]] <- dia.rect(limx[2]-.35,limy[2]-v /y.adj,yvars[v],xlim=limx,ylim=limy,...)} else {
                                     var.rect[[v+num.xvar]] <- dia.rect(v / y.adj,limx[2],yvars[v],xlim=limy,ylim=limx,...)}
                                  }
                          
                  }
   #we have drawn the y variables, now should we draw the Y factors
   if(!is.null(ymodel)){
   y.factors <- round(ymodel,digits)
   if(is.null(ncol(y.factors))) { y.factors <- matrix(y.factors) }
     
   num.y <- nrow(y.factors)
  
   	for (f in 1:num.yfactors) {if(lr) {
   		fact.rect[[f+num.xfactors]] <- dia.ellipse(2*limx[2]/scale.xaxis,(num.yfactors+1-f)*f.yscale +y.fadj,yfact[f],xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)} else {
   		fact.rect[[f+num.xfactors]] <- dia.ellipse(f*f.yscale+ y.fadj,2*limx[2]/scale.xaxis,yfact[f],ylim=c(0,nvar),xlim=c(0,nvar),e.size=e.size,...)}
   		  
     		for (v in 1:num.yvar) {if(is.numeric(y.factors[v,f])) {
     		{if(simple && (abs(y.factors[v,f]) == max(abs(y.factors[v,])) )  && (abs(y.factors[v,f]) > cut) | (!simple && (abs(factors[v,f]) > cut))) {
     		     if(lr) { dia.arrow(from=fact.rect[[f+num.xfactors]],to=var.rect[[v+num.xvar]]$left,labels =y.factors[v,f],col=((sign(y.factors[v,f])<0) +1),adj = f %% adj +1)} else {
     		               dia.arrow(from=fact.rect[[f+num.xfactors]],to=var.rect[[v+num.xvar]]$bottom,labels =y.factors[v,f],col=((sign(y.factors[v,f])<0) +1))}
     		     }
                                    }
                   } else {if(factors[v,f] !="0")  {if(lr) {dia.arrow(from=fact.rect[[f+num.xfactors]],to=var.rect[[v+num.xvar]]$left,labels =y.factors[v,f]) } else {
                                                            dia.arrow(from=fact.rect[[f+num.xfactors]],to=var.rect[[v+num.xvar]]$bottom,labels =y.factors[v,f])
                   }
                   }
             }}
                             } 
                       
  
  	if (num.yfactors ==1) {  lavaan[[num.xfactors +1 ]] <-   paste(fact[num.xfactors +1], "=~")      
    for (i in 1:num.y) { sem[k,1] <- paste(fact[1+num.xfactors],"->",yvars[i],sep="")
                      lavaan[[num.xfactors +1]] <-  paste0(lavaan[[num.xfactors +1]], ' + ', yvars[i]) 
                           if(is.numeric(y.factors[i] ) ) {sem[k,2] <- paste("Fy",yvars[i],sep="")} else {sem[k,2] <- y.factors[i]}
                         k <- k +1
                        }                      
                         } else {   #end of if num.yfactors ==1 
        
                   for (i in 1:num.y) {
                   for (f in 1:num.yfactors) { lavaan[[num.xfactors +f ]] <-   paste(fact[num.xfactors +f], "=~")  
                    if( abs(y.factors[i,f]) > cut ) {
                           sem[k,1] <- paste(fact[f+num.xfactors],"->",vars[i+num.xvar],sep="")
                            lavaan[[num.xfactors +f]] <-  paste0(lavaan[[num.xfactors +f]], ' + ', yvars[i]) 
                           if(is.numeric(y.factors[i,f])) { sem[k,2] <- paste("Fy",f,vars[i+num.xvar],sep="")} else {sem[k,2] <- y.factors[i,f]}
                         k <- k+1 }   #end of if 
                         }  #end of factor
                        } # end of variable loop      
                               } #end of else if 
                  # }
       
     if(errors) {  for (i in 1:num.y) { 
                      if(lr) {dia.self(var.rect[[i+num.xvar]],side=3) } else {dia.self(var.rect[[i+num.xvar]],side=3)}
                    sem[k,1] <- paste(vars[i+num.xvar],"<->",vars[i+num.xvar],sep="")
                    sem[k,2] <- paste("y",i,"e",sep="")
                    k <- k+1 }}
  
    } #end of if.null(ymodel)

# if(!is.null(Rx)) {#draw the correlations between the x variables
#                 for (i in 2:num.xvar) {
#                      	for (j in 1:(i-1)) {  
#                      	 if((!is.numeric(Rx[i,j] ) && ((Rx[i,j] !="0")||(Rx[j,i] !="0")))||  ((is.numeric(Rx[i,j]) && abs(Rx[i,j]) > cut ))) {
#                           if (lr) {if(abs(i-j) < 2) { dia.curve(from=var.rect[[j]]$left,to=var.rect[[i]]$left, labels = Rx[i,j],scale=-3*(i-j)/num.xvar)} else {        dia.curve(from=var.rect[[j]]$left,to=var.rect[[i]]$left, labels = Rx[i,j],scale=-3*(i-j)/num.xvar)} 
#                         	 } else {
#                         	 if(abs(i-j) < 2) { dia.curve(from=var.rect[[j]]$bottom,to=var.rect[[i]]$bottom, labels = Rx[i,j],scale=-3*(i-j)/num.xvar)} else {dia.curve(from=var.rect[[j]]$bottom,to=var.rect[[i]]$bottom, labels = Rx[i,j],scale=-3*(i-j)/num.xvar)}
#                          }	
#                            }}
#                            }
#                
#                  }
# if(!is.null(Ry)) {#draw the correlations between the y variables
#                 for (i in 2:num.yvar) {
#                      	for (j in 1:(i-1)) {  
#                      	 if((!is.numeric(Ry[i,j] ) && ((Ry[i,j] !="0")||(Ry[j,i] !="0")))||  ((is.numeric(Ry[i,j]) && abs(Ry[i,j]) > cut ))) {
#                           if (lr) {if(abs(i-j) < 2) { dia.curve(from=var.rect[[j+num.xvar]]$right,to=var.rect[[i+num.xvar]]$right, labels = Ry[i,j],scale=3*(i-j)/num.xvar)} else {dia.curve(from=var.rect[[j+num.xvar]]$right,to=var.rect[[i+num.xvar]]$right, labels = Ry[i,j],scale=3*(i-j)/num.xvar)} 
#                         	 } else {
#                         	 if(abs(i-j) < 2) { dia.curve(from=var.rect[[j+num.xvar]]$bottom,to=var.rect[[i+num.xvar]]$bottom, labels = Ry[i,j],scale=3*(i-j)/num.xvar)} else {dia.curve(from=var.rect[[j+num.xvar]]$bottom,to=var.rect[[i+num.xvar]]$bottom, labels = Ry[i,j],scale=3*(i-j)/num.xvar)}
#                          }	
#                            }}
#                            }
#                
#                  }
        
        Phi <- round(Phi,digits)
        
if(!regression) {
 if(!is.null(Phi)) {if (!is.matrix(Phi)) { if(!is.null(FALSE)) {Phi <- matrix(c(1,0,Phi,1),ncol=2)} else {Phi <- matrix(c(1,Phi,Phi,1),ncol=2)}} 
 
 if(num.xfactors>1) {for (i in 2:num.xfactors) { #first do the correlations within the f set
      for (j in 1:(i-1)) {
                 {if((!is.numeric(Phi[i,j] ) && ((Phi[i,j] !="0")||(Phi[j,i] !="0")))||  ((is.numeric(Phi[i,j]) && abs(Phi[i,j]) > cut ))) {
                           
                          if(Phi[i,j] == Phi[j,i] ) {
                                if(lr) {dia.curve(from=fact.rect[[i]]$right,to=fact.rect[[j]]$right, labels = Phi[i,j],scale=2*(i-j)/num.xfactors)} else {
                                        dia.curve(from=fact.rect[[i]]$top,to=fact.rect[[j]]$top, labels = Phi[i,j],scale=2*(i-j)/num.xfactors)}
                                                     	sem[k,1]  <- paste(fact[i],"<->",fact[j],sep="")
                                                     	sem[k,2] <-  paste("rF",i,"F",j,sep="")
                                                     	 lavaan[[num.xfactors +num.yfactors +1]] <- paste(fact[i], "~~", fact[j])} else {#directed arrows
                                                     	
                                                     	if(Phi[i,j] !="0") { if(lr) {	if(abs(i-j) < 2) {dia.arrow(from=fact.rect[[j]],to=fact.rect[[i]], labels = Phi[i,j],scale=2*(i-j)/num.xfactors)} else {
                                                     	                                                 dia.curved.arrow(from=fact.rect[[j]]$right,to=fact.rect[[i]]$right, labels = Phi[i,j],scale=2*(i-j)/num.xfactors)}
                                                     	                             } else {
                                                                                       if(abs(i-j) < 2) { dia.arrow(from=fact.rect[[j]],to=fact.rect[[i]], labels = Phi[i,j],scale=2*(i-j)/num.xfactors)} else {
                                                                                                          dia.curved.arrow(from=fact.rect[[j]]$top,to=fact.rect[[i]]$top, labels = Phi[i,j],scale=2*(i-j)/num.xfactors)}
                                                                                       }	
                                                      		sem[k,1]  <- paste(fact[j]," ->",fact[i],sep="")
                                                      		sem[k,2] <-  paste("rF",j,"F",i,sep="")} else {
                                                      		
                                                      		 if(lr) {	if(abs(i-j) < 2) {dia.arrow(from=fact.rect[[i]],to=fact.rect[[j]], labels = Phi[j,i],scale=2*(i-j)/num.xfactors)} else {
                                                     	                                                 dia.curved.arrow(from=fact.rect[[i]]$right,to=fact.rect[[j]]$right, labels = Phi[j,i],scale=2*(i-j)/num.xfactors)}
                                                     	                             } else {
                                                                                       if(abs(i-j) < 2) { dia.arrow(from=fact.rect[[i]],to=fact.rect[[j]], labels = Phi[j,i],scale=2*(i-j)/num.xfactors)} else {
                                                                                                          dia.curved.arrow(from=fact.rect[[i]]$top,to=fact.rect[[j]]$top, labels = Phi[j,i],scale=2*(i-j)/num.xfactors)}
                                                                                       }	
                                                      		
                                                      
                                                      sem[k,1]  <- paste(fact[i],"<-",fact[j],sep="")
                                                      lavaan[[num.xfactors +num.yfactors +k]] <- paste(fact[j], "~", fact[i])
                                                      sem[k,2] <-  paste("rF",i,"F",j,sep="")} 
                                                     	} } else { 
                           
                            sem[k,1]  <- paste(fact[i],"<->",fact[j],sep="")
                            
                             if (is.numeric(Phi[i,j])) {sem[k,2] <-  paste("rF",i,"F",j,sep="")} else {sem[k,2] <- Phi[i,j] } }
                            
                                   
                            k <- k + 1} }
                            
                            }
                            }  #end of correlations within the fx set
      if(!is.null(ymodel)) {
                  for (i in 1:num.xfactors) { 
                      for (j in 1:num.yfactors) {
                      		if((!is.numeric(Phi[j+num.xfactors,i] ) && (Phi[j+num.xfactors,i] !="0"))||  ((is.numeric(Phi[j+num.xfactors,i]) && abs(Phi[j+num.xfactors,i]) > cut ))) {
                      
                       dia.arrow(from=fact.rect[[i]],to=fact.rect[[j+num.xfactors]],Phi[j+num.xfactors,i],adj=i %% adj +1)
                       sem[k,1]  <- paste(fact[i],"->",fact[j+num.xfactors],sep="")
                       lavaan[[num.xfactors +num.yfactors +k]] <- paste(fact[j+num.xfactors], "~", fact[i]) } else {
                      
                       sem[k,1]  <- paste(fact[i],"<->",fact[j+num.xfactors],sep="")}
                       
                        if (is.numeric(Phi[j+num.xfactors,i])) {sem[k,2] <-  paste("rX",i,"Y",j,sep="")} else {sem[k,2] <- Phi[j+num.xfactors,i] } 
                       k <- k + 1 }
                        } }
               } }  
     if(num.factors > 0 ) { 
          
     for(f in 1:num.factors) {
       	 	sem[k,1]  <- paste(fact[f],"<->",fact[f],sep="")
        	sem[k,3] <-  "1"
        	 k <- k+1
     		} 
model=sem[1:(k-1),] 
class(model) <- "mod"   #suggested by John Fox to make the output cleaner
lavaan <- unlist(lavaan)
lavaan <- noquote(lavaan)
result <- list(sem=model,lavaan=lavaan) 
invisible(result) }
   }
   
 
