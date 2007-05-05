"multi.hist" <-
function(x) {nvar <- dim(x)[2]  #number of variables
     nsize=trunc(sqrt(nvar))+1   #size of graphic
     old.par <- par(no.readonly = TRUE) # all par settings which can be changed
     par(mfrow=c(nsize,nsize))       #set new graphic parameters
     for (i in 1:nvar) {
     name=names(x)[i]                #get the names for the variables
     hist(x[,i],main=name,xlab=name) }  #draw the histograms for each variable
     on.exit(par(old.par))   #set the graphic parameters back to the original
     }

