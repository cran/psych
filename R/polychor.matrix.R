"polychor.matrix" <-
function(x,y=NULL) {
		if (!require(polycor)) {stop("I am sorry, you need to have loaded the polychor package")}
 sizex <- dim(x)[2]
if (((is.data.frame(y))|(is.matrix(y))))  sizey<-dim(y)[2] 
   else  sizey <- dim(x)[2]
        
 result<-matrix(1,nrow=sizey,ncol=sizex)   #create the output array

 xnames<- names(x)
 colnames(result)<- names(x)

if (((is.data.frame(y))|(is.matrix(y))))  rownames(result) <- names(y)
    else  rownames(result) <- names(x)
 
    
if (!((is.data.frame(y))|(is.matrix(y)))) {     #default case returns a square matrix
   for (i in 2: sizex ) {
     for (j in 1:( i-1)) {
      result[j,i]<-polychor(table(x[,j],x[,i]) )
       result[i,j] <- result[j,i]
       }
      }
      }
     else {                   #if y is input, then return the rectangular array
        for (i in 1: sizex ) {
            for (j in 1:sizey) {
                 result[j,i]<-polychor(table(x[,i],y[,j]) )
                     }
                } }
   return (result) }

